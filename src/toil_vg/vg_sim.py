#!/usr/bin/env python
"""
vg_sim.py: this wrapper to run vg sim in parallel

"""

import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, struct, socket, threading
import string
import getpass
import pdb
import gzip
import logging

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_mapeval import run_gam_to_fastq

logger = logging.getLogger(__name__)

def sim_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("xg_indexes", nargs='+', type=make_url,
                        help="Path(s) to xg index(es) (separated by space)")
    parser.add_argument("num_reads", type=int,
                        help="Number of reads to simulate")
    parser.add_argument("out_store",
                        help="output store. All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--gam", action="store_true",
                        help="Output GAM file, annotated gam file and truth positions")
    parser.add_argument("--fastq_out", action="store_true",
                        help="Ouput fastq file (in addition to GAM)")
    parser.add_argument("--annotate_xg", type=make_url,
                        help="xg index used for gam annotation (if different from input indexes)")
    parser.add_argument("--gbwt_index", type=make_url,
                        help="GBWT index file containing individual to simulate from")
    parser.add_argument("--sample_name", type=str, default=None,
                        help="Sample name of the individual in the GBWT to simulate from")
    parser.add_argument("--drop_contigs_matching", default=[], action='append',
                        help="drop reads with true positions on contigs (e.g. decoys) matching the given regex(es) from GAM")
    parser.add_argument("--tag_bed", default=[], action='append', type=make_url,
                        help="tag reads overlapping features in the given BED(s) with the feature names")
    parser.add_argument("--path", default=[], action='append',
                        help="simulate reads from the given path name in each XG file")
    parser.add_argument("--sim_opts", type=str,
                        help="arguments for vg sim (wrapped in \"\"). Do not include -x, -n, -s, -a, or -P")
    parser.add_argument("--sim_chunks", type=int, default=1,
                        help="split simulation into this many chunks, run in parallel when possible")
    parser.add_argument("--seed", type=int, default=None,
                        help="random seed")
    parser.add_argument("--fastq", type=make_url,
                        help="use error profile derived from given reads (ignores -l, -N, -f from sim_opts)")
    parser.add_argument("--out_name", type=str,
                        help="prefix for output file names")
    parser.add_argument("--validate", action="store_true",
                        help="run vg validate on simulated GAM")
                        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def validate_sim_options(options):
    require(not options.sim_opts or all([i not in options.sim_opts for i in ['-x', '--xg-name',
        '-n', '--num-reads', '-a', '--align-out', '-s', '--random-seed', '-P', '--path',
        '-m', '--sample-name', '-g', '--gbwt-name']]),
        'sim-opts cannot contain -x, -n, -s, -a, -P, -m, or -g')
    require(options.sim_chunks > 0, '--sim_chunks must be >= 1')
    require(options.seed is None or options.seed > 0,
            'random seed must be greater than 0 (vg sim ignores seed 0)')
    require(options.gam or len(options.drop_contigs_matching) == 0,
            'can only drop reads annotated as from particular contigs if producing GAM')
    require(options.gam or not options.validate,
            '--validate on applicable to --gam output')
    require((options.gbwt_index is None) == (options.sample_name is None),
            '--gbwt_index and --sample_name must be used together')
    require(len(options.xg_indexes) == 1 or not options.gbwt_index,
            'only a single XG can be used with a --gbwt_index')
    require(not options.annotate_xg or not options.gbwt_index,
            'a --gbwt_index cannot be used with a separate --annotate_xg')
    
def run_sim(job, context, num_reads, gam, fastq_out, seed, sim_chunks,
            xg_file_ids, xg_annot_file_id, gbwt_file_id = None, sample_name = None,
            tag_bed_ids = [], paths = [], drop_contigs_matching = [],
            fastq_id = None, out_name = None, validate = False):
    """  
    run a bunch of simulation child jobs, merge up their output as a follow on
    """
    sim_out_id_infos = []

    # no seed specified, we choose one at random
    if seed is None:
        seed = random.randint(0, 2147483647)
        RealtimeLogger.info('No seed specifed, choosing random value = {}'.format(seed))

    # encapsulate follow-on
    child_job = Job()
    job.addChild(child_job)
        
    # we can have more than one xg file if we've split our input graphs up
    # into haplotypes
    for xg_i, xg_file_id in enumerate(xg_file_ids):
        file_reads = num_reads / len(xg_file_ids)
        if xg_file_id == xg_file_ids[-1]:
            file_reads += num_reads % len(xg_file_ids)
                    
        # Define a seed base for this set of chunks, leaving space for each chunk before the next seed base
        seed_base = seed + xg_i * sim_chunks 

        # each element is either reads_chunk_id or (gam_chunk_id, true_pos_chunk_id)
        # if --gam not specified
        for chunk_i in range(sim_chunks):
            chunk_reads = file_reads / sim_chunks
            if chunk_i == sim_chunks - 1:
                chunk_reads += file_reads % sim_chunks
            sim_out_id_info = child_job.addChildJobFn(run_sim_chunk, context, gam, seed_base, xg_file_id,
                                                      xg_annot_file_id, chunk_reads, chunk_i, xg_i,
                                                      gbwt_file_id=gbwt_file_id, sample_name=sample_name,
                                                      tag_bed_ids=tag_bed_ids, paths=paths,
                                                      drop_contigs_matching=drop_contigs_matching,
                                                      fastq_id=fastq_id, validate=validate,
                                                      cores=context.config.sim_cores, memory=context.config.sim_mem,
                                                      disk=context.config.sim_disk).rv()
            sim_out_id_infos.append(sim_out_id_info)
            
    merge_job = child_job.addFollowOnJobFn(run_merge_sim_chunks, context, gam,
                                           sim_out_id_infos, out_name,
                                           cores=context.config.sim_cores,
                                           memory=context.config.sim_mem,
                                           disk=context.config.sim_disk)
    
    merged_gam_id, true_id = merge_job.rv(0), merge_job.rv(1)

    if fastq_out:
        fastq_job = merge_job.addFollowOnJobFn(run_gam_to_fastq, context, merged_gam_id, False,
                                               out_name = out_name if out_name else 'sim',
                                               out_store = True,
                                               cores=context.config.sim_cores,
                                               memory=context.config.sim_mem,
                                               disk=context.config.sim_disk)
        merged_fq_id = fastq_job.rv(0)

    return merged_gam_id, true_id


def run_sim_chunk(job, context, gam, seed_base, xg_file_id, xg_annot_file_id, num_reads, chunk_i, xg_i,
                  gbwt_file_id = None, sample_name = None,
                  tag_bed_ids = [], paths = [], drop_contigs_matching = [], fastq_id = None, validate = False):
    """
    simulate some reads (and optionally gam)
    determine with true positions in xg_annot_file_id, into a true position TSV file.
    xg_annot_file_id may be none in which case xg_file_id is used.
    tag reads by overlap with tag_bed_ids BED files in the true position file
    if producing gam, drop reads whose annotated truth contigs match any of the regular expressions in drop_contigs_matching
    return either reads_chunk_id or (gam_chunk_id, tsv_chunk_id) if gam is true 
    """

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # read the xg file
    xg_file = os.path.join(work_dir, 'index.xg')
    job.fileStore.readGlobalFile(xg_file_id, xg_file)

    # read the fastq file for quality training
    if fastq_id:
        fastq_file = os.path.join(work_dir, 'error_template.fastq')
        job.fileStore.readGlobalFile(fastq_id, fastq_file)
    else:
        fastq_file = None
    
    # and the annotation xg file
    if xg_annot_file_id and xg_annot_file_id != xg_file_id:
        xg_annot_file = os.path.join(work_dir, 'annot_index.xg')
        job.fileStore.readGlobalFile(xg_annot_file_id, xg_annot_file)
    else:
        # We always need to do annotation, because the path position
        # annotations that come from vg sim are only initial path positions. vg
        # annotate actually searches a distance.
        xg_annot_file = xg_file

    # And the GCSA haplotypes
    if gbwt_file_id:
        gbwt_file = os.path.join(work_dir, 'haplotypes.gbwt')
        job.fileStore.readGlobalFile(gbwt_file_id, gbwt_file)
    else:
        gbwt_file = None

    # and the tag region BED files
    tag_beds = []
    for i, bed_id in enumerate(tag_bed_ids):
        bed_file = os.path.join(work_dir, 'tag_{}.bed'.format(i))
        job.fileStore.readGlobalFile(bed_id, bed_file)
        tag_beds.append(bed_file)

    # run vg sim
    sim_cmd = ['vg', 'sim', '-x', os.path.basename(xg_file), '-n', int(num_reads)] + context.config.sim_opts
    if seed_base is not None:
        sim_cmd += ['-s', seed_base + chunk_i]
    if fastq_id:
        sim_cmd += ['-F', os.path.basename(fastq_file)]
    if gbwt_file:
        sim_cmd += ['-g', os.path.basename(gbwt_file)]
    if sample_name:
        sim_cmd += ['-m', sample_name]
    if paths:
        for path in paths:
            # Restrict to just this path
            sim_cmd += ['-P', path]

    if not gam:
        # output reads
        reads_file = os.path.join(work_dir, 'sim_reads_{}_{}'.format(xg_i, chunk_i))

        # run vg sim
        with open(reads_file, 'wb') as output_reads:
            try:
                context.runner.call(job, sim_cmd, work_dir = work_dir, outfile=output_reads)
            except:
                # Dump everything we need to replicate the problem
                context.write_output_file(job, xg_file)
                raise

        # write to the store
        return context.write_intermediate_file(job, reads_file)
    else:
        # output gam
        gam_file = os.path.join(work_dir, 'sim_{}_{}.gam'.format(xg_i, chunk_i))
        gam_json = os.path.join(work_dir, 'sim_{}_{}.json'.format(xg_i, chunk_i))
        
        # run vg sim, write output gam, annotated gam if needed, and tagged gam json
        # (from vg/scripts/map-sim)
        cmd = [sim_cmd + ['-a']]
        
        # Current versions of vg sim annotate GAM output, but not well enough.
        # They only include initial path positions; they don't do a search like
        # vg annotate does. So we always need to run vg annotate.
        
        # We need to annotate with a different XG, or we have tags to
        # apply, or we're simulating with a GBWT and would get thread name
        # annotations from vg sim instead of the primary path.
        
        # Keep around unannotated and unfiltered GAM 
        unannotated_gam_file = os.path.join(work_dir, 'sim_{}_{}_unannotated.gam'.format(xg_i, chunk_i))
        cmd.append(['tee', os.path.basename(unannotated_gam_file)])
        
        annotate_command = ['vg', 'annotate', '-p', '-x', os.path.basename(xg_annot_file), '-a', '-']
        for tag_bed in tag_beds:
            # Make sure to annotate with overlap to tag regions
            annotate_command += ['-b', os.path.basename(tag_bed)]
        cmd.append(annotate_command)
        
        if len(drop_contigs_matching) > 0:
            # We need to filter out reads annotated as being on these contigs
            filter_command = ['vg', 'filter', '-']
            for exclude_regex in drop_contigs_matching:
                filter_command.append('-X')
                filter_command.append(exclude_regex)
            cmd.append(filter_command)
                
        cmd.append(['tee', os.path.basename(gam_file)])
        cmd.append(['vg', 'view', '-aj', '-'])
        with open(gam_json, 'wb') as output_json:
            try:
                context.runner.call(job, cmd, work_dir = work_dir, outfile=output_json)
                if validate:
                    context.runner.call(job, ['vg', 'validate',
                                              '--gam', os.path.basename(gam_file),
                                              os.path.basename(xg_file)], work_dir = work_dir)
            except:
                # Dump everything we need to replicate the problem
                context.write_output_file(job, xg_file)
                if xg_annot_file != xg_file:
                    context.write_output_file(job, xg_annot_file)
                context.write_output_file(job, unannotated_gam_file)
                context.write_output_file(job, gam_file)
                for tag_bed in tag_beds:
                    context.write_output_file(job, tag_bed)
                if fastq_file is not None:
                    context.write_output_file(job, fastq_file)
                raise

        # turn the annotated gam json into truth positions, as separate command since
        # we're going to use a different docker container.  (Note, would be nice to
        # avoid writing the json to disk)
        # note: in the following, we are writing the read name as the first column,
        # then a comma-separated tag list of overlapped features (or ".") in the second column,
        # then a path-name, path-offset in each successive pair of columns, and finally
        # a score and MAPQ so we match the extract_gam_read_stats/extract_bam_read_stats
        # format from toil-vg mapeval.
        jq_cmd = ['jq', '-c', '-r', '[.name] + '
                  'if (.annotation.features | length) > 0 then [.annotation.features | join(",")] else ["."] end + '
                  'if .refpos != null then [.refpos[] | .name, if .offset != null then .offset else 0 end] else [] end + '
                  '[.score] + '
                  'if .mapping_quality == null then [0] else [.mapping_quality] end | @tsv',
                   os.path.basename(gam_json)]
        # output truth positions
        true_pos_file = os.path.join(work_dir, 'true_{}_{}.pos'.format(xg_i, chunk_i))
        with open(true_pos_file, 'wb') as out_true_pos:
            context.runner.call(job, jq_cmd, work_dir = work_dir, outfile=out_true_pos)

        # get rid of that big json asap
        os.remove(gam_json)

        # write to file store, but don't dump all the chunks in the output store
        gam_chunk_id = context.write_intermediate_file(job, gam_file)
        true_pos_chunk_id = context.write_intermediate_file(job, true_pos_file)

        # return everything as a tuple.
        return gam_chunk_id, true_pos_chunk_id
        

def run_merge_sim_chunks(job, context, gam, sim_out_id_infos, out_name):
    """
    merge the sim output

    Takes a GAM flag for if we simulated a GAM, a list of sim chunk output
    values (which are single file IDs if not making GAM, or GAM file and true
    position TSV file ID pairs if making GAM), and an output file name part.

    Returns a merged file ID if not making GAMs, or a pair of merged GAM and truth TSV files if making GAM.

    """
    assert len(sim_out_id_infos) > 0

    work_dir = job.fileStore.getLocalTempDir()

    if out_name:
        reads_name = out_name
        pos_name = out_name
    else:
        reads_name = 'sim'
        pos_name = 'true'

    if not gam:
        # merge up the reads files
        merged_reads_file = os.path.join('{}.reads'.format(reads_name))
        with open(merged_reads_file, 'a') as out_reads:
            for i, reads_file_id in enumerate(sim_out_id_infos):
                reads_file = os.path.join(work_dir, 'sim_reads_{}'.format(i))
                job.fileStore.readGlobalFile(reads_file_id, reads_file)
                with open(reads_file) as rf:
                    shutil.copyfileobj(rf, out_reads)

        return context.write_output_file(job, merged_reads_file)
    
    else:
        # merge up the gam files
        merged_gam_file = os.path.join(work_dir, '{}.gam'.format(reads_name))
        merged_true_file = os.path.join(work_dir, '{}.pos.unsorted'.format(pos_name))
        
        with open(merged_gam_file, 'ab') as out_gam, \
             open(merged_true_file, 'ab') as out_true:
            
            for i, sim_out_id_info in enumerate(sim_out_id_infos):
                gam_file = os.path.join(work_dir, 'sim_{}.gam'.format(i))
                job.fileStore.readGlobalFile(sim_out_id_info[0], gam_file)
                with open(gam_file, 'rb') as rf:
                    shutil.copyfileobj(rf, out_gam)
                    
                true_file = os.path.join(work_dir, 'true_{}.pos'.format(i))
                job.fileStore.readGlobalFile(sim_out_id_info[1], true_file)
                with open(true_file, 'rb') as rf:
                    shutil.copyfileobj(rf, out_true)

        # sort the positions file
        sorted_true_file = os.path.join(work_dir, '{}.pos'.format(pos_name))
        sort_cmd = ['sort', os.path.basename(merged_true_file)]
        with open(sorted_true_file, 'wb') as out_true:
            context.runner.call(job, sort_cmd, work_dir = work_dir, outfile = out_true)

        
        merged_gam_id = context.write_output_file(job, merged_gam_file)
        true_id = context.write_output_file(job, sorted_true_file)

        return merged_gam_id, true_id
            
def sim_main(context, options):
    """
    Wrapper for vg sim. 
    """

    validate_sim_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            inputXGFileIDs = []
            for xg_index in options.xg_indexes:
                inputXGFileIDs.append(importer.load(xg_index))
            if options.annotate_xg:
                inputAnnotXGFileID = importer.load(options.annotate_xg)
            else:
                inputAnnotXGFileID = None
            if options.gbwt_index:
                inputGBWTIndexFileID = importer.load(options.gbwt_index)
            else:
                inputGBWTIndexFileID = None
            tagBEDIDs = []
            for url in options.tag_bed:
                tagBEDIDs.append(importer.load(url))
            if options.fastq:
                inputFastqFileID = importer.load(options.fastq)
            else:
                inputFastqFileID = None

            # can't make the fastq without going through gam
            if options.fastq_out:
                options.gam = True                

            importer.wait()

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            # FASTQs no longer need unzipping because vg can read .gz transparently.

            # Make a root job
            root_job = Job.wrapJobFn(run_sim, context, options.num_reads, options.gam,
                                     options.fastq_out,
                                     options.seed, options.sim_chunks,
                                     importer.resolve(inputXGFileIDs),
                                     importer.resolve(inputAnnotXGFileID),
                                     gbwt_file_id = importer.resolve(inputGBWTIndexFileID),
                                     sample_name = options.sample_name,
                                     tag_bed_ids = importer.resolve(tagBEDIDs),
                                     paths = options.path,
                                     drop_contigs_matching = options.drop_contigs_matching,
                                     fastq_id = importer.resolve(inputFastqFileID),
                                     out_name = options.out_name,
                                     validate = options.validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
                                     
            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
