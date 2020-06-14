#!/usr/bin/env python
"""
vg_map.py: map to vg graph producing gam for each chrom

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
from toil_vg.vg_surject import run_merge_bams, run_whole_surject 

logger = logging.getLogger(__name__)

def map_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("sample_name", type=str,
                        help="sample name (ex NA12878)")
    
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--kmer_size", type=int,
                        help="size of kmers to use in gcsa-kmer mapping mode")
        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add mapping options shared with mapeval and run
    map_parse_args(parser)

    # Add mapping options shared onyl with run
    map_parse_index_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)
    
def map_parse_index_args(parser):
    """
    Define map arguments shared with run but not mapeval
    """
    
    parser.add_argument("--xg_index", type=make_url,
                        help="Path to xg index")    
    parser.add_argument("--gcsa_index", type=make_url,
                        help="Path to GCSA index (for map and mpmap)")
    parser.add_argument("--minimizer_index", type=make_url,
                        help="Path to minimizer index (for giraffe)")
    parser.add_argument("--distance_index", type=make_url,
                        help="Path to distance index (for giraffe)")
    parser.add_argument("--gbwt_index", type=make_url,
                        help="Path to GBWT haplotype index")
    parser.add_argument("--snarls_index", type=make_url,
                        help="Path to snarls file")
                        
    parser.add_argument("--mapper", default="map", choices=["map", "mpmap", "giraffe"],
                        help="vg mapper to use")

def map_parse_args(parser, stand_alone = False):
    """
    Define map arguments shared with mapeval and run
    """
    
    parser.add_argument("--fastq", nargs='+', type=make_url,
                        help="Input fastq (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fq_split_cores", type=int,
                        help="number of threads used to split input FASTQs")   
    parser.add_argument("--gam_input_reads", type=make_url, default=None,
                        help="Input reads in GAM format")
    parser.add_argument("--bam_input_reads", type=make_url, default=None,
                        help="Input reads in BAM format")
    parser.add_argument("--single_reads_chunk", action="store_true", default=False,
                        help="do not split reads into chunks")
    parser.add_argument("--reads_per_chunk", type=int,
                        help="number of reads for each mapping job")
    parser.add_argument("--alignment_cores", type=int,
                        help="number of threads during the alignment step")
    parser.add_argument("--interleaved", action="store_true", default=False,
                        help="treat fastq as interleaved read pairs. overrides *_opts")
    parser.add_argument("--map_opts", type=str,
                        help="arguments for vg map (wrapped in \"\")")
    parser.add_argument("--mpmap_opts", type=str,
                        help="arguments for vg mpmap (wrapped in \"\")")
    parser.add_argument("--giraffe_opts", type=str,
                        help="arguments for vg giraffe (wrapped in \"\")")
    parser.add_argument("--bam_output", action="store_true",
                        help="write BAM output directly")
    parser.add_argument("--surject", action="store_true",
                        help="surject output, producing BAM in addition to GAM alignments")
    parser.add_argument("--validate", action="store_true",
                        help="run vg validate on ouput GAMs")
    parser.add_argument("--id_ranges", type=make_url, default=None,
                        help="Path to file with node id ranges for each chromosome in BED format.")

    
def validate_map_options(context, options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(options.xg_index is not None, 'All mappers require --xg_index')
    
    if options.mapper == 'map' or options.mapper == 'mpmap':
        require(options.gcsa_index, '--gcsa_index is required for map and mpmap')
    
    if options.mapper == 'giraffe':
        require(options.minimizer_index, '--minimizer_index is required for giraffe')
        require(options.distance_index, '--distance_index is required for giraffe')
        require(options.gbwt_index, '--gbwt_index is required for giraffe')
        require(not options.bam_input_reads, '--bam_input_reads is not supported with giraffe')
        require(not options.interleaved, '--interleaved is not supported with giraffe')
        require(options.fastq is None or len(options.fastq) < 2, 'Multiple --fastq files are not supported with giraffe')
    
    require(options.fastq is None or len(options.fastq) in [1, 2], 'Exacty 1 or 2 files must be'
            ' passed with --fastq')
    require(options.interleaved == False or options.fastq is None or len(options.fastq) == 1,
            '--interleaved cannot be used when > 1 fastq given')
    require(sum([1 if x else 0 for x in [options.fastq, options.gam_input_reads, options.bam_input_reads]]) == 1,
            'reads must be speficied with either --fastq or --gam_input_reads or --bam_input_reads')
    require(options.mapper == 'mpmap' or options.snarls_index is None,
            '--snarls_index can only be used with --mapper mpmap') 
    if options.mapper == 'mpmap':
        require('-S' in context.config.mpmap_opts or '--single-path-mode' in context.config.mpmap_opts,
                '-S must be used with mpmap mapper to produce GAM output')
        require(not options.bam_output,
                '--bam_output not currently supported with mpmap mapper')
    require (not options.bam_output or not options.surject,
             '--bam_output cannot be used in combination with --surject')
    require (not options.id_ranges or not options.surject,
             '--surject not currently supported with --id_ranges')
        

def run_split_reads_if_needed(job, context, fastq, gam_input_reads, bam_input_reads, reads_file_ids):
    """
    Return a list of lists of read chunk file IDs, one list per read files.
    
    If the workflow is in single_reads_chunk mode (according to
    context.options.single_read_chunk), produce one chunk per file.
    
    Otherwise, produce several chunks per file.
    """
    
    if not context.config.single_reads_chunk:
        reads_chunk_ids = job.addChildJobFn(run_split_reads, context, fastq, gam_input_reads, bam_input_reads,
                                            reads_file_ids,
                                            cores=context.config.misc_cores, memory=context.config.misc_mem,
                                            disk=context.config.misc_disk).rv()
    else:
        RealtimeLogger.info("Bypassing reads splitting because --single_reads_chunk enabled")
        reads_chunk_ids = [[r] for r in reads_file_ids]
        
    return reads_chunk_ids
    

def run_mapping(job, context, fastq, gam_input_reads, bam_input_reads, sample_name, interleaved, mapper,
                indexes, reads_file_ids=None, reads_chunk_ids=None,
                bam_output=False, surject=False, 
                gbwt_penalty=None, validate=False):
    """
    Split the fastq, then align each chunk.
    
    Exactly one of fastq, gam_input_reads, or bam_input_reads should be
    non-falsey, to indicate what kind of data the file IDs in reads_file_ids or
    reads_chunk_ids correspond to.
    
    Exactly one of reads_file_ids or read_chunks_ids should be specified.
    reads_file_ids holds a list of file IDs of non-chunked input read files,
    which will be chunked if necessary. reads_chunk_ids holds lists of chunk
    IDs for each read file, as produced by run_split_reads_if_needed.
    
    indexes is a dict from index type ('xg', 'gcsa', 'lcp', 'id_ranges',
    'gbwt', 'minimizer', 'distance', 'snarls') to index file ID. Some indexes
    are extra and specifying them will change mapping behavior. Some indexes
    are required for certain values of mapper.
    
    mapper can be 'map', 'mpmap', or 'giraffe'. For 'map' and 'mpmap', the 'gcsa'
    and 'lcp' indexes are required. For 'giraffe', the 'gbwt', 'minimizer' and
    'distance' indexes are required. All the mappers require the 'xg' index.
    
    If bam_output is set, produce BAMs. If surject is set, surject reads down
    to paths. 
    
    If the 'gbwt' index is present and gbwt_penalty is specified, the default
    recombination penalty will be overridden.
    
    returns output gams, one per chromosome, the total mapping time (excluding
    toil-vg overhead such as transferring and splitting files), and output
    BAMs, one per chromosome, if computed.
    """
    
    # Make sure we have exactly one type of input
    assert (bool(fastq) + bool(gam_input_reads) + bool(bam_input_reads) == 1)
    
    # Make sure we have exactly one kind of file IDs
    assert(bool(reads_file_ids) + bool(reads_chunk_ids) == 1)

    # We may have to have a job to chunk the reads
    chunk_job = None

    if reads_chunk_ids is None:
        # If the reads are not pre-chunked for us, we have to chunk them.
        chunk_job = job.addChildJobFn(run_split_reads_if_needed, context, fastq, gam_input_reads, bam_input_reads,
                                      reads_file_ids, cores=context.config.misc_cores, memory=context.config.misc_mem,
                                      disk=context.config.misc_disk)
        reads_chunk_ids = chunk_job.rv()
        
    # We need a job to do the alignment
    align_job = Job.wrapJobFn(run_whole_alignment, context, fastq, gam_input_reads, bam_input_reads, sample_name,
                              interleaved, mapper, indexes, reads_chunk_ids,
                              bam_output=bam_output, surject=surject,
                              gbwt_penalty=gbwt_penalty,
                              validate=validate,
                              cores=context.config.misc_cores,
                              memory=context.config.misc_mem, disk=context.config.misc_disk)
                 
    if chunk_job is not None:
        # Alignment must happen after chunking
        chunk_job.addFollowOn(align_job)
    else:
        # Alignment can happen now
        job.addChild(align_job)
                 
    return align_job.rv()

def run_split_reads(job, context, fastq, gam_input_reads, bam_input_reads, reads_file_ids):
    """
    split either fastq or gam input reads into chunks.  returns list of chunk file id lists
    (one for each input reads file)
    """

    # this is a list of lists: one list of chunks for each input reads file
    reads_chunk_ids = []
    if fastq and len(fastq):
        for fastq_i, reads_file_id in enumerate(reads_file_ids):
            reads_chunk_ids.append(job.addChildJobFn(run_split_fastq, context, fastq, fastq_i, reads_file_id,
                                                     cores=context.config.fq_split_cores,
                                                     memory=context.config.fq_split_mem, disk=context.config.fq_split_disk).rv())
    elif gam_input_reads:
        assert len(reads_file_ids) == 1
        reads_chunk_ids.append(job.addChildJobFn(run_split_gam_reads, context, gam_input_reads, reads_file_ids[0],
                                                  cores=context.config.fq_split_cores,
                                                  memory=context.config.fq_split_mem, disk=context.config.fq_split_disk).rv())
    else:
        assert bam_input_reads
        assert len(reads_file_ids) == 1
        reads_chunk_ids.append(job.addChildJobFn(run_split_bam_reads, context, bam_input_reads, reads_file_ids[0],
                                                  cores=context.config.fq_split_cores,
                                                  memory=context.config.fq_split_mem, disk=context.config.fq_split_disk).rv())
        
    return reads_chunk_ids


def run_split_fastq(job, context, fastq, fastq_i, sample_fastq_id):
   
    disk_required = job.fileStore.getGlobalFileSize(sample_fastq_id) * 2 + (2 * 1024**3)
    requeued = ensure_disk_bytes(job, run_split_fastq, disk_required)
    if requeued is not None:
        # If not, requeue the job with more disk.
        return requeued
   
    RealtimeLogger.info("Starting fastq split")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    fastq_name = os.path.basename(fastq[fastq_i])
    fastq_path = os.path.join(work_dir, fastq_name)
    fastq_gzipped = os.path.splitext(fastq_name)[1] == '.gz'
    fastq_name = os.path.splitext(fastq_name)[0]
    if fastq_gzipped:
        fastq_name = os.path.splitext(fastq_name)[0]
    job.fileStore.readGlobalFile(sample_fastq_id, fastq_path)

    # Split up the fastq into chunks

    # Make sure chunk size even in case paired interleaved
    chunk_size = context.config.reads_per_chunk
    if chunk_size % 2 != 0:
        chunk_size += 1

    # 4 lines per read
    chunk_lines = chunk_size * 4

    # Note we do this on the command line because Python is too slow
    if fastq_gzipped:
        cmd = [['gzip', '-d', '-c', os.path.basename(fastq_path)]]
    else:
        cmd = [['cat', os.path.basename(fastq_path)]]

    cmd.append(['split', '-l', str(chunk_lines),
                '--filter=pigz -p {} > $FILE.fq.gz'.format(max(1, int(context.config.fq_split_cores) - 1)),
                '-', '{}-chunk.'.format(fastq_name)])

    context.runner.call(job, cmd, work_dir = work_dir, tool_name='pigz')

    fastq_chunk_ids = []
    for chunk_name in sorted(os.listdir(work_dir)):
        if chunk_name.endswith('.fq.gz') and chunk_name.startswith('{}-chunk'.format(fastq_name)):
            fastq_chunk_ids.append(context.write_intermediate_file(job, os.path.join(work_dir, chunk_name)))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Split fastq into {} chunks. Process took {} seconds.".format(len(fastq_chunk_ids), run_time))

    return fastq_chunk_ids

def run_split_gam_reads(job, context, gam_input_reads, gam_reads_file_id):
    """ split up an input reads file in GAM format
    """
    
    disk_required = job.fileStore.getGlobalFileSize(gam_reads_file_id) * 2 + (2 * 1024**3)
    requeued = ensure_disk_bytes(job, run_split_gam_reads, disk_required)
    if requeued is not None:
        # If not, requeue the job with more disk.
        return requeued
    
    RealtimeLogger.info("Starting gam split")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    gam_path = os.path.join(work_dir, os.path.basename(gam_input_reads))
    job.fileStore.readGlobalFile(gam_reads_file_id, gam_path)

    # Split up the gam into chunks

    # Make sure chunk size even in case paired interleaved
    chunk_size = context.config.reads_per_chunk
    if chunk_size % 2 != 0:
        chunk_size += 1

    cmd = ['vg', 'chunk', '-a', os.path.basename(gam_path), '--gam-split-size', str(chunk_size),
           '--prefix', 'gam_reads_chunk']

    context.runner.call(job, cmd, work_dir = work_dir)

    gam_chunk_ids = []
    for chunk_name in os.listdir(work_dir):
        if chunk_name.endswith('.gam') and chunk_name.startswith('gam_reads_chunk'):
            gam_chunk_ids.append(context.write_intermediate_file(job, os.path.join(work_dir, chunk_name)))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Split gam into {} chunks. Process took {} seconds.".format(len(gam_chunk_ids), run_time))

    return gam_chunk_ids

def run_split_bam_reads(job, context, bam_input_reads, bam_reads_file_id):
    """ split up an input reads file in BAM format
    """
    
    disk_required = job.fileStore.getGlobalFileSize(bam_reads_file_id) * 2 + (2 * 1024**3)
    requeued = ensure_disk_bytes(job, run_split_bam_reads, disk_required)
    if requeued is not None:
        # If not, requeue the job with more disk.
        return requeued
    
    RealtimeLogger.info("Starting bam split")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    bam_path = os.path.join(work_dir, os.path.basename(bam_input_reads))
    job.fileStore.readGlobalFile(bam_reads_file_id, bam_path)

    # Split up the bam into chunks

    # Make sure chunk size even in case paired interleaved
    chunk_size = context.config.reads_per_chunk
    if chunk_size % 2 != 0:
        chunk_size += 1

    # 1 line per read
    chunk_lines = chunk_size * 1

    cmd = [['samtools', 'view', os.path.basename(bam_path)]]
    cmd.append(['split', '-l', str(chunk_lines),
                '--filter=bash -c \'cat <(samtools view -H {}) <(cat -)'.format(os.path.basename(bam_path)) +
                ' | samtools view -O BAM --threads {} -'.format(max(1, int(context.config.fq_split_cores) - 1)) +
                ' > $FILE.bam\'', '-', 'bam_reads_chunk.'])

    context.runner.call(job, cmd, work_dir = work_dir)

    bam_chunk_ids = []
    for chunk_name in sorted(os.listdir(work_dir)):
        if chunk_name.endswith('.bam') and chunk_name.startswith('bam_reads_chunk'):
            bam_chunk_ids.append(context.write_intermediate_file(job, os.path.join(work_dir, chunk_name)))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Split bam into {} chunks. Process took {} seconds.".format(len(bam_chunk_ids), run_time))

    return bam_chunk_ids

    
def run_whole_alignment(job, context, fastq, gam_input_reads, bam_input_reads, sample_name, interleaved, mapper,
                        indexes, reads_chunk_ids,
                        bam_output=False, surject=False, gbwt_penalty=None, validate=False):
    """
    align all fastq chunks in parallel
    
    Takes a dict from index type to index file ID. Some indexes are extra and
    specifying them will change mapping behavior.
    
    Returns a list of per-contig GAMs, the total allignment runtime, and a list
    of per-contig BAM file IDs (which is only nonempty when surject is true).
    
    """
    
    # this will be a list of lists.
    # gam_chunk_file_ids[i][j], will correspond to the jth path (from id_ranges)
    # for the ith gam chunk (generated from fastq shard i)
    gam_chunk_file_ids = []
    gam_chunk_running_times = []
    # depending on bam_output and surject options, we can make bam_output too
    bam_chunk_file_ids = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)

    for chunk_id, chunk_filename_ids in enumerate(zip(*reads_chunk_ids)):
        #Run graph alignment on each fastq chunk
        chunk_alignment_job = child_job.addChildJobFn(run_chunk_alignment, context, gam_input_reads, bam_input_reads,
                                                      sample_name,
                                                      interleaved, mapper, chunk_filename_ids, chunk_id,
                                                      indexes,
                                                      bam_output=bam_output,
                                                      gbwt_penalty=gbwt_penalty,
                                                      validate=validate,
                                                      cores=context.config.alignment_cores, memory=context.config.alignment_mem,
                                                      disk=context.config.alignment_disk)
        if not bam_output:
            gam_chunk_file_ids.append(chunk_alignment_job.rv(0))
        else:
            bam_chunk_file_ids.append(chunk_alignment_job.rv(0))
        gam_chunk_running_times.append(chunk_alignment_job.rv(1))


    if not bam_output:
        merge_gams_job = child_job.addFollowOnJobFn(run_merge_gams, context, sample_name, indexes.get('id_ranges'), gam_chunk_file_ids,
                                                    gam_chunk_running_times,
                                                    cores=context.config.misc_cores,
                                                    memory=context.config.misc_mem, disk=context.config.misc_disk)
        gam_chrom_ids = merge_gams_job.rv(0)
        gam_chunk_time = merge_gams_job.rv(1)
        bam_chrom_ids = []
    else:
        gam_chrom_ids = []
        gam_chunk_time = None
        merge_bams_job = child_job.addFollowOnJobFn(run_merge_bams, context, sample_name, bam_chunk_file_ids)
        bam_chrom_ids = [merge_bams_job.rv()]

    if surject:
        interleaved_surject = interleaved or (fastq and len(fastq) == 2)
        zip_job = child_job.addFollowOnJobFn(run_zip_surject_input, context, gam_chunk_file_ids)
        xg_id = indexes['xg-surject'] if 'xg-surject' in indexes else indexes['xg']
        bam_chrom_ids = [zip_job.addFollowOnJobFn(run_whole_surject, context, zip_job.rv(), sample_name + '-surject',
                                                  interleaved_surject, xg_id, []).rv()]

    return gam_chrom_ids, gam_chunk_time, bam_chrom_ids
    
def run_zip_surject_input(job, context, gam_chunk_file_ids):
    """
    run_whole_surject takes input in different format than what we have above, so we shuffle the 
    promised lists around here to avoid a (probably-needed) refactor of the existing interface
    """
    return list(zip(*gam_chunk_file_ids))


def run_chunk_alignment(job, context, gam_input_reads, bam_input_reads, sample_name, interleaved, mapper,
                        chunk_filename_ids, chunk_id, indexes,
                        bam_output=False, gbwt_penalty=None, always_check_population=True, validate=False):
                        
    """
    Align a chunk of reads.
    
    Takes a dict from index type to index file ID. Some indexes are extra and
    specifying them will change mapping behavior.
    """
    
    # TODO: Work out what indexes we will want and ensure_disk_bytes

    RealtimeLogger.info("Starting {} alignment on {} chunk {}".format(mapper, sample_name, chunk_id))

    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_file = os.path.join(work_dir, "graph.vg")

    # Work out what index files we need
    index_files = {}
    index_files['xg'] = graph_file + ".xg"
    if mapper == 'map' or mapper == 'mpmap':
        index_files['gcsa'] = graph_file + ".gcsa"
        index_files['lcp'] = index_files['gcsa'] + ".lcp"
        
        if 'gbwt' in indexes:
            # We have a GBWT haplotype index available.
            index_files['gbwt'] = graph_file + ".gbwt"
            
    if mapper == 'mpmap':
        if 'snarls' in indexes:
            # mpmap knows how to use the snarls, and we have them, so we should use them
            
            # Note that passing them will affect mapping, if using multiple
            # tracebacks. Since we only run single path mode, if multiple
            # tracebacks aren't used, mpmap will ignore the snarls.
            index_files['snarls'] = graph_file + ".snarls"
        
    if mapper == 'giraffe':
        index_files['minimizer'] = graph_file + ".min"
        index_files['distance'] = graph_file + ".dist"
        index_files['gbwt'] = graph_file + ".gbwt"
     
        
    for index_type in list(index_files.keys()):
        # Download each index file
        job.fileStore.readGlobalFile(indexes[index_type], index_files[index_type])
    
    # We need the sample reads (fastq(s) or gam) for alignment
    reads_files = []
    reads_ext = 'gam' if gam_input_reads else 'bam' if bam_input_reads else 'fq.gz'
    for j, chunk_filename_id in enumerate(chunk_filename_ids):
        reads_file = os.path.join(work_dir, 'reads_chunk_{}_{}.{}'.format(chunk_id, j, reads_ext))
        job.fileStore.readGlobalFile(chunk_filename_id, reads_file)
        reads_files.append(reads_file)
    
    # And a temp file for our aligner output
    output_file = os.path.join(work_dir, "{}_{}.gam".format(sample_name, chunk_id))

    # Open the file stream for writing
    with open(output_file, 'wb') as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = []
        
        if mapper == 'mpmap':
            vg_parts += ['vg', 'mpmap']
            vg_parts += context.config.mpmap_opts
            if '-S' not in vg_parts and '--single-path-mode' not in vg_parts:
                RealtimeLogger.warning('Adding --single-path-mode to mpmap options as only GAM output supported')
                vg_parts += ['--single-path-mode']
        elif mapper == 'map':
            vg_parts += ['vg', 'map'] 
            vg_parts += context.config.map_opts
        elif mapper == 'giraffe':
            vg_parts += ['vg', 'giraffe'] 
            vg_parts += context.config.giraffe_opts
        else:
            raise RuntimeError('Unimplemented mapper "{}"'.format(mapper))
            
        for reads_file in reads_files:
            input_flag = '-G' if gam_input_reads else '-b' if bam_input_reads else '-f'
            vg_parts += [input_flag, os.path.basename(reads_file)]
        vg_parts += ['-t', str(context.config.alignment_cores)]

        # Override the -i flag in args with the --interleaved command-line flag
        if interleaved is True and '-i' not in vg_parts and '--interleaved' not in vg_parts:
            vg_parts += ['-i']
        elif interleaved is False and 'i' in vg_parts:
            del vg_parts[vg_parts.index('-i')]
        if interleaved is False and '--interleaved' in vg_parts:
            del vg_parts[vg_parts.index('--interleaved')]

        # Override the --surject-to option
        if bam_output is True and '--surject-to' not in vg_parts:
            vg_parts += ['--surject-to', 'bam']
        elif bam_output is False and '--surject-to' in vg_parts:
            sidx = vg_parts.index('--surject-to')
            del vg_parts[sidx]
            del vg_parts[sidx]

        # Turn indexes into options
        type_to_option = {
            'gbwt': '--gbwt-name',
            'xg': '-x',
            'gcsa': '-g',
            'lcp': None,
            'distance': '-d',
            'minimizer': '-m',
            'snarls': '--snarls'
        }
        for index_type, index_file in list(index_files.items()):
            if type_to_option[index_type] is not None:
                vg_parts += [type_to_option[index_type], os.path.basename(index_file)]

        if 'gbwt' in index_files:
            # We may have a GBWT recombination rate/penalty override
            if gbwt_penalty is not None:
                # We have a recombination penalty value to apply
                if '--recombination-penalty' in vg_parts:
                    # Make sure to strip out the penalty if it is in args already
                    sidx = vg_parts.index('--recombination-penalty')
                    del vg_parts[sidx]
                    del vg_parts[sidx]
                    
                # Both map and mpmap take this option
                vg_parts += ['--recombination-penalty', str(gbwt_penalty)]
                
            if mapper == 'mpmap' and always_check_population:
                # Always try to population-score even unambiguous reads
                # mpmap can do this
                vg_parts += ['--always-check-population']

        RealtimeLogger.info(
            "Running VG for {} against {}: {}".format(sample_name, graph_file,
            " ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        command = vg_parts
        try:
            context.runner.call(job, command, work_dir = work_dir, outfile=alignment_file)
            end_time = timeit.default_timer()
            if validate:
                alignment_file.flush()
                context.runner.call(job, ['vg', 'validate', '--xg', os.path.basename(index_files['xg']),
                                          '--gam', os.path.basename(output_file)], work_dir = work_dir)
        except:
            # Dump everything we need to replicate the alignment
            end_time = timeit.default_timer()
            logging.error("Mapping failed. Dumping files.")
            for index_file in list(index_files.values()):
                context.write_output_file(job, index_file)
            for reads_file in reads_files:
                context.write_output_file(job, reads_file)
            raise
        
        # Mark when it's done
        run_time = end_time - start_time

    paired_end = '-i' in vg_parts or '--interleaved' in vg_parts or len(chunk_filename_ids) > 1
    RealtimeLogger.info("Aligned {}. Process took {} seconds with {} vg-{}".format(
        output_file, run_time, 'paired-end' if paired_end else 'single-end', mapper))

    if 'id_ranges' in indexes:
        # Break GAM into multiple chunks at the end. So we need the file
        # defining those chunks.
        id_ranges_file = os.path.join(work_dir, 'id_ranges.tsv')
        job.fileStore.readGlobalFile(indexes['id_ranges'], id_ranges_file)

        # Chunk the gam up by chromosome
        gam_chunks = split_gam_into_chroms(job, work_dir, context, index_files['xg'], id_ranges_file, output_file)

        # Write gam_chunks to store
        gam_chunk_ids = []
        for gam_chunk in gam_chunks:
            gam_chunk_ids.append(context.write_intermediate_file(job, gam_chunk))

        return gam_chunk_ids, run_time
        
    else:
        # We can just report one chunk of everything
        return [context.write_intermediate_file(job, output_file)], run_time

def split_gam_into_chroms(job, work_dir, context, xg_file, id_ranges_file, gam_file):
    """
    Create a sorted GAM index then use it to split up the given gam file
    (a local path) into a separate gam for each chromosome.  
    Return a list of filenames.  the ith filename will correspond
    to the ith path in the options.path_name list
    """

    # transfer id ranges to N:M format for new vg chunk interface
    cid_ranges_file = os.path.join(work_dir, 'cid_ranges.txt')
    with open(cid_ranges_file, 'w') as ranges, open(id_ranges_file) as in_ranges:
        for line in in_ranges:
            toks = line.strip().split()
            if len(toks) == 3:
                ranges.write('{}:{}\n'.format(toks[1], toks[2]))
 
    # Sort the GAM
    output_sorted = gam_file + '.sorted.gam'
    output_index = output_sorted + '.gai'
    sort_cmd = ['vg', 'gamsort', '-i', os.path.basename(output_index),
        '-t', str(context.config.alignment_cores), os.path.basename(gam_file)]
    with open(output_sorted, 'wb') as sorted_file:
        context.runner.call(job, sort_cmd, work_dir = work_dir, outfile = sorted_file)
 
    # Chunk the alignment into chromosomes using the id ranges
    # (note by using id ranges and 0 context and -a we avoid costly subgraph extraction)

    output_bed_name = os.path.join(work_dir, 'output_bed.bed')
    
    # Now run vg chunk on the gam index to get our gams
    # Note: using -a -c 0 -R bypasses subgraph extraction, which is important
    # as it saves a ton of time and memory
    chunk_cmd = ['vg', 'chunk', '-x', os.path.basename(xg_file),
                 '-a', os.path.basename(output_sorted), '-c', '0',
                 '-R', os.path.basename(cid_ranges_file),
                 '-b', os.path.splitext(os.path.basename(gam_file))[0],
                 '-t', str(context.config.alignment_cores),
                 '-E', os.path.basename(output_bed_name)]
    
    context.runner.call(job, chunk_cmd, work_dir = work_dir)
    
    # scrape up the vg chunk results into a list of paths to the output gam
    # chunks and return them.  we expect the filenames to be in the 4th column
    # of the output bed from vg chunk. 
    gam_chunks = []
    with open(output_bed_name) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                gam_chunks.append(os.path.join(work_dir, os.path.basename(toks[3].strip())))
                assert os.path.splitext(gam_chunks[-1])[1] == '.gam'

    return gam_chunks


def run_merge_gams(job, context, sample_name, id_ranges_file_id, gam_chunk_file_ids, gam_chunk_running_times):
    """
    Merge together gams, doing each chromosome in parallel
    
    Also totals up runtimes
    
    Returns the merged GAMs as a list, one per chromosome, and the total allignment runtime
    """
    
    if id_ranges_file_id is not None:
        # Get the real chromosome names
        id_ranges = parse_id_ranges(job, id_ranges_file_id)
        chroms = [x[0] for x in id_ranges]
    else:
        # Dump everything in a single default chromosome chunk with a default
        # name
        chroms = ["default"]
    
    chr_gam_ids = []

    for i, chr in enumerate(chroms):
        shard_ids = [gam_chunk_file_ids[j][i] for j in range(len(gam_chunk_file_ids))]
        assert(len(shard_ids) >= 1)
        chr_gam_id = job.addChildJobFn(run_merge_chrom_gam, context, sample_name, chr, shard_ids,
                                       cores=context.config.fq_split_cores,
                                       memory=context.config.fq_split_mem,
                                       disk=context.config.fq_split_disk).rv()
        chr_gam_ids.append(chr_gam_id)

    total_running_time = 0
    for running_time in gam_chunk_running_times:
        if running_time is None:
            total_running_time = None
            break
        else:
            total_running_time += float(running_time)
    
    return chr_gam_ids, total_running_time
    
def run_merge_chrom_gam(job, context, sample_name, chr_name, chunk_file_ids):
    """
    Make a chromosome gam by merging up a list of gam ids, one 
    for each  shard.
    """
    
    RealtimeLogger.info('For chrom {}, merge files: {}'.format(chr_name, chunk_file_ids))
    RealtimeLogger.info('Args: {} {} {} {} {}'.format(job, context, sample_name, chr_name, chunk_file_ids))
    
    # Check disk requirements to make sure we have enough room for 2 copies of everything, plus a bit.
    # Get file size in a way that is robust to strs
    disk_required = sum((job.fileStore.getGlobalFileSize(gam) for gam in chunk_file_ids)) * 2 + (2 * 1024**3)
    requeued = ensure_disk_bytes(job, run_merge_chrom_gam, disk_required)
    if requeued is not None:
        # If not, requeue the job with more disk.
        return requeued
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    output_file = os.path.join(work_dir, '{}_{}.gam'.format(sample_name, chr_name))

    assert(len(chunk_file_ids) >= 1)

    if len(chunk_file_ids) > 1:
        # Would be nice to be able to do this merge with fewer copies.. 
        with open(output_file, 'ab') as merge_file:
            for chunk_gam_id in chunk_file_ids:
                tmp_gam_file = os.path.join(work_dir, 'tmp_{}.gam'.format(uuid4()))
                RealtimeLogger.info('Download file ID {} to {}'.format(chunk_gam_id, tmp_gam_file))
                job.fileStore.readGlobalFile(chunk_gam_id, tmp_gam_file)
                with open(tmp_gam_file, 'rb') as tmp_f:
                    shutil.copyfileobj(tmp_f, merge_file)
                                
    # checkpoint to out store
    if len(chunk_file_ids) == 1:
        job.fileStore.readGlobalFile(chunk_file_ids[0], output_file)
    return context.write_output_file(job, output_file)

def map_main(context, options):
    """
    Wrapper for vg map. 
    """

    validate_map_options(context, options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Make an index collection
            indexes = {}
           
            # Upload each index we have
            if options.xg_index is not None:
                indexes['xg'] = importer.load(options.xg_index)
            if options.gcsa_index is not None:
                indexes['gcsa'] = importer.load(options.gcsa_index)
                indexes['lcp'] = importer.load(options.gcsa_index + ".lcp")
            if options.gbwt_index is not None:
                indexes['gbwt'] = importer.load(options.gbwt_index)
            if options.distance_index is not None:
                indexes['distance'] = importer.load(options.distance_index)
            if options.minimizer_index is not None:
                indexes['minimizer'] = importer.load(options.minimizer_index)
            if options.snarls_index is not None:
                indexes['snarls'] = importer.load(options.snarls_index)
            if options.id_ranges is not None:
                indexes['id_ranges'] = importer.load(options.id_ranges)
            
            # Upload other local files to the remote IO Store
            inputReadsFileIDs = []
            if options.fastq:
                for sample_reads in options.fastq:
                    inputReadsFileIDs.append(importer.load(sample_reads))
            elif options.gam_input_reads:
                inputReadsFileIDs.append(importer.load(options.gam_input_reads))
            else:
                assert options.bam_input_reads
                inputReadsFileIDs.append(importer.load(options.bam_input_reads))

            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_mapping, context, options.fastq,
                                     options.gam_input_reads, options.bam_input_reads,
                                     options.sample_name,
                                     options.interleaved, options.mapper, importer.resolve(indexes),
                                     reads_file_ids=importer.resolve(inputReadsFileIDs),
                                     bam_output=options.bam_output, surject=options.surject,
                                     validate=options.validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
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
    
