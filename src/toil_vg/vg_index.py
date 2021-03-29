#!/usr/bin/env python
"""
vg_index.py: index a graph so it can be mapped to

"""

import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit, distutils.util
import logging, logging.handlers, struct, socket, threading
import string
import getpass
import pdb
import logging

from distutils import util

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def index_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # Options specific to the toil-vg index driver
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    parser.add_argument("--graphs", nargs='+', default=[], type=make_url,
                        help="input graph(s). one per chromosome (separated by space)")

    parser.add_argument("--chroms", nargs='+',
                        help="name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " has multiple elements, must be same length/order as --chroms (not needed for xg_index)")

    parser.add_argument("--node_mapping", type=make_url,
                        help="node mapping file required for gbwt pruning. Created by toil-vg construct"
                        " (or vg ids -j)")
                        
    parser.add_argument("--bwa_index_fasta", type=make_url,
                        help="index the given FASTA for BWA MEM alignment")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add indexing options
    index_toggle_parse_args(parser)
    index_parse_args(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)

def index_toggle_parse_args(parser):
    """
    Common args we do not want to share with toil-vg run, and which toggle
    different index types on and off.
    
    Safe to use in toil-vg construct without having to import any files.
    """
    parser.add_argument("--gcsa_index", dest="indexes",  default=[], action="append_const", const="gcsa",
                        help="Make a gcsa index for each output graph")
    parser.add_argument("--xg_index",  dest="indexes", action="append_const", const="xg",
                        help="Make an xg index for each output graph")
    parser.add_argument("--gbwt_index", dest="indexes", action="append_const", const="gbwt",
                        help="Make a GBWT index alongside the xg index for each output graph")
    parser.add_argument("--snarls_index",  dest="indexes", action="append_const", const="snarls",
                        help="Make an snarls file for each output graph")
    parser.add_argument("--trivial_snarls_index",  dest="indexes", action="append_const", const="trivial_snarls",
                        help="Make a trivial-inclusive snarls file for each output graph")
    parser.add_argument("--distance_index", dest="indexes", action="append_const", const="distance",
                        help="Make a (minimum) distance index for each output graph")
    parser.add_argument("--minimizer_index", dest="indexes", action="append_const", const="minimizer",
                        help="Make a minimizer index for each output graph")
    parser.add_argument("--id_ranges_index",  dest="indexes", action="append_const", const="id_ranges",
                        help="Make chromosome id ranges tables (so toil-vg map can optionally split output by chromosome)")
    parser.add_argument("--alt_path_gam_index", dest="indexes", action="append_const", const="alt-gam",
                        help="Save alt paths from vg into an indexed GAM")
    parser.add_argument("--xg_alts", dest="indexes", action="append_const", const="xg_alts",
                        help="Include alt paths in xg index")
    parser.add_argument("--all_index",  dest="indexes", action="store_const",
                        const=["gcsa", "xg", "gbwt", "snarls", "trivial_snarls", "distance", "minimizer", "id_ranges"],
                        help="Equivalent to --gcsa_index --xg_index --gbwt_index --snarls_index --trivial_snarls_index "
                        "--distance_index --minimizer_index --id_ranges_index")
    
def index_parse_args(parser):
    """
    Indexing parameters which can be part of the main toil-vg run pipeline.
    """
    
    parser.add_argument("--gcsa_index_cores", type=int,
        help="number of threads during the gcsa indexing step")
    parser.add_argument("--xg_index_cores", type=int,
        help="number of threads during the xg indexing step")
    parser.add_argument("--gbwt_index_cores", type=int,
        help="number of threads during the gbwt indexing step")    

    parser.add_argument("--index_name", type=str, default='index',
                        help="name of index files. <name>.xg, <name>.gcsa etc.")

    parser.add_argument("--gcsa_opts", type=str,
                        help="Options to pass to gcsa indexing.")
                        
    parser.add_argument("--minimizer_opts", type=str,
                        help="Options to pass to minimizer indexing.")

    parser.add_argument("--vcf_phasing", nargs='+', type=make_url, default=[],
                        help="Import phasing information from VCF(s) into xg (or GBWT with --gbwt_index)")
    parser.add_argument("--vcf_phasing_regions", nargs='+', default=[],
                        help="Hint the relevant chrom:start-end regions to the GBWT indexer, for subregion graphs")
    parser.add_argument("--gbwt_input", type=make_url,
                        help="Use given GBWT for GCSA2 pruning")
    parser.add_argument("--gbwt_prune", action='store_true',
                        help="Use gbwt for gcsa pruning")
    parser.add_argument("--force_phasing", type=lambda x:bool(util.strtobool(x)), default=None,
                        help="If 'True', randomly phase unphased variants and discard unresolveable overlaps for GBWT")
                        
def validate_index_options(options):
    """
    Validate the index options semantics enforced by index only.
    Throw an error if an invalid combination of options has been selected.
    """
    if len(options.indexes) > 0:
        require(len(options.graphs) == 0 or options.chroms, '--chroms must be specified for --graphs')
        require(len(options.graphs) == 1 or len(options.chroms) == len(options.graphs),
                '--chroms and --graphs must have'
                ' same number of arguments if more than one graph specified if doing anything but xg indexing')
    require(any([len(options.indexes) > 0, 
                 options.bwa_index_fasta]),
            'one of --xg_index, --gcsa_index, --snarls_index, --trivial_snarls_index, --id_ranges_index, '
            '--gbwt_index, --minimizer_index, --distance_index, --all_index, --alt_path_gam_index or '
            '--bwa_index_fasta is required')
    require(not options.gbwt_prune or options.node_mapping,
            '--node_mapping required with --gbwt_prune')
    require('gbwt' not in options.indexes or not options.gbwt_input,
            'only one of --gbwt_index and --gbwt_input can be used at a time')
    if options.gbwt_input:
        require(options.gbwt_prune == 'gbwt', '--gbwt_prune required with --gbwt_input')
    validate_shared_index_options(options)

def validate_shared_index_options(options):
    """
    Validate the index options semantics enforced by index and construct.
    Throw an error if an invalid combination of options has been selected.
    """
    
    if options.vcf_phasing:
        require(all([vcf.endswith('.vcf.gz') for vcf in options.vcf_phasing]),
                'input phasing files must end with .vcf.gz')
    if 'gbwt' in options.indexes:
        require(len(options.vcf_phasing) > 0, 'generating a GBWT requires a VCF with phasing information')
    if options.gbwt_prune:
        require(('gbwt' in options.indexes) or options.gbwt_input, '--gbwt_index or --gbwt_input required for --gbwt_prune')
    if options.vcf_phasing_regions:
        require('gbwt' in options.indexes, "cannot hint regions to GBWT indexer without building a GBWT index")
    
def run_gcsa_prune(job, context, graph_name, input_graph_id, gbwt_id, mapping_id, remove_paths = []):
    """
    Make a pruned graph using vg prune.  If unfold_mapping_id is provided, use -u, else -r
    """
    RealtimeLogger.info("Starting GCSA graph-pruning {}...".format("using GBWT" if gbwt_id else ""))
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Intermediate output
    restored_filename = os.path.join(work_dir, "restored_{}".format(graph_name))
    # Final output
    pruned_filename = os.path.join(work_dir, "unfolded_{}".format(graph_name))
    # Node Mapping output
    mapping_filename = os.path.join(work_dir, 'node_mapping')

    # Download input 
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(input_graph_id, graph_filename)
    gbwt_filename = graph_filename + '.gbwt'
    if gbwt_id:
        job.fileStore.readGlobalFile(gbwt_id, gbwt_filename)    
    if mapping_id:
        job.fileStore.readGlobalFile(mapping_id, mapping_filename, mutable=True)

    if remove_paths:
        # Remove alt paths so they don't make the graph too complex when getting restored
        remove_cmd = ['vg', 'mod', os.path.basename(graph_filename), '-I']
        for remove_path in remove_paths:
            remove_cmd += ['--retain-path', remove_path]
        cmd = [remove_cmd, ['vg', 'prune', '-']]
    else:
        cmd = [['vg', 'prune', os.path.basename(graph_filename)]]
    
    cmd[-1] += ['--threads', str(job.cores)]
    if context.config.prune_opts:
        cmd[-1] += context.config.prune_opts
    if gbwt_id:
        cmd[-1] += ['--append-mapping', '--mapping', os.path.basename(mapping_filename), '--unfold-paths']
        cmd[-1] += ['--gbwt-name', os.path.basename(gbwt_filename)]
    else:
        cmd[-1] += ['--restore-paths']
        
    with open(pruned_filename, 'wb') as pruned_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=pruned_file)
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA pruning. Process took {} seconds.".format(run_time))

    pruned_graph_id = context.write_intermediate_file(job, pruned_filename)
    if gbwt_id:
        mapping_id = context.write_intermediate_file(job, mapping_filename)
    else:
        mapping_id = None

    return pruned_graph_id, mapping_id

def run_gcsa_prep(job, context, input_graph_ids,
                  graph_names, index_name,
                  chrom_gbwt_ids, node_mapping_id, skip_pruning=False,
                  remove_paths=[]):
    """
    Do all the preprocessing for gcsa indexing (pruning)
    Then launch the indexing as follow-on
    """    
    RealtimeLogger.info("Starting gcsa preprocessing...")
    start_time = timeit.default_timer()
    if chrom_gbwt_ids:
        assert len(chrom_gbwt_ids) <= len(input_graph_ids)

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
    prune_root_job = Job()
    child_job.addChild(prune_root_job)

    # keep these in lists for now just in case
    prune_jobs = []
    # todo: figure out how best to update file with toil without making copies
    mapping_ids = [node_mapping_id] if node_mapping_id and chrom_gbwt_ids else []

    # prune each input graph.
    prune_ids = []
    for graph_i, input_graph_id in enumerate(input_graph_ids):
        gbwt_id = chrom_gbwt_ids[graph_i] if chrom_gbwt_ids else None
        mapping_id = mapping_ids[-1] if mapping_ids and gbwt_id else None
        # toggle between parallel/sequential based on if we're unfolding or not
        add_fn = prune_jobs[-1].addFollowOnJobFn if prune_jobs and gbwt_id else prune_root_job.addChildJobFn
        if not skip_pruning:
            prune_job = add_fn(run_gcsa_prune, context, graph_names[graph_i],
                               input_graph_id, gbwt_id, mapping_id,
                               remove_paths = remove_paths,
                               cores=context.config.prune_cores,
                               memory=context.config.prune_mem,
                               disk=context.config.prune_disk)
            prune_id = prune_job.rv(0)
            if gbwt_id:
                prune_jobs.append(prune_job)
                mapping_ids.append(prune_job.rv(1))
        else:
            prune_id = input_graph_id
        prune_ids.append(prune_id)

    return child_job.addFollowOnJobFn(run_gcsa_indexing, context, prune_ids,
                                      graph_names, index_name, mapping_ids[-1] if mapping_ids else None,
                                      cores=context.config.gcsa_index_cores,
                                      memory=context.config.gcsa_index_mem,
                                      disk=context.config.gcsa_index_disk,
                                      preemptable=context.config.gcsa_index_preemptable).rv()
    
def run_gcsa_indexing(job, context, prune_ids, graph_names, index_name, mapping_id):
    """
    Make the gcsa2 index. Return its store id
    """
    
    RealtimeLogger.info("Starting gcsa indexing...")
    start_time = timeit.default_timer()     

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Scratch directory for indexing
    index_temp_dir = os.path.join(work_dir, 'index-temp')
    os.makedirs(index_temp_dir)

    # Download all the pruned graphs.  
    prune_filenames = []
    
    for graph_i, prune_id in enumerate(prune_ids):
        prune_filename = os.path.join(work_dir, remove_ext(os.path.basename(graph_names[graph_i]), '.vg') + '.prune.vg')
        job.fileStore.readGlobalFile(prune_id, prune_filename)
        prune_filenames.append(prune_filename)

    # Download the mapping_id
    mapping_filename = None
    if mapping_id:
        mapping_filename = os.path.join(work_dir, 'node_mapping')
        job.fileStore.readGlobalFile(mapping_id, mapping_filename)

    # Where do we put the GCSA2 index?
    gcsa_filename = "{}.gcsa".format(index_name)

    command = ['vg', 'index', '-g', os.path.basename(gcsa_filename)] + context.config.gcsa_opts
    command += ['--threads', str(job.cores)]
    command += ['--temp-dir', os.path.join('.', os.path.basename(index_temp_dir))]
    
    if mapping_id:
        command += ['--mapping', os.path.basename(mapping_filename)]

    for prune_filename in prune_filenames:
        command += [os.path.basename(prune_filename)]

    try:
        context.runner.call(job, command, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the index run
        logging.error("GCSA indexing failed. Dumping files.")
        for prune_filename in prune_filenames:
            context.write_output_file(job, prune_filename)
        if mapping_id:
            context.write_output_file(job, mapping_filename)
        raise

    # Checkpoint index to output store
    gcsa_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename))
    lcp_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename) + ".lcp")

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA index. Process took {} seconds.".format(run_time))

    return gcsa_file_id, lcp_file_id

def run_concat_vcfs(job, context, vcf_ids, tbi_ids):
    """
    concatenate a list of vcfs.  we do this because vg index -v only takes one vcf, and
    we may be working with a set of chromosome vcfs. 
    """

    work_dir = job.fileStore.getLocalTempDir()

    vcf_names = ['chrom_{}.vcf.gz'.format(i) for i in range(len(vcf_ids))]
    out_name = 'genome.vcf.gz'

    for vcf_id, tbi_id, vcf_name in zip(vcf_ids, tbi_ids, vcf_names):
        job.fileStore.readGlobalFile(vcf_id, os.path.join(work_dir, vcf_name))
        job.fileStore.readGlobalFile(tbi_id, os.path.join(work_dir, vcf_name + '.tbi'))

    cmd = ['bcftools', 'concat'] + [vcf_name for vcf_name in vcf_names] + ['-O', 'z']
    
    with open(os.path.join(work_dir, out_name), 'wb') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)

    cmd = ['tabix', '-f', '-p', 'vcf', out_name]
    context.runner.call(job, cmd, work_dir=work_dir)

    out_vcf_id = context.write_intermediate_file(job, os.path.join(work_dir, out_name))
    out_tbi_id = context.write_intermediate_file(job, os.path.join(work_dir, out_name + '.tbi'))

    return out_vcf_id, out_tbi_id

 

def run_combine_graphs(job, context, inputGraphFileIDs, graph_names, index_name, intermediate=False):
    """
    Merge a list of graph files. We do this because the haplotype index vg index
    produces can only be built for a single graph.
    
    Takes the file IDs to concatenate, the human-readable names for those
    files, and the base name of the index we are working on (which is used to
    derive the graph output name).
    
    Graph files to concatenate must already be in the same ID space.
    
    If intermediate is set to true, the concatenated graph is not written to
    the output store.
    
    """
    
    # Define work directory for local files
    work_dir = job.fileStore.getLocalTempDir()
    
    RealtimeLogger.info("Starting VG graph merge...")
    start_time = timeit.default_timer()
    
    # The file names we are given can be very long, so if we download and cat
    # everything we can run into maximum command line length limits.
    
    # Unfortuantely, we need to use vg to do the graph combining because who
    # knows what HandleGraph format each file is in.
    
    # So download the files to short names.
    filenames = []
    for number, in_id in enumerate(inputGraphFileIDs):
        # Determine where to save the graph
        filename = '{}.vg'.format(number)
        
        # Put it in the workdir
        full_filename = os.path.join(work_dir, filename)
        
        # Save to the given file
        got_filename = job.fileStore.readGlobalFile(in_id, full_filename)
        
        RealtimeLogger.info('Downloaded graph ID {} to {} (which should be {}) for joining'.format(in_id, got_filename, full_filename))
        
        # Keep the filename
        filenames.append(filename)
        
    # Work out the file name we want to report
    concatenated_basename = "{}.cat.vg".format(index_name)
    
    # Run vg to combine into that file
    cmd = ['vg', 'combine'] + filenames
    
    try:
        with open(os.path.join(work_dir, concatenated_basename), 'wb') as out_file:
            context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)
    except:
        # Dump everything we need to replicate the index run
        logging.error("Graph merging failed. Dumping files.")

        for graph_filename in filenames:
            context.write_output_file(job, os.path.join(work_dir, graph_filename))
        
        raise
    
    # Now we generate the concatenated file ID
    concatenated_file_id = None
    if intermediate:
        # Save straight to the file store
        concatenated_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, concatenated_basename))
    else:
        # Checkpoint concatednated graph file to output store
        concatenated_file_id = context.write_output_file(job, os.path.join(work_dir, concatenated_basename))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished VG graph merge. Process took {} seconds.".format(run_time))

    return (concatenated_file_id, concatenated_basename)

def run_xg_indexing(job, context, inputGraphFileIDs, graph_names, index_name,
                    vcf_phasing_file_id = None, tbi_phasing_file_id = None,
                    make_gbwt=False, gbwt_regions=[],
                    intermediate=False, include_alt_paths=False):
    """

    Make the xg index and optional GBWT haplotype index.
    
    Saves the xg in the outstore as <index_name>.xg and the GBWT, if requested,
    as <index_name>.gbwt.
    
    If gbwt_regions is specified, it is a list of chrom:start-end region
    specifiers, restricting, on each specified chromosome, the region of the
    VCF that GBWT indexing will examine.
            
    if make_gbwt is specified *and* a phasing VCF is specified, the GBWT will
    be generated. Otherwise it won't be (for example, for single-contig graphs
    where no VCF is available).
    
    Return a tuple of file IDs, (xg_id, gbwt_id, thread_db_id). The GBWT ID
    will be None if no GBWT is generated. The thread DB ID will be None if no
    thread DB is generated.
    
    If intermediate is set to true, do not save the produced files to the
    output store.
    """
    
    RealtimeLogger.info("Starting xg indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Scratch directory for indexing
    index_temp_dir = os.path.join(work_dir, 'index-temp')
    os.makedirs(index_temp_dir)
    
    RealtimeLogger.info("inputGraphFileIDs: {}".format(str(inputGraphFileIDs)))
    RealtimeLogger.info("graph_names: {}".format(str(graph_names)))
    # Our local copy of the graphs
    graph_filenames = []
    for i, graph_id in enumerate(inputGraphFileIDs):
        graph_filename = os.path.join(work_dir, graph_names[i])
        job.fileStore.readGlobalFile(graph_id, graph_filename)
        graph_filenames.append(os.path.basename(graph_filename))

    # If we have a separate GBWT it will go here
    gbwt_filename = os.path.join(work_dir, "{}.gbwt".format(index_name))
    # And if we ahve a separate thread db it will go here
    thread_db_filename = os.path.join(work_dir, "{}.threads".format(index_name))
    
    # Get the vcf file for loading phasing info
    if vcf_phasing_file_id:
        phasing_file = os.path.join(work_dir, 'phasing.{}.vcf.gz'.format(index_name))
        job.fileStore.readGlobalFile(vcf_phasing_file_id, phasing_file)
        job.fileStore.readGlobalFile(tbi_phasing_file_id, phasing_file + '.tbi')
        phasing_opts = ['-v', os.path.basename(phasing_file)]
        
        if make_gbwt:
            # Write the haplotype index to its own file
            phasing_opts += ['--gbwt-name', os.path.basename(gbwt_filename)]
                       
            for region in gbwt_regions:
                phasing_opts += ['--region', region]

            if context.config.force_phasing:
                # We need to discard overlaps also to really get rid of haplotype breaks.
                phasing_opts += ['--force-phasing', '--discard-overlaps']
    else:
        phasing_opts = []
        
    # Where do we put the XG index?
    xg_filename = "{}.xg".format(index_name)

    # Now run the indexer.
    RealtimeLogger.info("XG Indexing {}".format(str(graph_filenames)))

    command = ['vg', 'index', '--threads', str(job.cores), '--xg-name', os.path.basename(xg_filename)]
    command += phasing_opts + graph_filenames
    command += ['--temp-dir', os.path.join('.', os.path.basename(index_temp_dir))]

    if include_alt_paths:
        command += ['--xg-alts']
    
    try:
        context.runner.call(job, command, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the index run
        logging.error("XG indexing failed. Dumping files.")

        for graph_filename in graph_filenames:
            context.write_output_file(job, os.path.join(work_dir, graph_filename))
        if vcf_phasing_file_id:
            context.write_output_file(job, phasing_file)
            context.write_output_file(job, phasing_file + '.tbi')

        raise

    # Determine if we want to checkpoint index to output store
    write_function = context.write_intermediate_file if intermediate else context.write_output_file
    xg_file_id = write_function(job, os.path.join(work_dir, xg_filename))
    
    gbwt_file_id = None
    thread_db_file_id = None
    if make_gbwt and vcf_phasing_file_id:
        # Also save the GBWT if it was generated
        gbwt_file_id = write_function(job, gbwt_filename)
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished XG index. Process took {} seconds.".format(run_time))

    # TODO: convert to a dict
    return (xg_file_id, gbwt_file_id)

def run_cat_xg_indexing(job, context, inputGraphFileIDs, graph_names, index_name,
                        vcf_phasing_file_id = None, tbi_phasing_file_id = None,
                        make_gbwt=False, gbwt_regions=[], 
                        intermediate=False, intermediate_cat=True, include_alt_paths=False):
    """
    Encapsulates run_combine_graphs and run_xg_indexing job functions.
    Can be used for ease of programming in job functions that require running only
    during runs of the run_xg_indexing job function.

    Note: the resources assigned to indexing come from those assigned to this parent job
    (as they can get toggled between xg and gbwt modes in the caller)
    
    If intermediate is set to True, do not save the final .xg to the output store.
    
    If intermediate_cat is False and intermediate is also False, save the .cat.vg to the output store.
    """
    
    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)    
    
    # Concatenate the graph files.
    vg_concat_job = child_job.addChildJobFn(run_combine_graphs, context, inputGraphFileIDs,
                                            graph_names, index_name, intermediate=(intermediate or intermediate_cat),
                                            cores=job.cores,
                                            memory=job.memory,
                                            disk=job.disk)
    
    return child_job.addFollowOnJobFn(run_xg_indexing,
                                      context, [vg_concat_job.rv(0)],
                                      [vg_concat_job.rv(1)], index_name,
                                      vcf_phasing_file_id, tbi_phasing_file_id,
                                      make_gbwt=make_gbwt, gbwt_regions=gbwt_regions,
                                      intermediate=intermediate,
                                      include_alt_paths=include_alt_paths,
                                      cores=job.cores,
                                      memory=job.memory,
                                      disk=job.disk,
                                      preemptable=job.preemptable).rv()
                                      
def run_snarl_indexing(job, context, inputGraphFileIDs, graph_names, index_name=None, include_trivial=False):
    """
    Compute the snarls of the graph.
    
    Saves the snarls file in the outstore as <index_name>.snarls or
    <index_name>.trivial.snarls, as appropriate, unless index_name is None.
    
    If incluse_trivial is set to true, include trivial snarls, which mpmap
    cannot yet filter out itself for snarl cutting, but which are needed for
    distance indexing.
    
    Return the file ID of the snarls file.
    """
    
    assert(len(inputGraphFileIDs) == len(graph_names))
    
    # Decide on an index output extension.
    extension = '.trivial.snarls' if include_trivial else '.snarls'
    
    if len(inputGraphFileIDs) > 1:
        # We have been given multiple chromosome graphs. Since snarl indexing
        # can take a lot of memory, we are going to process each one separately
        # and then concatenate the results.
        
        RealtimeLogger.info("Breaking up snarl computation for {}".format(str(graph_names)))
        
        snarl_jobs = []
        for file_id, file_name in zip(inputGraphFileIDs, graph_names):
            # For each input graph, make a child job to index it.
            snarl_jobs.append(job.addChildJobFn(run_snarl_indexing, context, [file_id], [file_name],
                                                include_trivial=include_trivial,
                                                cores=context.config.snarl_index_cores,
                                                memory=context.config.snarl_index_mem,
                                                disk=context.config.snarl_index_disk))
                                                
        # Make a job to concatenate the indexes all together                                        
        concat_job = snarl_jobs[0].addFollowOnJobFn(run_concat_files, context, [job.rv() for job in snarl_jobs],
                                                    index_name + extension if index_name is not None else None,
                                                    cores=context.config.snarl_index_cores,
                                                    memory=context.config.snarl_index_mem,
                                                    disk=context.config.snarl_index_disk)
        
        for i in range(1, len(snarl_jobs)):
            # And make it wait for all of them
            snarl_jobs[i].addFollowOn(concat_job)
            
        return concat_job.rv()
        
    else:
        # Base case: single graph
   
        RealtimeLogger.info("Starting snarl computation {} trivial snarls...".format('with' if include_trivial else 'without'))
        start_time = timeit.default_timer()
        
        # Define work directory for docker calls
        work_dir = job.fileStore.getLocalTempDir()

        # Download the one graph
        graph_id = inputGraphFileIDs[0]
        graph_filename = graph_names[0]
        job.fileStore.readGlobalFile(graph_id, os.path.join(work_dir, graph_filename))

        # Where do we put the snarls?
        snarl_filename = os.path.join(work_dir, (index_name if index_name is not None else "part") + extension)

        # Now run the indexer.
        RealtimeLogger.info("Computing snarls for {}".format(graph_filename))

        cmd = ['vg', 'snarls', graph_filename]
        if include_trivial:
            cmd += ['--include-trivial']
        with open(snarl_filename, "wb") as snarl_file:
            try:
                # Compute snarls to the correct file
                context.runner.call(job, cmd, work_dir=work_dir, outfile=snarl_file)
            except:
                # Dump everything we need to replicate the indexing
                logging.error("Snarl indexing failed. Dumping files.")
                context.write_output_file(job, os.path.join(work_dir, graph_filename))
                raise
        
        if index_name is not None:
            # Checkpoint index to output store
            snarl_file_id = context.write_output_file(job, snarl_filename)
        else:
            # Just save the index as an intermediate
            snarl_file_id = context.write_intermediate_file(job, snarl_filename)
            
        
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Finished computing snarls. Process took {} seconds.".format(run_time))

        return snarl_file_id
        
def run_distance_indexing(job, context, input_xg_id, input_trivial_snarls_id, index_name=None, max_distance_threshold=0):
    """
    Make a distance index from the given XG index and the given snarls file,
    including the trivial snarls.
    
    TODO: also support a single VG file and snarls.
    
    If index_name is not None, saves it as <index_name>.dist to the output
    store.
    
    If max_distance_threshold is strictly positive, also includes a max
    distance index with the given limit. By default includes only a min
    distance index.
    
    Returns the file ID of the resulting distance index.
    """
    
    RealtimeLogger.info("Starting distance indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download the input files.
    xg_filename = os.path.join(work_dir, 'graph.xg')
    job.fileStore.readGlobalFile(input_xg_id, xg_filename)
    trivial_snarls_filename = os.path.join(work_dir, 'graph.trivial.snarls')
    job.fileStore.readGlobalFile(input_trivial_snarls_id, trivial_snarls_filename)

    # Where do we put the distance index?
    distance_filename = os.path.join(work_dir, (index_name if index_name is not None else 'graph') + '.dist')

    cmd = ['vg', 'index', '-t', max(1, int(job.cores)), '-j', os.path.basename(distance_filename),
        '-x', os.path.basename(xg_filename), '-s', os.path.basename(trivial_snarls_filename)]
    
    if max_distance_threshold > 0:
        # Add a max distance index with this limit.
        cmd.append('-w')
        cmd.append(str(max_distance_threshold))
    
    try:
        # Compute the index to the correct file
        context.runner.call(job, cmd, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the indexing
        logging.error("Distance indexing failed. Dumping files.")
        context.write_output_file(job, xg_filename)
        context.write_output_file(job, trivial_snarls_filename)
        if os.path.exists(distance_filename):
            context.write_output_file(job, distance_filename)
        raise
    
    if index_name is not None:
        # Checkpoint index to output store
        distance_file_id = context.write_output_file(job, distance_filename)
    else:
        # Just save the index as an intermediate
        distance_file_id = context.write_intermediate_file(job, distance_filename)
        
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished computing distance index. Process took {} seconds.".format(run_time))

    return distance_file_id
        
def run_minimizer_indexing(job, context, input_xg_id, input_gbwt_id, index_name=None):
    """
    Make a minimizer index file for the graph and haplotypes described by the
    given input XG and GBWT indexes.
    
    If index_name is not None, saves it as <index_name>.min to the output
    store.
    
    Returns the file ID of the resulting minimizer index.
    """
    
    RealtimeLogger.info("Starting minimizer indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Check our setup.
    # TODO: Be able to build a minimizer index for all paths somehow if the
    # GBWT isn't available, for e.g. a linear graph.
    assert input_xg_id is not None
    assert input_gbwt_id is not None

    # Download the input files.
    xg_filename = os.path.join(work_dir, 'graph.xg')
    job.fileStore.readGlobalFile(input_xg_id, xg_filename)
    gbwt_filename = os.path.join(work_dir, 'graph.gbwt')
    job.fileStore.readGlobalFile(input_gbwt_id, gbwt_filename)

    # Where do we put the minimizer index?
    minimizer_filename = os.path.join(work_dir, (index_name if index_name is not None else 'graph') + '.min')

    cmd = ['vg', 'minimizer', '-t', max(1, int(job.cores)), '-i', os.path.basename(minimizer_filename),
        '-g', os.path.basename(gbwt_filename)] + context.config.minimizer_opts + [os.path.basename(xg_filename)]
    try:
        # Compute the index to the correct file
        context.runner.call(job, cmd, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the indexing
        logging.error("Minimizer indexing failed. Dumping files.")
        context.write_output_file(job, xg_filename)
        context.write_output_file(job, gbwt_filename)
        context.write_output_file(job, minimizer_filename)
        raise
    
    if index_name is not None:
        # Checkpoint index to output store
        minimizer_file_id = context.write_output_file(job, minimizer_filename)
    else:
        # Just save the index as an intermediate
        minimizer_file_id = context.write_intermediate_file(job, minimizer_filename)
        
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished computing minimizer index. Process took {} seconds.".format(run_time))

    return minimizer_file_id


def run_id_ranges(job, context, inputGraphFileIDs, graph_names, index_name, chroms):
    """ Make a file of chrom_name <tab> first_id <tab> last_id covering the 
    id ranges of all chromosomes.  This is to speed up gam splitting down the road.
    
    Chroms is a list of contig names or sets of contig names that correspond to each graph.
    """
    
    RealtimeLogger.info("Starting id ranges...")
    start_time = timeit.default_timer()
    
    # Our id ranges (list of triples)
    id_ranges = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)    

    # Get the range for one graph per job. 
    for graph_id, graph_name, chrom in zip(inputGraphFileIDs, graph_names, chroms):
        id_range = child_job.addChildJobFn(run_id_range, context, graph_id, graph_name, chrom,
                                           cores=context.config.prune_cores,
                                           memory=context.config.prune_mem, disk=context.config.prune_disk).rv()
        
        id_ranges.append(id_range)

    # Merge them into a file and return its id
    return child_job.addFollowOnJobFn(run_merge_id_ranges, context, id_ranges, index_name,
                                      cores=context.config.misc_cores, memory=context.config.misc_mem,
                                      disk=context.config.misc_disk).rv()

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished id ranges. Process took {} seconds.".format(run_time))
    
def run_id_range(job, context, graph_id, graph_name, chrom):
    """
    Compute a node id range for a graph (which should be an entire contig/chromosome with
    contiguous id space -- see vg ids) using vg stats
    
    Chrom is a contig name or set of contig names we expect to be in the graph.
    
    If multiple contigs are in the graph, we comma-separate them.
    
    Returns a tuple of 3 strings: chromosome name, first base number, last base number.
    """
    work_dir = job.fileStore.getLocalTempDir()

    # download graph
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(graph_id, graph_filename)

    #run vg stats
    #expect result of form node-id-range <tab> first:last
    command = ['vg', 'stats', '--node-id-range', os.path.basename(graph_filename)]
    stats_out = context.runner.call(job, command, work_dir=work_dir, check_output = True).decode('utf-8').strip().split()
    assert stats_out[0] == 'node-id-range'
    first, last = stats_out[1].split(b':')
    
    if isinstance(chrom, set):
        chrom = ','.join(sorted(chrom))
    
    return chrom, first, last
    
def run_merge_id_ranges(job, context, id_ranges, index_name):
    """ create a BED-style file of id ranges
    """
    work_dir = job.fileStore.getLocalTempDir()

    # Where do we put the id ranges tsv?
    id_range_filename = os.path.join(work_dir, '{}_id_ranges.tsv'.format(index_name))

    with open(id_range_filename, 'w') as f:
        for id_range in id_ranges:
            f.write('{}\t{}\t{}\n'.format(*id_range))

    # Checkpoint index to output store
    return context.write_output_file(job, id_range_filename)

def run_merge_gbwts(job, context, chrom_gbwt_ids, index_name):
    """ merge up some gbwts
    """
    work_dir = job.fileStore.getLocalTempDir()

    gbwt_chrom_filenames = []

    for i, gbwt_id in enumerate(chrom_gbwt_ids):
        if gbwt_id:
            gbwt_filename = os.path.join(work_dir, '{}.gbwt'.format(i))
            job.fileStore.readGlobalFile(gbwt_id, gbwt_filename)
            gbwt_chrom_filenames.append(gbwt_filename)

    if len(gbwt_chrom_filenames) == 0:
        return None
    elif len(gbwt_chrom_filenames) == 1:
        return context.write_output_file(job, gbwt_chrom_filenames[0],
                                         out_store_path = index_name + '.gbwt')
    else:
        # Merge the GBWT files together
        cmd = ['vg', 'gbwt', '--merge', '--fast', '--output', index_name + '.gbwt']
        cmd += [os.path.basename(f) for f in gbwt_chrom_filenames]
        
        try:
            context.runner.call(job, cmd, work_dir=work_dir)
        except:
            # Dump everything we need to replicate the merge
            logging.error("GBWT merge failed. Dumping files.")
            for f in gbwt_chrom_filenames:
                context.write_output_file(job, f)
            
            raise

        return context.write_output_file(job, os.path.join(work_dir, index_name + '.gbwt'))
        
def run_bwa_index(job, context, fasta_file_id, bwa_index_ids=None, intermediate=False, copy_fasta=False):
    """
    Make a bwa index for a fast sequence if not given in input.
    
    If intermediate is set to True, do not output them. Otherwise, output them
    as bwa.fa.<index type>.
    
    Returns a dict from index extension to index file ID.
    
    Note that BWA produces 'amb', 'ann',  'bwt', 'pac', and 'sa' index files.
    
    If such a nonempty dict is passed in already, return that instead (and
    don't output any files).
    
    If copy_fasta is True (and intermediate is False), also output the input FASTA to the out store.
    
    """
    if not bwa_index_ids:
        bwa_index_ids = dict()
        work_dir = job.fileStore.getLocalTempDir()
        # Download the FASTA file to be indexed
        # It would be nice to name it the same as the actual input FASTA but we'd have to peek at the options
        fasta_file = os.path.join(work_dir, 'bwa.fa')
        job.fileStore.readGlobalFile(fasta_file_id, fasta_file)
        cmd = ['bwa', 'index', os.path.basename(fasta_file)]
        context.runner.call(job, cmd, work_dir = work_dir)
        
        # Work out how to output the files
        write_file = context.write_intermediate_file if intermediate else context.write_output_file
        
        for idx_file in glob.glob('{}.*'.format(fasta_file)):
            # Upload all the index files created, and store their IDs under their extensions
            bwa_index_ids[idx_file[len(fasta_file):]] = write_file(job, idx_file)
            
        if copy_fasta and not intermediate:
            # We ought to upload the FASTA also.
            context.write_output_file(job, fasta_file)

    return bwa_index_ids
    
def run_minimap2_index(job, context, fasta_file_id, minimap2_index_id=None, intermediate=False, copy_fasta=False):
    """
    Make a minimap2 index for a fasta sequence if not given in input.
    
    If intermediate is set to True, do not output it. Otherwise, output it
    as minimap2.fa.mmi.
    
    Returns the index file ID.
    
    If copy_fasta is True (and intermediate is False), also output the input FASTA to the out store.
    
    """
    if not minimap2_index_id:
        work_dir = job.fileStore.getLocalTempDir()
        # Download the FASTA file to be indexed
        # It would be nice to name it the same as the actual input FASTA but we'd have to peek at the options
        fasta_file = os.path.join(work_dir, 'minimap2.fa')
        job.fileStore.readGlobalFile(fasta_file_id, fasta_file)
        
        # Say where the index should go
        index_file = os.path.join(work_dir, 'minimap2.fa.mmi')
        
        # Make the index
        cmd = ['minimap2', '-d', os.path.basename(index_file), os.path.basename(fasta_file)]
        context.runner.call(job, cmd, work_dir = work_dir)
        
        # Work out how to output the files
        write_file = context.write_intermediate_file if intermediate else context.write_output_file
        
        minimap2_index_id = write_file(job, index_file)
        
        if copy_fasta and not intermediate:
            # We ought to upload the FASTA also.
            context.write_output_file(job, fasta_file)

    return minimap2_index_id
        
def run_alt_path_extraction(job, context, inputGraphFileIDs, graph_names, index_name):
    """
    Pull the alt paths out of the graph and into a GAM index.

    The application is to be able to get them into vg call by way of vg chunk and vg augment.

    Hopefully, this is a stopgap and we can eventually fix up xg to handle them efficiently.
    
    Return the file ID of the GAM
    """
    
    assert(len(inputGraphFileIDs) == len(graph_names))
    
    if len(inputGraphFileIDs) > 1:
        # We have been given multiple chromosome graphs. 
        
        RealtimeLogger.info("Breaking up alt path GAM computation for {}".format(str(graph_names)))
        
        sub_jobs = []
        for i, (file_id, file_name) in enumerate(zip(inputGraphFileIDs, graph_names)):
            # For each input graph, make a child job to index it.
            sub_jobs.append(job.addChildJobFn(run_alt_path_extraction, context, [file_id], [file_name],
                                              index_name + '.{}'.format(i) if index_name else None,
                                              cores=context.config.chunk_cores,
                                              memory=context.config.chunk_mem,
                                              disk=context.config.chunk_disk))
        
        # Make a job to concatenate the indexes all together                                        
        concat_job = sub_jobs[0].addFollowOnJobFn(run_concat_files, context, [job.rv() for job in sub_jobs],
                                                  index_name + '_alts.gam' if index_name is not None else None,
                                                  memory=context.config.chunk_mem,
                                                  disk=context.config.chunk_disk)
        
        for i in range(1, len(sub_jobs)):
            # And make it wait for all of them
            sub_jobs[i].addFollowOn(concat_job)
            
        return concat_job.rv()
        
    else:
        # Base case: single graph
   
        start_time = timeit.default_timer()
        
        # Define work directory for docker calls
        work_dir = job.fileStore.getLocalTempDir()

        # Download the one graph
        graph_id = inputGraphFileIDs[0]
        graph_filename = graph_names[0]
        job.fileStore.readGlobalFile(graph_id, os.path.join(work_dir, graph_filename))

        # Where do we put the gam?
        gam_filename = os.path.join(work_dir, "{}_alts.gam".format(index_name if index_name is not None else "part"))

        cmd = ['vg', 'paths', '-v', graph_filename, '-Q', '_alt_', '-X']
        with open(gam_filename, 'wb') as gam_file:
            try:
                # Compute snarls to the correct file
                context.runner.call(job, cmd, work_dir=work_dir, outfile=gam_file)
            except:
                # Dump everything we need to replicate the indexing
                logging.error("Alt path gam extraction failed. Dumping files.")
                context.write_output_file(job, os.path.join(work_dir, graph_filename))
                raise
        
        if index_name is not None:
            # Checkpoint index to output store
            gam_file_id = context.write_output_file(job, gam_filename)
        else:
            # Just save the index as an intermediate
            gam_file_id = context.write_intermediate_file(job, gam_filename)
            
        
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Finished GAM extraction. Process took {} seconds.".format(run_time))

        return gam_file_id

def run_gam_indexing(job, context, gam_id, index_name):
    """ Index a gam.  Return the sorted gam and its .gai index.
    """
    work_dir = job.fileStore.getLocalTempDir()
    gam_filename = os.path.join(work_dir, "{}_alts.gam".format(index_name if index_name is not None else "index"))
    job.fileStore.readGlobalFile(gam_id, gam_filename + '.unsorted')

    cmd = ['vg', 'gamsort', os.path.basename(gam_filename) + '.unsorted',
           '-i', os.path.basename(gam_filename) + '.gai', '-t', str(job.cores)]
    with open(gam_filename, 'wb') as gam_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=gam_file)

    if index_name is not None:
        # Checkpoint index to output store
        sorted_gam_id = context.write_output_file(job, gam_filename)
        gai_id = context.write_output_file(job, gam_filename + '.gai')
    else:
        # Just save the index as an intermediate
        sorted_gam_id = context.write_intermediate_file(job, gam_filename)
        gai_id = context.write_output_file(job, gam_filename + '.gai')
    
    return sorted_gam_id, gai_id
    
    
def run_indexing(job, context, inputGraphFileIDs,
                 graph_names, index_name, chroms,
                 vcf_phasing_file_ids = [], tbi_phasing_file_ids = [],
                 bwa_fasta_id=None,
                 gbwt_id = None, node_mapping_id = None,
                 wanted = set(),
                 gbwt_prune=False, gbwt_regions=[],
                 dont_restore_paths=[],
                 coalesce_regions=[]):
    """
    
    Run indexing logic by itself.
    
    vcf_phasing_file_ids and tbi_phasing_file_ids are phasing data VCFs. There
    can be 0 of them, 1 for all chromosomes, or one for each chromosome in
    chroms order.
    
    gbwt_regions is a list of chrom:start-end regions pecifiers to restrict, on
    those chromosomes, the regions examined in the VCF by the GBWT indexing.
    
    wanted is a set of the index type strings ('xg', 'gcsa', 'gbwt',
    'id_ranges', 'snarls', 'trivial_snarls', 'minimizer', 'distance',
    'alt-gam') that should be created. Each of them becomes a key in the output
    dict, except that:
    
    * The 'bwa' index is produced if bwa_fasta_id is set instead of if 'bwa' is
    in wanted.
    
    * The 'lcp' index is produced whenever the 'gcsa' index is produced.
    
    * The 'id_ranges' index is only produced if multiple chromosomes are used.
    
    * The 'chrom_...' versions of indexes are created when the overall index is
    created.
    
    Return a dict from index type ('xg','chrom_xg', 'gcsa', 'lcp', 'gbwt',
    'chrom_gbwt', 'chrom_thread', 'id_ranges', 'snarls', 'trivial_snarls',
    'alt-gam', 'bwa') to index file ID(s) if created.
    
    For 'chrom_xg' and 'chrom_gbwt' the value is a list of one XG or GBWT or
    thread DB per chromosome in chroms, to support `vg prune`. For
    'chrom_thread', we have a value per chromosome that actually has any
    threads (instead of padding out with Nones). For 'bwa', the result is
    itself a dict from BWA index extension to file ID. The others are all
    single file IDs.
    
    If gbwt_id is specified, and the gbwt index is not built, the passed ID is
    re-used.
    
    If the 'gbwt' index is requested and gbwt_id is not specified, the 'xg'
    index will also be computed. If no phasing VCFs are provided, computing
    this index will be skipped.
    
    If the 'minimizer' index is requested, the 'xg' index will also be
    computed, and the 'gbwt' index will either be computed or sourced from
    gbwt_id. If the 'gbwt' index is not available, computing this index will be
    skipped.
    
    If the 'distance' index is requested, the 'trivial_snarls' and 'xg' indexes
    will also be computed.
    
    If coalesce_regions is set, it must be a list of sets of 'chroms' region
    names. Each set of region names will be expected to be together in a graph
    file, instead of in separate graph files.
    
    """
    
    # Coalesce the chroms, so we have some sets of chroms that live in the same
    # graph file.
    chroms, chrom_names = apply_coalesce(chroms, coalesce_regions=coalesce_regions)
    
    # Make a master child job
    child_job = Job()
    job.addChild(child_job)
    
    # And one job for all the per-chromosome xg jobs
    chrom_xg_root_job = Job()
    child_job.addChild(chrom_xg_root_job)
    
    # And inside it make one job for the main whole-graph xg construction that has to come after it
    xg_root_job = Job()
    chrom_xg_root_job.addFollowOn(xg_root_job)
    
    RealtimeLogger.info("Running indexing: {}.".format({
         'graph_names': graph_names,
         'index_name': index_name,
         'chroms': chroms if len(chroms) < 100 else f'{len(chroms)} items',
         'vcf_phasing_file_ids': vcf_phasing_file_ids,
         'tbi_phasing_file_ids': tbi_phasing_file_ids,
         'gbwt_id': gbwt_id,
         'node_mapping_id': node_mapping_id,
         'wanted': wanted,
         'gbwt_prune': gbwt_prune,
         'bwa_fasta_id': bwa_fasta_id,
         'coalesce_regions': coalesce_regions if max([len(x) for x in coalesce_regions] + [0]) < 100 else '(many)'
    }))

    # This will hold the index to return
    indexes = {}
    if gbwt_id:
        indexes['gbwt'] = gbwt_id
    elif 'gbwt' in wanted:
        # We need to do the xg so we can make the GBWT.
        # TODO: write a codepath that makes the GBWT without making the XG
        wanted.add('xg')

    # We shouldn't accept any phasing files when not making a GBWT index with them.
    assert(len(vcf_phasing_file_ids) == 0 or ('gbwt' in wanted))
    
    if 'minimizer' in wanted:
        # The minimizer index has some dependencies
        wanted.add('xg')
        if not gbwt_id:
            wanted.add('gbwt')
            
    if 'distance' in wanted:
        # The distance index also has some dependencies
        wanted.add('xg')
        wanted.add('trivial_snarls')
        
    # We guarantee that if 'gbwt' is in indexes, then there is (a promise for)
    # an actual GBWT.

    if 'xg' in wanted or 'gcsa' in wanted:
        indexes['chrom_xg'] = []
        indexes['chrom_gbwt'] = []
        
        if 'gbwt' in wanted and len(vcf_phasing_file_ids) > 0:
            # We want to make a GBWT, and we can in fact make a GBWT.
            # That's the only case we want per-chromosome XGs for anymore.
        
            # In its current state, vg prune requires chromosomal xgs, so we must make
            # these xgs if we're doing any kind of gcsa indexing.  Also, if we're making
            # a gbwt, we do that at the same time (merging later if more than one graph).
            #
            # TODO: This is the *only* way to make the GBWT, actually. Write
            # code to make it separately.
            if not chroms or len(chroms) == 1:
                chroms = [index_name]
            indexes['chrom_xg'] = []
            indexes['chrom_gbwt'] = []
            
            # Check our input phasing VCF set for plausibility
            if len(vcf_phasing_file_ids) != len(tbi_phasing_file_ids):
                # Each VCF needs an index
                raise RuntimeError("Found {} phasing VCFs and {} indexes; counts must match!".format(
                    len(vcf_phasing_file_ids), len(tbi_phasing_file_ids)))
                    
            if len(vcf_phasing_file_ids) > len(chroms):
                # We can only handle no VCFs, one VCF, or one VCF per chromosome until we run out of VCFs.
                # So what we can't handle is more VCFs than chromosomes
                RealtimeLogger.error("Chromosomes: {}".format(chroms))
                RealtimeLogger.error("VCFs: {}".format(vcf_phasing_file_ids))
                raise RuntimeError("Found too many ({}) phasing VCFs for {} chromosomes".format(
                    len(vcf_phasing_file_ids), len(chroms)))
            
            
            for i, chrom in enumerate(chroms):
                # For each chromosome
                
                # Find the phasing VCF
                if len(vcf_phasing_file_ids) == 0:
                    # There may be 0
                    vcf_id = None
                    tbi_id = None
                elif len(vcf_phasing_file_ids) == 1:
                    # There may be one for all chromosomes
                    vcf_id = vcf_phasing_file_ids[0]
                    tbi_id = tbi_phasing_file_ids[0]
                elif i < len(vcf_phasing_file_ids):
                    # Otherwise the VCFs and chromosomes correspond in order, until the VCFs are depleted.
                    # There is one for this chromosome
                    vcf_id = vcf_phasing_file_ids[i]
                    tbi_id = tbi_phasing_file_ids[i]
                else:
                    # We have run out of VCFs for chromosomes to be in
                    vcf_id = None
                    tbi_id = None
                
                # Make a job to index just this chromosome and produce a
                # per-chromosome xg, gbwt, and threads file. Since there may be
                # thousands of chromosomes (including e.g. decoys) in a
                # whole-genome reference, keep these files as intermediates and
                # don't put them in the outstore, unless we're only doing one contig.
                xg_chrom_index_job = chrom_xg_root_job.addChildJobFn(run_cat_xg_indexing,
                                                                     context, [inputGraphFileIDs[i]],
                                                                     [graph_names[i]], chrom_names[i],
                                                                     vcf_id, tbi_id,
                                                                     make_gbwt=('gbwt' in wanted),
                                                                     gbwt_regions=gbwt_regions, intermediate=(len(chroms) > 1),
                                                                     include_alt_paths=('xg_alts' in wanted),
                                                                     cores=context.config.gbwt_index_cores,
                                                                     memory=context.config.gbwt_index_mem,
                                                                     disk=context.config.gbwt_index_disk,
                                                                     preemptable='gbwt' not in wanted or context.config.gbwt_index_preemptable)
                indexes['chrom_xg'].append(xg_chrom_index_job.rv(0))
                indexes['chrom_gbwt'].append(xg_chrom_index_job.rv(1))

            if len(chroms) > 1:
                # Once all the per-chromosome GBWTs are done and we are ready to make the whole-graph GBWT, merge them up
                indexes['gbwt'] = xg_root_job.addChildJobFn(run_merge_gbwts, context, indexes['chrom_gbwt'],
                                                            index_name,
                                                            cores=context.config.xg_index_cores,
                                                            memory=context.config.xg_index_mem,
                                                            disk=context.config.xg_index_disk).rv()
            else:
                # There's only one chromosome, so the one per-chromosome GBWT becomes the only GBWT
                indexes['gbwt'] = indexes['chrom_gbwt'][0]
                
        # now do the whole genome xg (without any gbwt)
        if 'chrom_xg' in indexes and len(indexes['chrom_xg']) == 1:
            # We made per-chromosome XGs and we have exactly one.
            # our first chromosome is effectively the whole genome (note that above we
            # detected this and put in index_name so it's saved right (don't care about chrom names))
            indexes['xg'] = indexes['chrom_xg'][0]
        elif 'xg' in wanted:
            # Build an xg index for the whole genome.
            
            xg_index_job = xg_root_job.addChildJobFn(run_cat_xg_indexing,
                                                     context, inputGraphFileIDs,
                                                     graph_names, index_name,
                                                     None, None,
                                                     make_gbwt=False,
                                                     include_alt_paths=('xg_alts' in wanted),
                                                     cores=context.config.xg_index_cores,
                                                     memory=context.config.xg_index_mem,
                                                     disk=context.config.xg_index_disk)
            
            indexes['xg'] = xg_index_job.rv(0)


    gcsa_root_job = Job()
    # gcsa follows from chrom_xg jobs only if per-chromosome gbwts are needed for per-chromosome pruning
    if gbwt_prune:
        chrom_xg_root_job.addFollowOn(gcsa_root_job)
    else:
        child_job.addChild(gcsa_root_job)
    
    if 'gcsa' in wanted:
        # We know we made the per-chromosome indexes already, so we can use them here to make the GCSA                                               
        # todo: we're only taking in a genome gbwt as input, because that's all we write
        if ('chrom_gbwt' not in indexes or indexes['chrom_gbwt'] == []) and 'gbwt' in indexes:
            # We lack per-chromosome GBWTs but we have a whole genome one we can use
            indexes['chrom_gbwt'] = indexes['gbwt'] * len(inputGraphFileIDs)
        gcsa_job = gcsa_root_job.addChildJobFn(run_gcsa_prep, context, inputGraphFileIDs,
                                               graph_names, index_name, 
                                               indexes.get('chrom_gbwt', []) if gbwt_prune else [],
                                               node_mapping_id,
                                               remove_paths=dont_restore_paths,
                                               cores=context.config.misc_cores,
                                               memory=context.config.misc_mem,
                                               disk=context.config.misc_disk)
        indexes['gcsa'] = gcsa_job.rv(0)
        indexes['lcp'] = gcsa_job.rv(1)
    
    if len(inputGraphFileIDs) > 1 and 'id_ranges' in wanted:
        # Also we need an id ranges file in parallel with everything else
        indexes['id_ranges'] = child_job.addChildJobFn(run_id_ranges, context, inputGraphFileIDs,
                                                       graph_names, index_name, chroms,
                                                       cores=context.config.misc_cores,
                                                       memory=context.config.misc_mem,
                                                       disk=context.config.misc_disk).rv()
                                                 
    if 'snarls' in wanted:
        # Also we need a snarl index in parallel with everything else
        indexes['snarls'] = child_job.addChildJobFn(run_snarl_indexing, context, inputGraphFileIDs,
                                                    graph_names, index_name,
                                                    cores=context.config.snarl_index_cores,
                                                    memory=context.config.snarl_index_mem,
                                                    disk=context.config.snarl_index_disk).rv()
                                                    
    if 'trivial_snarls' in wanted:
        # Also we need a snarl index with trivial snarls in parallel with everything else.
        # Make sure to save the job so things can wait on it.
        trivial_snarls_job = child_job.addChildJobFn(run_snarl_indexing, context, inputGraphFileIDs,
                                                     graph_names, index_name, include_trivial=True,
                                                     cores=context.config.snarl_index_cores,
                                                     memory=context.config.snarl_index_mem,
                                                     disk=context.config.snarl_index_disk)
                                                     
        indexes['trivial_snarls'] = trivial_snarls_job.rv()
                                                            
    if 'distance' in wanted:
        # We need a distance index, based on the XG and the trivial snarls, which we know are being computed.
        # Run it after our XG
        distance_job = xg_root_job.addFollowOnJobFn(run_distance_indexing, context, indexes['xg'],
                                                    indexes['trivial_snarls'], index_name,
                                                    cores=context.config.distance_index_cores,
                                                    memory=context.config.distance_index_mem,
                                                    disk=context.config.distance_index_disk)
        # Make sure it waits for trivial snarls
        trivial_snarls_job.addFollowOn(distance_job)
        
        indexes['distance'] = distance_job.rv()
        
    if 'minimizer' in wanted and 'gbwt' in indexes:
        # We need a minimizer index, based on the GBWT (either provided or
        # computed) and the XG (which we know is being computed).
        
        # If there's no GBWT available, we can't compute a minimizer index.
        
        # Run it after our XG.
        # We know that, if the GBWT is being computed, it also happens under the XG job.
        # TODO: change that.
        minimizer_job = xg_root_job.addFollowOnJobFn(run_minimizer_indexing, context, indexes['xg'],
                                                     indexes['gbwt'], index_name,
                                                     cores=context.config.minimizer_index_cores,
                                                     memory=context.config.minimizer_index_mem,
                                                     disk=context.config.minimizer_index_disk)
                                                     
        indexes['minimizer'] = minimizer_job.rv()
    

    if bwa_fasta_id:
        # We need to index a reference FASTA for BWA
        indexes['bwa'] = child_job.addChildJobFn(run_bwa_index, context, bwa_fasta_id,
                                                 cores=context.config.bwa_index_cores, memory=context.config.bwa_index_mem,
                                                 disk=context.config.bwa_index_disk).rv()

    if 'alt-gam' in wanted:
        alt_extract_job = child_job.addChildJobFn(run_alt_path_extraction, context, inputGraphFileIDs,
                                                  graph_names, None,
                                                  cores=context.config.chunk_cores,
                                                  memory=context.config.chunk_mem,
                                                  disk=context.config.chunk_disk)
        
        indexes['alt-gam'] = alt_extract_job.addFollowOnJobFn(run_gam_indexing, context, alt_extract_job.rv(),
                                                              index_name,
                                                              cores=context.config.snarl_index_cores,
                                                              memory=context.config.snarl_index_mem,
                                                              disk=context.config.snarl_index_disk).rv()
    
    return indexes

def index_main(context, options):
    """
    Wrapper for vg indexing. 
    """

    # check some options
    validate_index_options(options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            inputGraphFileIDs = []
            for graph in options.graphs:
                inputGraphFileIDs.append(importer.load(graph))
            inputPhasingVCFFileIDs = []
            inputPhasingTBIFileIDs = []
            for vcf in options.vcf_phasing:
                inputPhasingVCFFileIDs.append(importer.load(vcf))
                inputPhasingTBIFileIDs.append(importer.load(vcf + '.tbi', wait_on = inputPhasingVCFFileIDs[-1]))
            # like construct, if there's one vcf and many graphs, we apply the vcf everywhere
            if len(inputPhasingTBIFileIDs) == 1 and len(options.graphs) > 1:
                inputPhasingVCFFileIDs = [inputPhasingVCFFileIDs[0]] * len(options.graphs)
                inputPhasingTBIFileIDs = [inputPhasingTBIFileIDs[0]] * len(options.graphs)
            inputGBWTID = None
            if options.gbwt_input:
                inputGBWTID = importer.load(options.gbwt_input)
            inputNodeMappingID = None
            if options.node_mapping:
                inputNodeMappingID = importer.load(options.node_mapping)
            inputBWAFastaID = None
            if options.bwa_index_fasta:
                inputBWAFastaID = importer.load(options.bwa_index_fasta)
            
            # Handy to have meaningful filenames throughout, so we remember
            # the input graph names
            graph_names = [os.path.basename(i) for i in options.graphs]

            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, context, importer.resolve(inputGraphFileIDs),
                                     graph_names, options.index_name, options.chroms,
                                     vcf_phasing_file_ids=importer.resolve(inputPhasingVCFFileIDs),
                                     tbi_phasing_file_ids=importer.resolve(inputPhasingTBIFileIDs),
                                     gbwt_id=importer.resolve(inputGBWTID),
                                     node_mapping_id=importer.resolve(inputNodeMappingID),
                                     bwa_fasta_id=importer.resolve(inputBWAFastaID),
                                     wanted=set(options.indexes), gbwt_prune=options.gbwt_prune,
                                     gbwt_regions=options.vcf_phasing_regions,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            index_key_and_id = toil.start(init_job)
        else:
            index_key_and_id = toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
