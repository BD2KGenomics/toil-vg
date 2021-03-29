#!/usr/bin/env python
"""
vg_construct.py: construct a graph from a vcf and fasta

"""

import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, struct, socket, threading
import string
import getpass
import pdb
import logging
import os

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_index import run_xg_indexing, run_indexing, run_bwa_index, index_parse_args, index_toggle_parse_args, validate_shared_index_options
from toil_vg.vg_msga import run_msga, msga_parse_args
logger = logging.getLogger(__name__)

# from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/IlluminaPlatinumGenomes-user-guide.pdf
CEPH_SAMPLES="NA12889 NA12890 NA12891 NA12892 NA12877 NA12878 NA12879 NA12880 NA12881 NA12882 NA12883 NA12884 NA12885 NA12886 NA12887 NA12888 NA12893".split()

def construct_subparser(parser):
    """
    Create a subparser for construction.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    parser.add_argument("--fasta", default=[], type=make_url, nargs='+',
                        help="Reference sequence in fasta or fasta.gz (single fasta or 1/region in same order as --regions)")
    parser.add_argument("--vcf", default=[], nargs='+',
                        help="Variants to make graph from (single vcf or 1/region IN SAME ORDER AS REGIONS (as passed in by --regions"
                        " or scanned in from --fasta_regions or --regions_file)).  "
                        "VCFs separated by commas will be merged together and treated as one using bcftools merge "
                        "(so must be over distinct sample sets).")
    parser.add_argument("--regions", default=[], nargs='+',
                        help="1-based inclusive VCF coordinates in the form of SEQ or SEQ:START-END")
    parser.add_argument("--regions_file", default=None,
                        help="List of regions (replaces --regions). Only first column of each line considered (so .fai acceptable)")
    parser.add_argument("--fasta_regions", action="store_true",
                        help="Infer regions from fasta file.  If multiple vcfs specified, any regions found that are not in --regions will be added without variants (useful for decoy sequences)")
    parser.add_argument("--regions_regex", default=[], nargs='+',
                        help="Ignore sequence names not fully matching (union of) given regexes when using --fasta_regions or --regions_file"
                        " (ex: --regions_regex \'chr[1-9,M,X,Y,EBV][0-9]{0,1}\' \'chr.*decoy\' to keep only chroms and decoys from hs38d1)")
    parser.add_argument("--alt_regions_bed", type=make_url,
                        help="BED file mapping alt regions (cols 1-3) to sequence names (col 4) from the FASTA. "
                        "Alt regions will be aligned to the graph using vg msga")
    parser.add_argument("--coalesce_regions", type=make_url,
                        help="File of tab-separated sets of sequence names for batching construction jobs, one per line.")
    parser.add_argument("--max_node_size", type=int, default=32,
                        help="Maximum node length")
    parser.add_argument("--alt_paths", action="store_true",
                        help="Save paths for alts with variant ID")
    parser.add_argument("--flat_alts", action="store_true",
                        help="flat alts")
    parser.add_argument("--handle_svs", action="store_true",
                        help="pass --handle-sv to vg construct to parse symbolic SV alts")
    parser.add_argument("--construct_cores", type=int,
                        help="Number of cores for vg construct")
    parser.add_argument("--out_name", default='graph',
                        help="Name used for output graphs and indexes")
    parser.add_argument("--merge_graphs", action="store_true",
                        help="Merge all regions into one graph")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize the graphs")
    parser.add_argument("--validate", action="store_true",
                        help="Run vg validate on constructed graphs")
    # useful for mixing, say, UCSC references with 1KG VCFs
    parser.add_argument("--add_chr_prefix", action="store_true",
                        help="add \"chr\" prefix to chromosome names if not already present")
    parser.add_argument("--remove_chr_prefix", action="store_true",
                        help="remove \"chr\" prefix from chromosome names")
    parser.add_argument("--keep_vcfs", action="store_true",
                        help="write the VCFs created to make the filtered and control graphs to the output store") 
    # Toggles for the different types of graph(s) that can be made.  Indexing and above options
    # will be applied to each one.  The output names will be prefixed with out_name. 
    parser.add_argument("--primary", action="store_true",
                        help="Make the primary graph (no variants) using just the FASTA")
    parser.add_argument("--pangenome", action="store_true",
                        help="Make the pangenome graph using the input VCF(s)")
    parser.add_argument("--pos_control", type=str,
                        help="Make a positive control (ref path plus sample variants) using this sample")
    parser.add_argument("--neg_control", type=str,
                        help="Make a negative control (exclude all sample variants) using this sample")
    parser.add_argument("--sample_graph", type=str,
                        help="Make a sample graph (only contains sample haplotypes) using this sample.  Only "
                        " phased variants will be included.  Will also make a _withref version that includes reference")
    parser.add_argument("--haplo_sample", type=str,
                        help="Make two haplotype thread graphs (for simulating from) for this sample.  Phasing"
                        " information required in the input vcf.")    
    parser.add_argument("--filter_ceph", action="store_true",
                        help="Make a graph where all variants private to the CEPH pedigree, which includes "
                        "NA12878 are excluded")
    parser.add_argument("--filter_samples", nargs='+',
                        help="Make a graph where all variants private to the the listed smaples are excluded")
    parser.add_argument("--min_af", type=float, default=[], nargs='+',
                        help="Create a graph including only variants with given minium allele frequency."
                        " If multiple frequencies given, a graph will be made for each one")
    parser.add_argument("--bwa_reference", type=make_url,
                        help="Make a BWA reference (set of indexes) from the given FASTA (not the --fasta FASTAs).")

    parser.add_argument("--pre_min_af", type=float, default=None,
                        help="Run minimum allele frequency filter as preprocessing step on each input VCF.  "
                        "Unlike --min_af, this will be applied before merging and any control graph construction")
    parser.add_argument("--mask_ambiguous", action="store_true",
                        help="Convert IUPAC ambiguous characters in FASTA to Ns")
    
    # Add common indexing options shared with vg_index
    index_toggle_parse_args(parser)
    index_parse_args(parser)

    # Add common msga options shared with vg_msga
    msga_parse_args(parser)

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def re_fullmatch(regex, string, flags=0):
    """Emulate python-3.4 re.fullmatch().
    https://stackoverflow.com/questions/30212413/backport-python-3-4s-regular-expression-fullmatch-to-python-2
    """
    return re.match("(?:" + regex + r")\Z", string, flags=flags)
    
def validate_construct_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(options.regions or options.fasta_regions or options.regions_file,
            '--regions or --fasta_regions required')
    require(not options.regions_file or not (options.fasta_regions or options.regions),
            '--regions_file cannot be used with --regions or --fasta_regions')
    require(not options.regions_regex or (options.fasta_regions or options.regions_file),
            '--regions_regex can only be used with --fasta_regions or --regions_file')
    require(not options.add_chr_prefix or not options.remove_chr_prefix,
            '--add_chr_prefix cannot be used with --remove_chr_prefix')
    require(options.vcf == [] or len(options.vcf) == 1 or not options.regions or
            len(options.vcf) <= len(options.regions),
            'if many vcfs specified, cannot have more vcfs than --regions')
    require(len(options.fasta) == 1 or len(options.fasta) == len(options.regions),
            'if many fastas specified, must be same number as --regions')
    require(len(options.fasta) == 1 or not options.fasta_regions,
            '--fasta_regions currently only works when single fasta specified with --fasta')
    require(len(options.fasta) > 0 or options.bwa_reference,
            'either --fasta or --bwa_reference must be set to give something to construct')
    require('gbwt' not in options.indexes or 'xg' in options.indexes,
            '--xg_index required with --gbwt_index')
    # TODO: It seems like some of this code is designed to run multiple regions
    # in parallel, but the indexing code always indexes them together.
    require('gbwt' not in options.indexes or (not options.pangenome and not options.pos_control and
        not options.neg_control and not options.sample_graph and not options.haplo_sample and
        not options.min_af) or len(options.vcf) >= 1,
            '--gbwt_index with any graph other than --primary requires --vcf')
    require(not options.sample_graph or (options.regions or options.fasta_regions),
            '--regions or --fasta_regions required with --sample_graph')
    require(options.primary or options.pangenome or options.pos_control or options.neg_control or
            options.sample_graph or options.haplo_sample or options.filter_ceph or options.filter_samples or
            options.min_af or options.bwa_reference,
            'At least one kind of graph or reference must be specified for construction')
    require(not options.vcf or options.pangenome or options.pos_control or options.neg_control or
            options.sample_graph or options.haplo_sample or options.filter_ceph or options.filter_samples or
            options.min_af,
            'At least one kind of non-primary graph must be specified for construction with --vcf')
    require(options.vcf or not (options.pangenome or options.pos_control or options.neg_control or
                                options.sample_graph or options.haplo_sample or options.filter_ceph or
                                options.filter_samples or options.min_af),
            '--vcf required for construction of non-primary graph')
    # TODO: support new, more general CLI properly
    require(options.pos_control is None or options.neg_control is None or
            options.pos_control == options.neg_control,
            '--pos_control_sample and --neg_control_sample must be the same')
    require(not options.haplo_sample or not options.sample_graph or
            (options.haplo_sample == options.sample_graph),
            '--haplo_sample and --sample_graph must be the same')

    validate_shared_index_options(options)

def chr_name_map(to_ucsc, max_chrom=22):
    """
    Return a name map for chromosome name conversion in dict and string format,
    and a TSV string version of the same.
    
    Will contain mappings for chromosomes 1 to max_chrom, inclusive.
    """
    name_map = {}
    name_str = ''
    # TODO:  Should we do something with MT <==> chrM ?
    for i in list(range(1, max_chrom + 1)) + ['X', 'Y']:
        if to_ucsc:
            name_str += '{}\tchr{}\n'.format(i, i)
            name_map[str(i)] = 'chr{}'.format(i)
        else:
            name_str += 'chr{}\t{}\n'.format(i, i)
            name_map['chr{}'.format(i)] = str(i)
    return name_map, name_str    

def run_merge_all_vcfs(job, context, vcf_file_ids_list, vcf_names_list, tbi_file_ids_list):
    """
    takes a lists of lists of input vcfs.  make child merge job for each list
    """
    out_vcf_ids_list = []
    out_names_list = []
    out_tbi_ids_list = []
    for vcf_file_ids, vcf_names, tbi_file_ids in zip(vcf_file_ids_list, vcf_names_list, tbi_file_ids_list):
        if len(vcf_file_ids) > 1:
            merge_job = job.addChildJobFn(run_merge_vcfs, context, vcf_file_ids, vcf_names, tbi_file_ids,
                                          cores=context.config.construct_cores,
                                          memory=context.config.construct_mem,
                                          disk=context.config.construct_disk)
            out_vcf_ids_list.append(merge_job.rv(0))
            out_names_list.append(merge_job.rv(1))
            out_tbi_ids_list.append(merge_job.rv(2))
        else:
            out_vcf_ids_list.append(vcf_file_ids[0])
            out_names_list.append(vcf_names[0])
            out_tbi_ids_list.append(tbi_file_ids[0])
    return out_vcf_ids_list, out_names_list, out_tbi_ids_list
        
def run_merge_vcfs(job, context, vcf_file_ids, vcf_names, tbi_file_ids):
    """
    run bctools merge on a list of vcfs and return just one.  note that 
    bcftools merge expectes non-overlapping sample sets
    """
    assert len(vcf_file_ids) == len(vcf_names) == len(tbi_file_ids)
    if len(vcf_file_ids) == 1:
        return vcf_file_ids[0], vcf_names[0], tbi_file_ids[0]

    work_dir = job.fileStore.getLocalTempDir()

    names = []
    for vcf_id, vcf_name, tbi_id in zip(vcf_file_ids, vcf_names, tbi_file_ids):
        job.fileStore.readGlobalFile(vcf_id, os.path.join(work_dir, vcf_name))
        job.fileStore.readGlobalFile(tbi_id, os.path.join(work_dir, vcf_name) + '.tbi')
        names.append(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
    if len(names) != len(set(names)):
        raise RuntimeError('vcf merging expects unique filenames')

    merged_name = '_'.join(names) + '.vcf.gz'
    with open(os.path.join(work_dir, merged_name), 'wb') as merged_file:
        cmd = [['bcftools', 'merge', '--missing-to-ref', '--force-samples'] + vcf_names]
        # phase the ref/ref calls added by --missing-to-ref
        cmd.append(['sed', '-e', 's/0\/0/0\|0/g'])
        cmd.append(['bcftools', 'view', '-', '--output-type', 'z'])
        context.runner.call(job, cmd, work_dir = work_dir, outfile = merged_file)
    context.runner.call(job, ['tabix', '--preset', 'vcf', merged_name], work_dir = work_dir)

    return (context.write_intermediate_file(job, os.path.join(work_dir, merged_name)),
            merged_name,
            context.write_intermediate_file(job, os.path.join(work_dir, merged_name) + '.tbi'))

    
def run_unzip_fasta(job, context, fasta_id, fasta_name):
    """
    vg construct doesn't work with zipped fasta, so we run this on input fasta that end in .gz
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file, mutable=True)
    context.runner.call(job, ['bgzip', '-d', os.path.basename(fasta_file)], work_dir=work_dir)

    return context.write_intermediate_file(job, fasta_file[:-3])

def run_mask_ambiguous(job, context, fasta_id, fasta_name):
    """
    Replace IUPAC characters (of any case) with Ns.  That's how they end up in the XG anyway,
    and it will prevent some errors in vg construct if the VCF has N's but the fasta has something else. 
    (todo: would need to apply same thing to the VCF to be more robust, but will hold off until having
    a use case)
    """

    work_dir = job.fileStore.getLocalTempDir()

    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    mask_file = os.path.splitext(fasta_file)[0] + '-mask.fa'
    job.fileStore.readGlobalFile(fasta_id, fasta_file, mutable=True)

    fa_mask_cmd = ['awk',  'BEGIN{FS=\" \"}{if(!/>/){ gsub(/[YRWSKMDVHBXyrwskmdvhbx]/,"N"); print }else{print $1}}',
                   os.path.basename(fasta_file)]
    with open(mask_file, 'wb') as mf:
        context.runner.call(job, fa_mask_cmd, outfile=mf, work_dir=work_dir)

    return context.write_intermediate_file(job, mask_file), os.path.basename(mask_file)

def run_scan_fasta_sequence_names(job, context, fasta_id, fasta_name, regions = None, regions_regex = None):
    """
    scrape regions out of the (uncompressed) fasta, appending them to given regions list if provided
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file)
    
    # reluctant to use slow python library, so just running grep instead
    cmd = ['grep', '>', os.path.basename(fasta_file)]
    grep_output = context.runner.call(job, cmd, work_dir = work_dir,
                                      check_output = True, tool_name='bgzip')

    # just taking first whitespace-separated token.  that's what corresponds to hs37d5 vcf
    seq_names = [] if not regions else regions
    for line in grep_output.decode().split('\n'):
        if len(line) > 1:
            name = line.split()[0]
            if name.startswith('>') and (not regions or name[1:] not in regions) and \
               (not regions_regex or re_fullmatch(regions_regex, name[1:])):
                seq_names.append(name[1:])
    
    return seq_names

def run_scan_regions_file(job, context, regions_id, regions_regex = None):
    """
    Read a list of regions
    """
    work_dir = job.fileStore.getLocalTempDir()
    regions_path = os.path.join(work_dir, 'regions.tsv')
    job.fileStore.readGlobalFile(regions_id, regions_path)
    out_regions = []
    with open(regions_path) as regions_file:
        for line in regions_file:
            region_name = line.strip().split()[0]
            if len(region_name) > 0 and (not regions_regex or re_fullmatch(regions_regex, region_name)):
                out_regions.append(region_name)
    return out_regions

def run_fix_chrom_names(job, context, to_ucsc, regions, fasta_ids, fasta_names,
                        vcf_ids_list, vcf_names_list, tbi_ids_list, alt_regions_id):
    """
    Apply name mappings to regions list, fasta files and vcf files.  if to_ucsc is true we convert
    1 -> chr1 etc.  otherwise, we go the other way.  
    """

    work_dir = job.fileStore.getLocalTempDir()
    
    # How many chromosomes should we generate name mappings for?
    # No more than there are regions certainly.
    # But also never less than the 22 we expect in humans.
    max_chrom = max(22, len(regions))

    name_map, name_str = chr_name_map(to_ucsc, max_chrom)
    out_regions = []
    something_to_rename = False

    # Map the regions.
    # Some regions may map to the same region as a previous region if e.g. we
    # fetched region names out of the FASTA. So deduplicate here.
    done_regions = set()
    for region in regions:
        region_name = region.split(':')[0]
        if region_name in name_map:
            something_to_rename = True
            renamed_region = name_map[region_name] + region[len(region_name):]
            if renamed_region not in done_regions:
                out_regions.append(renamed_region)
                done_regions.add(renamed_region)
        else:
            something_to_rename = something_to_rename or region_name in list(name_map.values())
            if region not in done_regions:
                out_regions.append(region)
                done_regions.add(region)
        
    # map the vcf
    out_vcf_ids = []
    out_vcf_names = []
    out_tbi_ids = []
    if something_to_rename:
        # make our name mapping file
        name_map_path = os.path.join(work_dir, 'name_map.tsv')
        with open(name_map_path, 'w') as name_map_file:
            name_map_file.write(name_str)
        name_map_id = context.write_intermediate_file(job, name_map_path)
        
        for vcf_ids, vcf_names, tbi_ids in zip(vcf_ids_list, vcf_names_list, tbi_ids_list):
            out_vcf_ids.append([])
            out_vcf_names.append([])
            out_tbi_ids.append([])
            for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
                vcf_rename_job = job.addChildJobFn(run_fix_vcf_chrom_names, context, vcf_id, vcf_name, tbi_id, name_map_id,
                                                   cores=context.config.construct_cores,
                                                   memory=context.config.construct_mem,
                                                   disk=context.config.construct_disk)
                out_vcf_ids[-1].append(vcf_rename_job.rv(0))
                out_vcf_names[-1].append(vcf_rename_job.rv(1))
                out_tbi_ids[-1].append(vcf_rename_job.rv(2))
    else:
        out_vcf_ids = vcf_ids_list
        out_vcf_names = vcf_names_list
        out_tbi_ids = tbi_ids_list

    # map the fasta
    out_fasta_ids = []    
    out_fasta_names = []
    if something_to_rename:
        for fasta_id, fasta_name in zip(fasta_ids, fasta_names):
            assert not fasta_name.endswith('.gz')
            in_fasta_name = os.path.basename(fasta_name)
            job.fileStore.readGlobalFile(fasta_id, os.path.join(work_dir, in_fasta_name))
            out_fasta_name = os.path.splitext(fasta_name)[0] + '-renamed' + os.path.splitext(fasta_name)[1]
            with open(os.path.join(work_dir, out_fasta_name), 'w') as out_fasta_file, \
                 open(os.path.join(work_dir, in_fasta_name)) as in_fasta_file:
                # TODO: is this too slow in python?
                for line in in_fasta_file:
                    if line.startswith('>'):
                        region_name = line[1:].split()[0]
                        if region_name in name_map:
                            out_fasta_file.write('>{}\n'.format(name_map[region_name]))
                        else:
                            out_fasta_file.write(line)
                    else:
                        out_fasta_file.write(line)
            out_fasta_ids.append(context.write_intermediate_file(job, os.path.join(work_dir, out_fasta_name)))
            out_fasta_names.append(out_fasta_name)
    else:
        out_fasta_ids = fasta_ids
        out_fasta_names = fasta_names

    # map the alt regions
    if alt_regions_id:
        alt_regions_path = os.path.join(work_dir, 'alt-regions.bed')
        alt_regions_out_path = os.path.join(work_dir, 'alt-regions-fix.bed')
        job.fileStore.readGlobalFile(alt_regions_id, alt_regions_path)
        with open(alt_regions_path) as in_regions, open(alt_regions_out_path, 'w') as out_alt_regions:
            for line in in_regions:
                toks = line.strip().split('\t')
                if len(toks) >= 4 and toks[0] != '#':
                    if toks[0] in name_map:
                            out_alt_regions.write('{}\t{}\t{}\t{}\n'.format(name_map[toks[0]], toks[1], toks[2], toks[3]))
                    else:
                        out_alt_regions.write(line)
        out_alt_regions_id = context.write_intermediate_file(job, alt_regions_out_path)
    else:
        out_alt_regions_id = None

    return out_regions, out_fasta_ids, out_fasta_names, out_vcf_ids, out_vcf_names, out_tbi_ids, out_alt_regions_id

def run_fix_vcf_chrom_names(job, context, vcf_id, vcf_name, tbi_id, name_file_id):
    """
    use bcftools annotate to rename chromosomes in a vcf
    """
    work_dir = job.fileStore.getLocalTempDir()
    name_map_path = os.path.join(work_dir, 'name_map.tsv')
    job.fileStore.readGlobalFile(name_file_id, name_map_path)

    assert vcf_name.endswith('.vcf.gz')
    in_vcf_name = os.path.basename(vcf_name)
    job.fileStore.readGlobalFile(vcf_id, os.path.join(work_dir, in_vcf_name))
    job.fileStore.readGlobalFile(tbi_id, os.path.join(work_dir, in_vcf_name + '.tbi'))
    out_vcf_name = in_vcf_name[:-7] + '-renamed.vcf.gz'
    context.runner.call(job, ['bcftools', 'annotate', '--rename-chrs', os.path.basename(name_map_path),
                              '--output-type', 'z', '--output', out_vcf_name, os.path.basename(in_vcf_name)],
                        work_dir = work_dir)
    context.runner.call(job, ['tabix', '--force', '--preset', 'vcf', out_vcf_name], work_dir = work_dir)
    return (context.write_intermediate_file(job, os.path.join(work_dir, out_vcf_name)),
            out_vcf_name,
            context.write_intermediate_file(job, os.path.join(work_dir, out_vcf_name + '.tbi')))

def run_subtract_alt_regions(job, context, alt_regions_id, regions):
    """
    make sure that alt contigs don't wind up in our regions names, as we want them
    to get aligned into chromosomes rather than form their own components
    """
    work_dir = job.fileStore.getLocalTempDir()
    alt_regions_path = os.path.join(work_dir, 'alt-regions.bed')
    job.fileStore.readGlobalFile(alt_regions_id, alt_regions_path)
    alt_regions = set()
    with open(alt_regions_path) as in_regions:
        for line in in_regions:
            toks = line.strip().split('\t')
            if len(toks) >= 4 and toks[0] != '#':
                alt_regions.add(toks[3])
                
    return [region for region in regions if region not in alt_regions], list(alt_regions)
    
def run_read_coalesce_list(job, context, coalesce_regions_id):
    """
    Read the given input file.
    Produce a list of sets of region/contig names to coalesce into single jobs.
    Treats the input file as tab-separated, one set per line.
    """
    
    work_dir = job.fileStore.getLocalTempDir()
    coalesce_regions_path = os.path.join(work_dir, 'coalesce-regions.tsv')
    job.fileStore.readGlobalFile(coalesce_regions_id, coalesce_regions_path)
    coalesce_regions = []
    with open(coalesce_regions_path) as in_regions:
        for line in in_regions:
            toks = line.strip().split('\t')
            if len(toks) > 0 and toks[0] != '' and toks[0][0] != '#':
                coalesce_regions.append(set(toks))
    
    return coalesce_regions
                
def run_generate_input_vcfs(job, context, vcf_ids, vcf_names, tbi_ids,
                            regions, output_name,
                            do_primary = False,
                            do_pan = False,
                            pos_control_sample = None,
                            neg_control_sample = None,
                            sample_graph = None,
                            haplo_sample = None,
                            filter_samples = [],
                            min_afs = [],
                            vcf_subdir = None):
    """
    Preprocessing step to make a bunch of vcfs if wanted:
    - positive control
    - negative control
    - family filter
    - primary
    - thresholded by a given minimum allele frequency
    returns a dictionary of name -> (vcf_id, vcf_name, tbi_id, merge_name, region_names) tuples
    where name can be used to, ex, tell the controls apart
    if vcf_subdir is specified, the various created vcfs will be stored in a subfolder of that
    name in the output store.  if it's not specified, then these intermediate vcfs will not be saved
    """

    output = dict()

    # primary graph, containing just the input fasta
    if do_primary:
        if regions:
            primary_region_names = ['primary'  + '_' + c.replace(':','-') for c in regions]
        else:
            primary_region_names = None

        primary_output_name = 'primary.vg' if '_' not in output_name else 'primary' + output_name[output_name.find('_')+1:]
        output['primary'] = [[], [], [], primary_output_name, primary_region_names]

    # a straight-up pangenome graph from the input vcf
    if do_pan:
        output[output_name] = [vcf_ids, vcf_names, tbi_ids, output_name,
                               [output_name + '_' + c.replace(':','-') for c in regions] if regions else None]

    # our positive control consists of the reference path and any variant in the sample
    if pos_control_sample or neg_control_sample:
        control_sample = pos_control_sample if pos_control_sample else neg_control_sample
        assert not neg_control_sample or neg_control_sample == control_sample
        assert not pos_control_sample or pos_control_sample == control_sample
        pos_control_vcf_ids, pos_control_tbi_ids = [], []
        neg_control_vcf_ids, neg_control_tbi_ids = [], []
        pos_control_vcf_names, neg_control_vcf_names = [], []
        
        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id,
                                              control_sample,
                                              pos_only = not neg_control_sample,
                                              vcf_subdir = vcf_subdir,
                                              cores=context.config.construct_cores,
                                              memory=context.config.construct_mem,
                                              disk=context.config.construct_disk)
            pos_control_vcf_ids.append(make_controls.rv(0))
            pos_control_tbi_ids.append(make_controls.rv(1))
            neg_control_vcf_ids.append(make_controls.rv(2))
            neg_control_tbi_ids.append(make_controls.rv(3))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            pos_control_vcf_names.append('{}_{}.vcf.gz'.format(vcf_base, control_sample))
            neg_control_vcf_names.append('{}_minus_{}.vcf.gz'.format(vcf_base, control_sample))
        if regions:
            pos_region_names = [output_name + '_{}'.format(control_sample)  + '_' + c.replace(':','-') for c in regions]
            neg_region_names = [output_name + '_minus_{}'.format(control_sample) + '_' + c.replace(':','-')  for c in regions]
        else:
            pos_region_names = None
            neg_region_names = None
        pos_output_name = remove_ext(output_name, '.vg') + '_{}.vg'.format(control_sample)
        neg_output_name = remove_ext(output_name, '.vg') + '_minus_{}.vg'.format(control_sample)

        if pos_control_sample:
            output['pos-control'] = [pos_control_vcf_ids, pos_control_vcf_names, pos_control_tbi_ids,
                                     pos_output_name, pos_region_names]

        if neg_control_sample:
            output['neg-control'] = [neg_control_vcf_ids, neg_control_vcf_names, neg_control_tbi_ids,
                                     neg_output_name, neg_region_names]

    # For our sample graph, we're going to need to start by making someing like the positive control, but
    # filtering for phased variants.  Note that making the actual graphs from these vcfs is a two step
    # process, where a graph is constructed then haplotypes extracted.
    if sample_graph:
        sample_graph_vcf_ids, sample_graph_tbi_ids = [], []
        sample_graph_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            make_sample = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, sample_graph,
                                            pos_only = True,
                                            vcf_subdir = vcf_subdir,
                                            cores=context.config.construct_cores,
                                            memory=context.config.construct_mem,
                                            disk=context.config.construct_disk)
            sample_graph_vcf_ids.append(make_sample.rv(0))
            sample_graph_tbi_ids.append(make_sample.rv(1))
            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            sample_graph_vcf_names.append('{}_{}_sample_withref.vcf.gz'.format(vcf_base, sample_graph))
        if regions:
            sample_graph_region_names = [output_name + '_{}_sample_withref'.format(sample_graph)  + '_' + c.replace(':','-') for c in regions]
        else:
            sample_graph_region_names = None
        sample_graph_output_name = remove_ext(output_name, '.vg') + '_{}_sample_withref.vg'.format(sample_graph)
        
        output['sample-graph'] = [sample_graph_vcf_ids, sample_graph_vcf_names, sample_graph_tbi_ids,
                                  sample_graph_output_name, sample_graph_region_names]


    # we want a vcf to make a gbwt out of for making haplo graphs
    if haplo_sample:
        hap_control_vcf_ids, hap_control_tbi_ids = [], []
        hap_control_vcf_names = []
        
        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, haplo_sample,
                                              pos_only = True, 
                                              vcf_subdir = vcf_subdir,
                                              cores=context.config.construct_cores,
                                              memory=context.config.construct_mem,
                                              disk=context.config.construct_disk)
            hap_control_vcf_ids.append(make_controls.rv(0))
            hap_control_tbi_ids.append(make_controls.rv(1))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            hap_control_vcf_names.append('{}_{}_haplo.vcf.gz'.format(vcf_base, haplo_sample))
            
        if regions:
            hap_region_names = [output_name + '_{}_haplo'.format(haplo_sample)  + '_' + c.replace(':','-') for c in regions]
        else:
            hap_region_names = None
        hap_output_name = remove_ext(output_name, '.vg') + '_{}_haplo.vg'.format(haplo_sample)
        
        output['haplo'] = [hap_control_vcf_ids, hap_control_vcf_names, hap_control_tbi_ids,
                           hap_output_name, hap_region_names]
        
    # our family filter
    if filter_samples:
        filter_vcf_ids, filter_tbi_ids = [], []
        filter_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            filter_job = job.addChildJobFn(run_filter_vcf_samples, context, vcf_id, vcf_name, tbi_id,
                                           filter_samples,
                                           vcf_subdir = vcf_subdir,
                                           cores=context.config.construct_cores,
                                           memory=context.config.construct_mem,
                                           disk=context.config.construct_disk)

            filter_vcf_ids.append(filter_job.rv(0))
            filter_tbi_ids.append(filter_job.rv(1))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            filter_vcf_names.append('{}_filter.vcf.gz'.format(vcf_base))
        if regions:
            filter_region_names = [output_name + '_filter'  + '_' + c.replace(':','-') for c in regions]
        else:
            filter_region_names = None
        filter_output_name = remove_ext(output_name, '.vg') + '_filter.vg'

        output['filter'] = [filter_vcf_ids, filter_vcf_names, filter_tbi_ids,
                            filter_output_name, filter_region_names]

    # and one for each minimum allele frequency filter
    for min_af in min_afs:
        af_vcf_ids, af_tbi_ids = [], []
        af_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            af_job = job.addChildJobFn(run_min_allele_filter_vcf_samples, context, vcf_id, vcf_name, tbi_id,
                                       min_af,
                                       vcf_subdir = vcf_subdir,
                                       cores=context.config.construct_cores,
                                       memory=context.config.construct_mem,
                                       disk=context.config.construct_disk)

            af_vcf_ids.append(af_job.rv(0))
            af_tbi_ids.append(af_job.rv(1))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            af_vcf_names.append('{}_minaf_{}.vcf.gz'.format(vcf_base, min_af))
        if regions:
            af_region_names = [output_name + '_minaf_{}'.format(min_af) + '_' + c.replace(':','-') for c in regions]
        else:
            af_region_names = None
        af_output_name = remove_ext(output_name, '.vg') + '_minaf_{}.vg'.format(min_af)

        output['minaf-{}'.format(min_af)] = [af_vcf_ids, af_vcf_names, af_tbi_ids,
                           af_output_name, af_region_names]

    # pad out vcf lists with nones so they are the same size as regions
    # since we allow fasta regions that dont have corresponding vcf
    # note: single vcf , multiple regions case not handled here as it's
    # treated below (the same vcf is given to each region)
    if regions and len(regions) > len(vcf_ids) and len(vcf_ids) != 1:
        padding = [None] * (len(regions) - len(vcf_ids))
        for key, val in list(output.items()):
            val[0] += padding
            val[1] += padding
            val[2] += padding

    return output

    
def run_construct_all(job, context, fasta_ids, fasta_names, vcf_inputs, 
                      max_node_size, alt_paths, flat_alts, handle_svs, regions,
                      merge_graphs = False, sort_ids = False, join_ids = False,
                      wanted_indexes = set(), 
                      haplo_extraction_sample = None, haplotypes = [0,1], gbwt_prune = False,
                      normalize = False, validate = False, alt_regions_id = None,
                      alt_regions = [], coalesce_regions = []):
    """ 
    construct many graphs in parallel, optionally doing indexing too. vcf_inputs
    is a list of tuples as created by run_generate_input_vcfs
    
    Returns a list of tuples of the form (vg_ids, vg_names, indexes), where
    indexes is the index dict from index type to file ID.
    
    If coalesce_regions is set (and alt_regions is not), it must be a list of
    sets of region names. Each set of region names will be run together as a
    single construction job, replacing the individual jobs for those regions.
    
    """

    output = []
    
    for name, (vcf_ids, vcf_names, tbi_ids, output_name, region_names) in list(vcf_inputs.items()):
        merge_output_name = output_name if merge_graphs or not regions or len(regions) < 2 else None
        output_name_base = remove_ext(output_name, '.vg')
        # special case that need thread indexes no matter what
        haplo_extraction = name in ['haplo', 'sample-graph']
        construct_job = job.addChildJobFn(run_construct_genome_graph, context, fasta_ids,
                                          fasta_names, vcf_ids, vcf_names, tbi_ids,
                                          max_node_size, ('gbwt' in wanted_indexes) or haplo_extraction or alt_paths,
                                          flat_alts, handle_svs, regions,
                                          region_names, sort_ids, join_ids, name, merge_output_name,
                                          normalize and name != 'haplo', validate, alt_regions_id,
                                          coalesce_regions=coalesce_regions)
                                          
        mapping_id = construct_job.rv('mapping')
        
        # Find the joined VG files, which always exist
        joined_vg_ids = construct_job.rv('joined')
        # And give them names
        joined_vg_names = construct_job.rv('joined_names')
        
        if merge_graphs or not regions or len(regions) < 2: 
            # Sometimes we will have a single VG file also
            single_vg_id = construct_job.rv('merged')
            single_vg_name = remove_ext(merge_output_name, '.vg') + '.vg'
        else:
            # But sometimes not
            single_vg_id = None
            single_vg_name = None
            
        # Now the graphs are ready

        if not regions:
            chroms = []
            gbwt_regions = []
        else:
            # We have regions specified to restrict to.
            
            # Get the chromosome names
            chroms = [p.split(':')[0] for p in regions]
            
            RealtimeLogger.info('Operating on chromosomes: %s', chroms)
            
            # Make sure we have no more than 1 region per chromosome.
            # Otherwise GBWT region restriction will mess things up.
            assert(len(chroms) == len(set(chroms)))
            
            # Get the regions that are restrictions smaller than a whole chromosome to hint the GBWT.
            # Otherwise running a small region of a big VCF means a very slow GBWT construction step.
            gbwt_regions = [p for p in regions if ':' in p]

        # strip nones out of vcf list            
        input_vcf_ids = []
        input_tbi_ids = []
        if haplo_extraction or ('gbwt' in wanted_indexes):
            for vcf_id, tbi_id in zip(vcf_ids, tbi_ids):
                if vcf_id and tbi_id:
                    input_vcf_ids.append(vcf_id)
                    input_tbi_ids.append(tbi_id)
                else:
                    assert vcf_id == None and tbi_id == None

        index_prev_job = construct_job
        if haplo_extraction:
            haplo_index_job = construct_job.addFollowOnJobFn(run_make_haplo_indexes, context,
                                                             input_vcf_ids, input_tbi_ids,
                                                             vcf_names, joined_vg_ids, joined_vg_names,
                                                             output_name_base, regions, haplo_extraction_sample,
                                                             intermediate=merge_graphs)
            haplo_xg_ids = haplo_index_job.rv(0)
            gbwt_ids = haplo_index_job.rv(1)
            
            if name == 'sample-graph':
                # ugly hack to distinguish the graphs with reference and our extracted sample graph
                sample_name_base = output_name_base.replace('_withref', '')
                sample_merge_output_name = merge_output_name.replace('_withref', '') if merge_output_name else None
                region_names = [r.replace('_withref', '') for r in region_names]
                # Extract out our real sample graph          
                sample_job = haplo_index_job.addFollowOnJobFn(run_make_sample_graphs, context,
                                                              joined_vg_ids, joined_vg_names,
                                                              haplo_xg_ids, sample_name_base, regions,
                                                              haplo_extraction_sample, gbwt_ids)
                # Put them back together again with a no-op join, producing
                # many graphs again and maybe a single merged graph with the
                # reference removed.
                join_job = sample_job.addFollowOnJobFn(run_join_graphs, context, sample_job.rv(),
                                                       False, region_names, name, sample_merge_output_name,
                                                       cores=context.config.construct_cores,
                                                       memory=context.config.xg_index_mem,
                                                       disk=context.config.xg_index_disk)

                # Want to keep a whole-genome withref xg index around for mapeval purposes
                if len(regions) > 1 and ('xg' in wanted_indexes):
                    wanted = set('xg')
                    construct_job.addFollowOnJobFn(run_indexing, context, joined_vg_ids,
                                                   joined_vg_names, output_name_base, chroms, [], [], 
                                                   wanted=wanted, coalesce_regions=coalesce_regions)
                
                # In the indexing step below, we want to index our haplo-extracted sample graph
                # So replace the withref graph IDs and names with these
                
                # Find the joined VG files, which always exist
                joined_vg_ids = join_job.rv('joined')
                # And give them names
                rename_job = join_job.addFollowOnJobFn(run_remove_withref, context, joined_vg_names,
                                                       cores=context.config.misc_cores,
                                                       memory=context.config.misc_mem,
                                                       disk=context.config.misc_disk)
                joined_vg_names = rename_job.rv()
                index_prev_job = rename_job
                
                if sample_merge_output_name:
                    # We expect a single output graph too
                    single_vg_id = join_job.rv('merged')
                    single_vg_name = remove_ext(sample_merge_output_name, '.vg') + '.vg'
                else:
                    # No single merged graph
                    single_vg_id = None
                    single_vg_name = None
                
                output_name_base = sample_name_base
            
            elif name == 'haplo':
                assert haplo_extraction_sample is not None            
                haplo_job = haplo_index_job.addFollowOnJobFn(run_make_haplo_graphs, context,
                                                             joined_vg_ids, joined_vg_names, haplo_xg_ids,
                                                             output_name_base, regions,
                                                             haplo_extraction_sample, haplotypes, gbwt_ids,
                                                             intermediate = merge_graphs)

                # we want an xg index from our thread graphs to pass to vg sim for each haplotype
                for haplotype in haplotypes:
                    haplo_xg_job = haplo_job.addFollowOnJobFn(run_xg_indexing, context, haplo_job.rv(haplotype),
                                                              joined_vg_names,
                                                              output_name_base + '_thread_{}'.format(haplotype),
                                                              include_alt_paths = 'xg_alts' in wanted_indexes,
                                                              cores=context.config.xg_index_cores,
                                                              memory=context.config.xg_index_mem,
                                                              disk=context.config.xg_index_disk)
                    
        # some indexes should never get built for haplo/sample graphs.
        # So work out what indexes to build.
        wanted = set(wanted_indexes)
        if name == 'haplo':
            wanted.discard('gcsa')
        if haplo_extraction:
            wanted.discard('snarls')
            wanted.discard('trivial_snarls')
            wanted.discard('gbwt')
        
        indexing_job = index_prev_job.addFollowOnJobFn(run_indexing, context, joined_vg_ids,
                                                       joined_vg_names, output_name_base, chroms,
                                                       input_vcf_ids if ('gbwt' in wanted) else [],
                                                       input_tbi_ids if ('gbwt' in wanted) else [],
                                                       node_mapping_id=mapping_id,
                                                       wanted=wanted,
                                                       gbwt_prune=gbwt_prune and 'gbwt' in wanted,
                                                       gbwt_regions=gbwt_regions,
                                                       dont_restore_paths=alt_regions,
                                                       coalesce_regions=coalesce_regions)
        indexes = indexing_job.rv()    

        output.append((joined_vg_ids, joined_vg_names, indexes))
    return output
    
def run_remove_withref(job, context, joined_vg_names):
    """
    Return the names in joined_vg_names with '_withref' removed from them.
    """
    
    return [n.replace('_withref', '') for n in joined_vg_names]
                

def run_construct_genome_graph(job, context, fasta_ids, fasta_names, vcf_ids, vcf_names, tbi_ids,
                               max_node_size, alt_paths, flat_alts, handle_svs, regions, region_names,
                               sort_ids, join_ids, name, merge_output_name, normalize, validate,
                               alt_regions_id, coalesce_regions=[]):
    """
    
    Construct graphs from one or more FASTA files and zero or more VCFs.
    
    If regions and region_names are set, constructs only for the specified
    regions, and constructs one graph per region (except if coalesce_regions is
    set; see below). Otherwise, constructs one graph overall on a single
    default region.
    
    If coalesce_regions is set (and alt_regions_id is not), it must be a list
    of sets of region names. Each set of region names will be run together as a
    single construction job, replacing the individual jobs for those regions.
    
    If merge_output_name is set, merges all constructed graphs together and
    outputs them under that name. Otherwise, outputs each graph constructed
    under its own name, but in a unified ID space.
    
    Returns a dict containing:
    
    'joined': a list of the unmerged, id-joined graph file IDs for each region.
    
    'merged': the merged graph file ID, if merge_output_name is set, or the
    only graph, if there is only one. None otherwise.
    
    'mapping': the file ID of the .mapping file produced by `vg ids --join`, if
    id joining had to happen. None otherwise.
    
    """

    # encapsulate follow-on
    child_job = Job()
    job.addChild(child_job)

    work_dir = job.fileStore.getLocalTempDir()

    if not regions:
        regions, region_names = [None], ['genome']
        
    if not alt_regions_id:
        # Coalesce regions (which we can't yet do if also running MSGA)
        regions, region_names = apply_coalesce(regions, region_names=region_names, coalesce_regions=coalesce_regions)
            
    region_graph_ids = []    
    for i, (region, region_name) in enumerate(zip(regions, region_names)):
        if not vcf_ids or (len(vcf_ids) > 1 and i >= len(vcf_ids)):
            # no vcf for region
            vcf_id = None
            tbi_id = None
            vcf_name = None
        elif len(vcf_ids) == 1:
            # special case: 1 vcf given, so assumed for all regions
            vcf_id = vcf_ids[0]
            tbi_id = tbi_ids[0]
            vcf_name = vcf_names[0]
        else:
            # one vcf per region
            vcf_id = vcf_ids[i]
            tbi_id = tbi_ids[i]
            vcf_name = vcf_names[i]
        fasta_id = fasta_ids[0] if len(fasta_ids) == 1 else fasta_ids[i]
        fasta_name = fasta_names[0] if len(fasta_names) == 1 else fasta_names[i]
        construct_region_job = child_job.addChildJobFn(run_construct_region_graph, context,
                                                       fasta_id, fasta_name,
                                                       vcf_id, vcf_name, tbi_id, region, region_name,
                                                       max_node_size, alt_paths, flat_alts, handle_svs,
                                                       # todo: bump as command line option?
                                                       #       also, needed if we update vg docker image?
                                                       is_chrom=not region or ':' not in region,
                                                       sort_ids=sort_ids,
                                                       normalize=normalize,
                                                       validate=validate,
                                                       cores=context.config.construct_cores,
                                                       memory=context.config.construct_mem,
                                                       disk=context.config.construct_disk)
        if alt_regions_id:
            region_graph_ids.append(construct_region_job.addFollowOnJobFn(run_msga, context, region_name + '.vg',
                                                                          construct_region_job.rv(),
                                                                          fasta_id,
                                                                          alt_regions_id,
                                                                          region,
                                                                          normalize=normalize,
                                                                          max_node_size=max_node_size,
                                                                          validate=validate,
                                                                          cores=context.config.alignment_cores,
                                                                          memory=context.config.alignment_mem,
                                                                          disk=context.config.alignment_disk).rv())
        else:
            region_graph_ids.append(construct_region_job.rv())

    return child_job.addFollowOnJobFn(run_join_graphs, context, region_graph_ids, join_ids,
                                      region_names, name, merge_output_name,
                                      cores=context.config.construct_cores,
                                      memory=context.config.xg_index_mem,
                                      disk=context.config.xg_index_disk).rv()

def run_join_graphs(job, context, region_graph_ids, join_ids, region_names, name, merge_output_name = None):
    """
    Join the ids of some graphs. If a merge_output_name is given, merge the
    graph files all together as well.
    
    Saves the unmerged, id-joined graphs, or the single merged graph if its
    name is given, to the output store. Also saves the node mapping file,
    produced from the `vg ids --join` call, to the output store.
    
    Skips doing any joining or merging if there is only one input graph.
    
    If join_ids is false, assumes the input graphs are already id-joined, and
    passes them through, merging if requested.
    
    Returns a dict containing:
    
    'joined': a list of the unmerged, id-joined graph file IDs (or the input
    graph(s) re-uploaded as output files if no joining occurred)
    
    'joined_names': a list of .vg filenames for those graphs 
    
    'merged': the merged graph file ID, if merging occurred, or the only input
    graph ID, if there was only one. None otherwise.
    
    'mapping': the file ID of the .mapping file produced by `vg ids --join`, if
    run. None otherwise.
    
    """
        
    work_dir = job.fileStore.getLocalTempDir()

    # Download graph for each region.
    # To keep command line lengths short we name the files by numbers.
    region_files = []
    # But track their human-readable names
    human_names = []
    for number, (region_graph_id, region_name) in enumerate(zip(region_graph_ids, region_names)):
        region_file = '{}.vg'.format(number)
        job.fileStore.readGlobalFile(region_graph_id, os.path.join(work_dir, region_file), mutable=True)
        region_files.append(region_file)
        human_names.append(remove_ext(region_name, '.vg') + '.vg')

    if merge_output_name:
        merge_output_name = remove_ext(merge_output_name, '.vg') + '.vg'
        
    # This is our return value. Initialize it as empty but with all the keys
    # set to make asking for things with .rv() easier.
    to_return = {
        'joined': [],
        'joined_names': human_names,
        'merged': None,
        'mapping': None
    }
    
    if join_ids and len(region_files) != 1:
        # The graphs aren't pre-joined, and we have more than one.
        # Do the actual joining
        
        mapping_file = merge_output_name[:-3] if merge_output_name else name
        mapping_file = os.path.join(work_dir, mapping_file + '.mapping')

        # join the ids
        cmd = ['vg', 'ids', '--join', '--mapping', os.path.basename(mapping_file)] + region_files
        context.runner.call(job, cmd, work_dir=work_dir)
        
        # save the mapping file
        to_return['mapping'] = context.write_intermediate_file(job, mapping_file)
    
    if merge_output_name is not None:
        # We want a single merged output file, so merge the graphs that we now know are in a joined ID space.
        
        # Make sure we aren't writing to an input file
        assert merge_output_name not in region_files
        
        # Run vg to combine into that file
        cmd = ['vg', 'combine'] + region_files
        with open(os.path.join(work_dir, merge_output_name), 'wb') as merge_file:
            context.runner.call(job, cmd, work_dir=work_dir, outfile = merge_file)
                    
        # And write the merged graph as an output file
        to_return['merged'] = context.write_output_file(job, os.path.join(work_dir, merge_output_name))
        
        if join_ids and len(region_files) != 1:
            # If we do all the merging, and we made new joined graphs, write the joined graphs as intermediates
            to_return['joined'] = [context.write_intermediate_file(job, os.path.join(work_dir, f), dest)
                                   for f, dest in zip(region_files, human_names)]
        else:
            # We can just pass through the existing intermediate files without re-uploading
            to_return['joined'] = region_graph_ids
    else:
        # No merging happened, so the id-joined files need to be output files.
        # We assume they came in as intermediate files, even if we didn't join them.
        # So we defintiely have to write them.
        # And we need to make sure to write them under their assigned
        # region-based names, even if we downloaded them to shorter names.
        to_return['joined'] = [context.write_output_file(job, os.path.join(work_dir, f), dest)
                               for f, dest in zip(region_files, human_names)]
                    
    return to_return 
        
    
def run_construct_region_graph(job, context, fasta_id, fasta_name, vcf_id, vcf_name, tbi_id,
                               region, region_name, max_node_size, alt_paths, flat_alts, handle_svs,
                               is_chrom = False, sort_ids = True, normalize = False, validate = False):
    """
    Construct a graph from the vcf for a given region and return its file id.
    
    region may be a FASTA contig id, a set of FASTA contig IDs, or None for
    using the whole FASTA.
    
    If is_chrom is set, pass along that fact to the constructor so it doesn't
    try to pass a region out of the chromosome name.
    
    If sort_ids is set (the default), do a sort pass after construction to make
    sure the IDs come in topological order.
    
    If normalize is set, try to normalize the graph and merge splits that
    produce identical sequences.
    
    If validate is true, subject the graph to a `vg validate` pass after
    construct. This is off by default because vg currently does internal
    validation during construction.
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file)
    if vcf_id:
        vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
        job.fileStore.readGlobalFile(vcf_id, vcf_file)
        job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    cmd = ['vg', 'construct', '--reference', os.path.basename(fasta_file)]
    if vcf_id:
        cmd += ['--vcf', os.path.basename(vcf_file)]
    if region:
        if isinstance(region, str):
            cmd += ['--region', region]
        else:
            for fasta_id in region:
                cmd += ['--region', fasta_id]
        if is_chrom:
            cmd += ['--region-is-chrom']
    if max_node_size:
        cmd += ['--node-max', max_node_size]
    if alt_paths:
        cmd += ['--alt-paths']
    if flat_alts:
        cmd += ['--flat-alts']
    if handle_svs:
        cmd += ['--handle-sv']
    if job.cores:
        cmd += ['--threads', job.cores]

    if normalize or sort_ids:
        cmd = [cmd]

    if normalize:
        cmd.append(['vg', 'mod', '--until-normal', str(context.config.normalize_iterations), '-'])
        # can be done in single mod command, but weary of being sensitive to order of operations
        cmd.append(['vg', 'mod', '--chop', str(max_node_size), '-'])

    if sort_ids:
        cmd.append(['vg', 'ids', '--sort', '-'])

    vg_path = os.path.join(work_dir, region_name)
    try:
        with open(vg_path, 'wb') as vg_file:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = vg_file)
    except:
        # Dump everything we need to replicate the construction
        logging.error("Construction failed. Dumping files.")

        context.write_output_file(job, fasta_file)
        if vcf_id:
            context.write_output_file(job, vcf_file)
            context.write_output_file(job, vcf_file + '.tbi')
            
        raise
        
    if validate:
        # Check the constructed and possibly modified graph for errors
        context.runner.call(job, ['vg', 'validate', os.path.basename(vg_path)], work_dir = work_dir)

    return context.write_intermediate_file(job, vg_path)

def run_filter_vcf_samples(job, context, vcf_id, vcf_name, tbi_id, samples, vcf_subdir = None):
    """ 
    
    Use vcflib to remove all variants specifc to a set of samples.
    Keep all the sample data in the VCF except that for the sample that was removed.
    
    """
    if not samples:
        # We can exclude nothing with a no-op
        return vcf_id, tbi_id
    
    work_dir = job.fileStore.getLocalTempDir()

    # Download the original VCF
    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
    # Where will the final filtered VCF go?
    filter_vcf_name = '{}_filter.vcf.gz'.format(vcf_base)
    # What intermediate VCF will we use for variants to drop?
    private_vcf_name = '{}_private.vcf.gz'.format(vcf_base)

    # Make a VCF with only the variants for the sample we want gone
    # TODO: if none of the samples listed are present, we get *all* variants instead of no variants.
    # Then we proceed to remove all the variants in the isec step.
    # Can we detect/avoid this?
    cmd = ['bcftools', 'view', os.path.basename(vcf_file), '--private',
           '--samples', ','.join(samples), '--force-samples', '--output-type', 'z']
    with open(os.path.join(work_dir, private_vcf_name), 'wb') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)
        
    # bcftools isec demands indexed input, so index the itnermediate file.
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', private_vcf_name],
                        work_dir=work_dir)
                        
    # Now make a VCF that excludes those variants and also excludes the filtered-out samples.
    # We subtract the private variants from the original VCF, and then remove the samples we're excluding.
    cmd = [['bcftools', 'isec', '--complement', os.path.basename(vcf_file), os.path.basename(private_vcf_name),
            '--write', '1'],
           ['bcftools', 'view', '-', '--samples', '^' + (','.join(samples)), '--trim-alt-alleles',
            '--force-samples', '--output-type', 'z']]
    with open(os.path.join(work_dir, filter_vcf_name), 'wb') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)

    if vcf_subdir:
        write_fn = lambda x: context.write_output_file(job, x, out_store_path = os.path.join(vcf_subdir, os.path.basename(x)))
    else:
        write_fn = lambda x: context.write_intermediate_file(job, x)
        
    # Upload the final VCF
    out_vcf_id = write_fn(os.path.join(work_dir, filter_vcf_name))

    # Index it
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', filter_vcf_name],
                        work_dir=work_dir)
                                        
    # And upload the index
    out_tbi_id = write_fn(os.path.join(work_dir, filter_vcf_name) + '.tbi')
    
    return out_vcf_id, out_tbi_id
    
def run_make_control_vcfs(job, context, vcf_id, vcf_name, tbi_id, sample, pos_only = False,
                          vcf_subdir = None, no_filter_if_sample_not_found = False):
    """ make a positive and negative control vcf 
    The positive control has only variants in the sample, the negative
    control has only variants not in the sample
    """

    assert sample is not None
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    # In some cases, our sample may be missing from a chromosome (ex NA12878 from Y in 1000 Genomes)
    # bcftools -s won't work so we handle here as a special case, assuming no sample means no variants
    cmd = ['bcftools', 'query', '--list-samples', os.path.basename(vcf_file)]
    found_samples = context.runner.call(job, cmd, work_dir=work_dir, check_output=True)
    found_sample = sample in found_samples.decode().strip().split('\n')

    # Hacky interface to not do anything if we can't find the sample.
    # By default, we'd return an empty VCF in this case
    if not found_sample and no_filter_if_sample_not_found:
        return vcf_id, tbi_id, vcf_id, tbi_id

    # filter down to sample in question
    cmd = [['bcftools', 'view', os.path.basename(vcf_file), '--samples', sample, '--trim-alt-alleles']]
    
    if found_sample:
        # remove anything that's not alt (probably cleaner way to do this)
        gfilter = 'GT="0" || GT="0|0" || GT="0/0"'
        gfilter += ' || GT="." || GT=".|." || GT="./."'
        gfilter += ' || GT=".|0" || GT="0/."'
        gfilter += ' || GT="0|." || GT="./0"'

        cmd.append(['bcftools', 'view', '-', '--output-type', 'z', '--exclude', gfilter])
    else:
        # if the sample isn't in the vcf, then there are no variants of interest, so
        # we report a header field without any samples
        cmd[0] += ['--force-samples', '--header-only', '--output-type', 'z']            

    out_pos_name = remove_ext(remove_ext(os.path.basename(vcf_name), '.gz'), '.vcf')
    out_neg_name = out_pos_name + '_minus_{}.vcf.gz'.format(sample)
    out_pos_name += '_{}.vcf.gz'.format(sample)

    with open(os.path.join(work_dir, out_pos_name), 'wb') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)

    context.runner.call(job, ['tabix', '--force', '--preset', 'vcf', out_pos_name], work_dir=work_dir)

    # we don't write vcfs to the output store unless we have a subdir to dump them in
    if vcf_subdir:
        def write_fn(local_path, out_store_path = None):
            os_name = os.path.basename(local_path) if not out_store_path else os.path.basename(out_store_path)
            return context.write_output_file(job, local_path, os.path.join(vcf_subdir, os_name))
    else:
        def write_fn(local_path, out_store_path = None):
            return context.write_intermediate_file(job, local_path)

    pos_control_vcf_id = write_fn(os.path.join(work_dir, out_pos_name))
    pos_control_tbi_id = write_fn(os.path.join(work_dir, out_pos_name + '.tbi'))

    if pos_only:
        return pos_control_vcf_id, pos_control_tbi_id, None, None

    # subtract the positive control to make the negative control
    cmd = ['bcftools', 'isec', os.path.basename(vcf_file), out_pos_name, '-p', 'isec', '-O', 'z']
    context.runner.call(job, cmd, work_dir=work_dir)

    context.runner.call(job, ['tabix', '--force', '--preset', 'vcf', 'isec/0000.vcf.gz'], work_dir=work_dir)

    neg_control_vcf_id = write_fn(os.path.join(work_dir, 'isec', '0000.vcf.gz'), out_store_path = out_neg_name)
    neg_control_tbi_id = write_fn(os.path.join(work_dir, 'isec', '0000.vcf.gz.tbi'), out_store_path = out_neg_name + '.tbi')

    return pos_control_vcf_id, pos_control_tbi_id, neg_control_vcf_id, neg_control_tbi_id

def run_min_allele_filter_vcf_samples(job, context, vcf_id, vcf_name, tbi_id, min_af, vcf_subdir = None):
    """
    filter a vcf by allele frequency using bcftools --min-af
    """
    if not min_af:
        return vcf_id, tbi_id
    
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
    af_vcf_name = '{}_minaf_{}.vcf.gz'.format(vcf_base, min_af)

    cmd = ['bcftools', 'view', '--min-af', min_af, '-O', 'z', os.path.basename(vcf_file)]
    with open(os.path.join(work_dir, af_vcf_name), 'wb') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile=out_file)

    if vcf_subdir:
        write_fn = lambda x: context.write_output_file(job, x, out_store_path = os.path.join(vcf_subdir, os.path.basename(x)))
    else:
        write_fn = lambda x: context.write_intermediate_file(job, x)

    out_vcf_id = write_fn(os.path.join(work_dir, af_vcf_name))

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', af_vcf_name],
                        work_dir=work_dir)
                                        
    out_tbi_id = write_fn(os.path.join(work_dir, af_vcf_name) + '.tbi')
    
    return out_vcf_id, out_tbi_id

def run_make_haplo_indexes(job, context, vcf_ids, tbi_ids, vcf_names, vg_ids, vg_names,
                           output_name, regions, sample, intermediate = False):
    """
    return xg/gbwt for each chromosome for extracting haplotype thread graphs
    (to simulate from) or sample graphs (as positive control)
    """
    assert(sample is not None)

    # make sure we're only dealing with chrom names (should probably be error otherwise)
    chroms = [region[0:region.find(':')] if ':' in region else region for region in regions]

    # Drop Nones from the VCF names; for some reason it is getting padded with Nones.
    # TODO: Work out where/why that is happening and stop it.
    vcf_names = [v for v in vcf_names if v is not None]

    logger.debug('Making gbwt for {} vgs, {} chroms, {} vcfs, {} tbis, and {} vcf names'.format(
        len(vg_ids), len(chroms), len(vcf_ids), len(tbi_ids), len(vcf_names)))

    # validate options should enforce this but check to be sure assumptions met to avoid
    # returning nonsense
    assert len(vg_ids) == len(regions)
    assert len(vcf_ids) == 1 or len(vcf_ids) <= len(regions)
    assert len(tbi_ids) == len(vcf_ids)
    assert len(vcf_names) == len(vcf_ids)
    
    logger.info('Making gbwt for chromosomes {}'.format(chroms))

    xg_ids = []
    gbwt_ids = []
    
    for i, (vg_id, vg_name, region) in enumerate(zip(vg_ids, vg_names, chroms)):
        if len(vcf_names) == 1: 
            # One VCF for all contigs
            vcf_name = vcf_names[0]
            vcf_id = vcf_ids[0]
            tbi_id = tbi_ids[0]
        elif i < len(vcf_names):
            # One VCF for this contig
            vcf_name = vcf_names[i]
            vcf_id = vcf_ids[i]
            tbi_id = tbi_ids[i]
        else:
            # No VCF for this contig
            vcf_name = None
            vcf_id = None
            tbi_id = None
            
        # index the graph and vcf to make the gbwt
        xg_name = remove_ext(vg_name, '.vg')
        xg_job = job.addChildJobFn(run_xg_indexing, context, [vg_id], [vg_name],
                                   xg_name, vcf_id, tbi_id,
                                   make_gbwt=True,
                                   intermediate=intermediate,
                                   cores=context.config.xg_index_cores,
                                   memory=context.config.xg_index_mem,
                                   disk=context.config.xg_index_disk)
        xg_ids.append(xg_job.rv(0))
        gbwt_ids.append(xg_job.rv(1))

    return xg_ids, gbwt_ids

def run_make_haplo_graphs(job, context, vg_ids, vg_names, xg_ids,
                          output_name, regions, sample, haplotypes, gbwt_ids,
                          intermediate = False):
    """
    Make some haplotype graphs for threads in a gbwt. regions must be defined
    since we use the chromosome name to get the threads. Also, gbwt_ids must be
    specified (one genome gbwt or one per region).
    
    Returns a list of haplotypes, where each haplotype is a list of vg graphs subset to that haplotype.
    """

    assert(sample is not None)

    # ith element will be a list of graphs (1 list / region) for haplotype i
    thread_vg_ids = []
    for h in haplotypes:
        thread_vg_ids.append([])

    # make sure we're only dealing with chrom names (should probably be error otherwise)
    chroms = [region[0:region.find(':')] if ':' in region else region for region in regions]

    # validate options should enforce this but check to be sure assumptions met to avoid
    # returning nonsense
    assert len(vg_ids) == len(regions)
    
    logger.info('Making haplo graphs for chromosomes {}'.format(chroms))
    
    for i, (vg_id, vg_name, region, xg_id) in enumerate(zip(vg_ids, vg_names, chroms, xg_ids)):
        # make a thread graph from the xg
        assert not gbwt_ids or len(gbwt_ids) in [1, len(xg_ids)]
        # support whole-genome or chromosome gbwts
        gbwt_id = None if not gbwt_ids else gbwt_ids[0] if len(gbwt_ids) == 1 else gbwt_ids[i]
        hap_job = job.addChildJobFn(run_make_haplo_thread_graphs, context, vg_id, vg_name,
                                    output_name, [region], xg_id, sample, haplotypes, gbwt_id,
                                    intermediate = intermediate,
                                    cores=context.config.construct_cores,
                                    memory=context.config.construct_mem,
                                    disk=context.config.construct_disk)
        for j in range(len(haplotypes)):
            thread_vg_ids[j].append(hap_job.rv(j))

    return thread_vg_ids

def run_make_haplo_thread_graphs(job, context, vg_id, vg_name, output_name, chroms, xg_id,
                                 sample, haplotypes, gbwt_id, intermediate = False):
    """
    make some haplotype graphs for threads in a gbwt
    
    If there are no haplotypes in the gbwt, passes through the portion of the input graph covered by paths.
    """
    work_dir = job.fileStore.getLocalTempDir()

    xg_path = os.path.join(work_dir, vg_name[:-3] + '.xg')
    job.fileStore.readGlobalFile(xg_id, xg_path)

    vg_path = os.path.join(work_dir, vg_name)
    job.fileStore.readGlobalFile(vg_id, vg_path)

    if gbwt_id:
        gbwt_path = os.path.join(work_dir, vg_name[:-3] + '.gbwt')
        job.fileStore.readGlobalFile(gbwt_id, gbwt_path)
    
        # Check if there are any threads in the index
        # TODO: Won't be useful if the index covers multiple contigs because we aren't indexing one contig graph at a time.
        try:
            thread_count = int(context.runner.call(job,
                [['vg', 'paths', '--list', '--gbwt', os.path.basename(gbwt_path)], 
                ['wc', '-l']], work_dir = work_dir, check_output = True))
        except:
            # TODO: vg paths really needs to be fixed to be able to check for 0 threads without failing
            RealtimeLogger.warning("No GBWT threads found in {}.  Using reference path for haplotype extraction".format(
                os.path.basename(gbwt_path)))
            thread_count = 0
            
    else:
        # No gbwt means no threads
        thread_count = 0
    

    thread_vg_ids = []

    for hap in haplotypes:
        
        # This can't work if the sample is None and we want any haplotypes
        assert(sample is not None)

        try:
            # Work out a tag for this graph, depending on whether it belongs to one chromosome or not
            tag = '_{}'.format(chroms[0]) if len(chroms) == 1 else ''
        
            if thread_count == 0:
                # We have no haplotype data on this contig. This is something
                # like chrM, and we want to pass through the ref version.
                vg_with_thread_as_path_path = vg_path
            else:
                # We know we have haplotype data on this contig.
                # Pull out the graph with just the haplotype thread as the only path to vg_with_thread_as_path_path
                
                # To accomplish this we now need to make sure to use vg combine
                # to combine the path-only vg Protobuf and the actual graph. So
                # first get them in different files.
                
                base_graph_filename = '{}{}_thread_{}_base.vg'.format(output_name, tag, hap)
                
                # strip paths from our original graph            
                cmd = ['vg', 'paths', '-d', '-v', os.path.basename(vg_path)]
                with open(os.path.join(work_dir, base_graph_filename), 'wb') as out_file:
                    context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)
                    
                path_graph_filename = '{}{}_thread_{}_path.vg'.format(output_name, tag, hap)

                # get haplotype thread paths from the gbwt
                cmd = ['vg', 'paths', '--gbwt', os.path.basename(gbwt_path), '--extract-vg', '-x', os.path.basename(xg_path)]
                for chrom in chroms:
                    cmd += ['-q', '_thread_{}_{}_{}'.format(sample, chrom, hap)]
                with open(os.path.join(work_dir, path_graph_filename), 'wb') as out_file:
                    context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)
                
                # Now combine the two files, adding the paths to the graph
                vg_with_thread_as_path_path = os.path.join(work_dir, '{}{}_thread_{}_merge.vg'.format(output_name, tag, hap))
                logger.info('Creating thread graph {}'.format(vg_with_thread_as_path_path))
                cmd = ['vg', 'combine', '-c', base_graph_filename, path_graph_filename]
                with open(vg_with_thread_as_path_path, 'wb') as out_file:
                    context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)
                    
                # Now delete the intermediates
                os.unlink(os.path.join(work_dir, base_graph_filename))
                os.unlink(os.path.join(work_dir, path_graph_filename))
                    
            # Now trim the graph vg_with_thread_as_path_path into vg_trimmed_path, dropping anything not covered by a path
            vg_trimmed_path = os.path.join(work_dir, '{}{}_thread_{}.vg'.format(output_name, tag, hap))
            logger.info('Creating trimmed thread graph {}'.format(vg_trimmed_path))
            with open(vg_trimmed_path, 'wb') as trimmed_file:
                # Then we trim out anything other than our thread path
                cmd = [['vg', 'mod', '-N', os.path.basename(vg_with_thread_as_path_path)]]
                # And get rid of our thread paths since they take up lots of space when re-indexing
                filter_cmd = ['vg', 'paths', '-v', '-']
                for chrom in chroms:
                    filter_cmd += ['--retain-paths', chrom]
                cmd.append(filter_cmd)
                context.runner.call(job, cmd, work_dir = work_dir, outfile = trimmed_file)

            write_fn = context.write_intermediate_file if intermediate else context.write_output_file
            thread_vg_ids.append(write_fn(job, vg_trimmed_path))
            
        except:
            # Dump everything we need to replicate the thread extraction
            logging.error("Thread extraction failed. Dumping files.")

            context.write_output_file(job, vg_path)
            context.write_output_file(job, xg_path)
            if gbwt_id:
                context.write_output_file(job, gbwt_path)
            
            raise

    logger.info('Got {} thread file IDs'.format(len(thread_vg_ids)))

    return thread_vg_ids

def run_make_sample_graphs(job, context, vg_ids, vg_names, xg_ids,
                           output_name, regions, sample, gbwt_ids):
    """
    Make some sample graphs for threads in a gbwt. regions must be defined
    since we use the chromosome name to get the threads. Also, gbwt_ids must be
    specified (one genome gbwt or one per region).
    """

    assert(sample is not None)

    # ith element will be a sample graph for region i
    sample_vg_ids = []

    # make sure we're only dealing with chrom names (should probably be error otherwise)
    chroms = [region[0:region.find(':')] if ':' in region else region for region in regions]

    # validate options should enforce this but check to be sure assumptions met to avoid
    # returning nonsense
    assert len(vg_ids) == len(regions)
    
    logger.info('Making sample graphs for chromosomes {}'.format(chroms))
    
    for i, (vg_id, vg_name, region, xg_id) in enumerate(zip(vg_ids, vg_names, chroms, xg_ids)):
        # make a thread graph from the xg
        assert not gbwt_ids or len(gbwt_ids) in [1, len(xg_ids)]
        # support whole-genome or chromosome gbwts
        gbwt_id = gbwt_ids[0] if len(gbwt_ids) == 1 else gbwt_ids[i]        
        hap_job = job.addChildJobFn(run_make_sample_region_graph, context, vg_id, vg_name,
                                    output_name, region, xg_id, sample, [0,1], gbwt_id,
                                    cores=context.config.construct_cores,
                                    memory=context.config.construct_mem,
                                    disk=context.config.construct_disk)
        sample_vg_ids.append(hap_job.rv())

    return sample_vg_ids

def run_make_sample_region_graph(job, context, vg_id, vg_name, output_name, chrom, xg_id,
                                 sample, haplotypes, gbwt_id, leave_thread_paths=True, validate=True):
    """
    make a sample graph using the gbwt.
    
    Extract the subgraph visited by threads for the requested sample, if it is nonempty.
    Otherwise (for cases like chrM where there are no variant calls and no threads) pass through
    the primary path of the graph.
    
    If validate is True (the default), makes sure the final graph passes
    `vg validate` before sending it on.
    """

    # This can't work if the sample is None and we want any haplotypes
    assert(sample is not None)

    work_dir = job.fileStore.getLocalTempDir()

    xg_path = os.path.join(work_dir, vg_name[:-3] + '.xg')
    job.fileStore.readGlobalFile(xg_id, xg_path)

    vg_path = os.path.join(work_dir, vg_name)
    job.fileStore.readGlobalFile(vg_id, vg_path)

    gbwt_path = os.path.join(work_dir, vg_name[:-3] + '.gbwt')
    if gbwt_id:
        job.fileStore.readGlobalFile(gbwt_id, gbwt_path)
        
    # Check if there are any threads in the index
    assert(gbwt_id)
    thread_count = int(context.runner.call(job,
        [['vg', 'paths', '--threads', '--list', '--gbwt', os.path.basename(gbwt_path), '-x',  os.path.basename(xg_path)], 
        ['wc', '-l']], work_dir = work_dir, check_output = True))
    if thread_count == 0:
        # There are no threads in our GBWT index (it is empty).
        # This means that we have no haplotype data for this graph.
        # This means the graph's contigs contig probably should be included,
        # in at least their reference versions, in all graphs.
        # Use the whole graph as our "extracted" graph, which we 
        # will then pare down to the part covered by paths (i.e. the primary path)
        extract_graph_path = vg_path
    else:
        # We have actual thread data for the graph. Go extract the relevant threads.
        extract_graph_path = os.path.join(work_dir, '{}_{}_extract.vg'.format(output_name, chrom))
        logger.info('Creating sample extraction graph {}'.format(extract_graph_path))
        with open(extract_graph_path, 'wb') as extract_graph_file:
            # strip paths from our original graph            
            cmd = ['vg', 'paths', '-d', '-v', os.path.basename(vg_path)]
            context.runner.call(job, cmd, work_dir = work_dir, outfile = extract_graph_file)

            for hap in haplotypes:
                # get haplotype thread paths from the index
                if gbwt_id:
                    cmd = ['vg', 'paths', '--gbwt', os.path.basename(gbwt_path), '--extract-vg']
                else:
                    cmd = ['vg', 'find']
                cmd += ['-x', os.path.basename(xg_path)]
                cmd += ['-q', '_thread_{}_{}_{}'.format(sample, chrom, hap)]
                context.runner.call(job, cmd, work_dir = work_dir, outfile = extract_graph_file)

    sample_graph_path = os.path.join(work_dir, '{}_{}.vg'.format(output_name, chrom))
    logger.info('Creating sample graph {}'.format(sample_graph_path))
    with open(sample_graph_path, 'wb') as sample_graph_file:
        # Then we trim out anything other than our thread paths
        cmd = [['vg', 'mod', '-N', os.path.basename(extract_graph_path)]]
        if not leave_thread_paths:
            cmd.append(['vg', 'paths', '-v', '-', '-d'])
        context.runner.call(job, cmd, work_dir = work_dir, outfile = sample_graph_file)
        
    if validate:
        # Make sure that the resulting graph passes validation before returning it.
        # This is another whole graph load and so will take a while.
        context.runner.call(job, ['vg', 'validate', os.path.basename(sample_graph_path)], work_dir = work_dir)

    sample_vg_id = context.write_intermediate_file(job, sample_graph_path)
            
    return sample_vg_id

def construct_main(context, options):
    """
    Wrapper for vg constructing. 
    """

    # check some options
    validate_construct_options(options)

    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    # Merge up all filter samples
    filter_samples = []
    if options.filter_samples:
        filter_samples += options.filter_samples
    if options.filter_ceph:
        filter_samples += CEPH_SAMPLES
    filter_samples = list(set(filter_samples))

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            inputFastaFileIDs = [importer.load(fasta) for fasta in options.fasta]
            inputFastaNames = [os.path.basename(fasta) for fasta in options.fasta]

            inputVCFFileIDs = []
            inputVCFNames = []
            inputTBIFileIDs = []
            for vcf_batch in options.vcf:
                inputVCFFileIDs.append([importer.load(make_url(vcf)) for vcf in vcf_batch.split(',')])
                inputVCFNames.append([os.path.basename(vcf) for vcf in vcf_batch.split(',')])
                inputTBIFileIDs.append([importer.load(make_url(vcf + '.tbi'), wait_on = inputVCFFileIDs[-1][i]) \
                                        for i, vcf in enumerate(vcf_batch.split(','))])
            
            inputBWAFastaID=None
            if options.bwa_reference:
                inputBWAFastaID = importer.load(options.bwa_reference)

            inputRegionsFileID = None
            if options.regions_file:
                inputRegionsFileID = importer.load(options.regions_file)

            alt_regions_id = importer.load(options.alt_regions_bed) if options.alt_regions_bed else None
            coalesce_regions_id = importer.load(options.coalesce_regions) if options.coalesce_regions else None 

            importer.wait()
            inputFastaFileIDs = importer.resolve(inputFastaFileIDs)
            inputVCFFileIDs = importer.resolve(inputVCFFileIDs)
            inputTBIFileIDs = importer.resolve(inputTBIFileIDs)
            inputBWAFastaID = importer.resolve(inputBWAFastaID)
            inputRegionsFileID = importer.resolve(inputRegionsFileID)
            alt_regions_id = importer.resolve(alt_regions_id)
            coalesce_regions_id = importer.resolve(coalesce_regions_id)

            # We only support one haplotype extraction sample (enforced by validate) despire what CLI implies
            haplo_extraction_sample = options.haplo_sample if options.haplo_sample else options.sample_graph       
                   
            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
            # Current job in follow-on chain
            cur_job = init_job

            # Unzip the fasta
            for i, fasta in enumerate(options.fasta):
                if fasta.endswith('.gz'):
                    inputFastaFileIDs[i] = init_job.addChildJobFn(run_unzip_fasta, context, inputFastaFileIDs[i], 
                                                                  os.path.basename(fasta),
                                                                  disk=context.config.construct_disk).rv()
                    inputFastaNames[i] = inputFastaNames[i][:-3]

            # Mask out ambigous bases
            if options.mask_ambiguous:
                mask_root = Job()
                cur_job.addFollowOn(mask_root)
                cur_job = mask_root
                for i, (fasta_id, fasta_name) in enumerate(zip(inputFastaFileIDs, inputFastaNames)):
                    inputFastaFileIDs[i] = mask_root.addChildJobFn(run_mask_ambiguous, context, inputFastaFileIDs[i], inputFastaNames[i],
                                                                   disk=context.config.construct_disk).rv(0)

            # do minimum allele frequency filter as preprocessing step
            if options.pre_min_af:
                min_af_job = Job()
                cur_job.addFollowOn(min_af_job)
                cur_job = min_af_job
                af_vcf_ids_list, af_tbi_ids_list = [], []
                for vcf_ids, vcf_names, tbi_ids in zip(inputVCFFileIDs, inputVCFNames, inputTBIFileIDs):
                    af_vcf_ids, af_tbi_ids = [], []
                    for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
                        af_job = min_af_job.addChildJobFn(run_min_allele_filter_vcf_samples, context, vcf_id,
                                                          vcf_name, tbi_id, options.pre_min_af,
                                                          cores=context.config.construct_cores,
                                                          memory=context.config.construct_mem,
                                                          disk=context.config.construct_disk)
                        af_vcf_ids.append(af_job.rv(0))
                        af_tbi_ids.append(af_job.rv(1))
                    af_vcf_ids_list.append(af_vcf_ids)
                    af_tbi_ids_list.append(af_tbi_ids)
                inputVCFFileIDs, inputTBIFileIDs = af_vcf_ids_list, af_tbi_ids_list
                    
            regions_regex = None if not options.regions_regex else '|'.join(options.regions_regex)

            # Parse the regions from file
            if options.regions_file:
                cur_job = cur_job.addFollowOnJobFn(run_scan_regions_file, context, inputRegionsFileID, regions_regex,
                                                   memory=context.config.misc_mem,
                                                   disk=context.config.misc_disk)
                regions = cur_job.rv()
            elif options.fasta_regions:
                # Extract fasta sequence names and append them to regions
                # Make sure we have a plausible amount of disk for downloading it
                cur_job = cur_job.addFollowOnJobFn(run_scan_fasta_sequence_names, context,
                                                   inputFastaFileIDs[0],
                                                   inputFastaNames[0],
                                                   options.regions,
                                                   regions_regex,
                                                   memory=context.config.misc_mem,
                                                   disk=context.config.construct_disk)
                regions = cur_job.rv()
            else:
                regions = options.regions          

            # Preproces chromosome names everywhere to be consistent,
            # either mapping from 1-->chr1 etc, or going the other way.
            # Deduplicate regions that become the same after renaming.
            if options.add_chr_prefix or options.remove_chr_prefix:
                cur_job = cur_job.addFollowOnJobFn(run_fix_chrom_names, context,
                                                   options.add_chr_prefix,
                                                   regions,
                                                   inputFastaFileIDs,
                                                   inputFastaNames,
                                                   inputVCFFileIDs,
                                                   inputVCFNames,
                                                   inputTBIFileIDs,
                                                   alt_regions_id,
                                                   cores=context.config.construct_cores,
                                                   memory=context.config.construct_mem,
                                                   disk=context.config.construct_disk)
                regions = cur_job.rv(0)
                inputFastaFileIDs, inputFastaFileNames = cur_job.rv(1), cur_job.rv(2)
                inputVCFFileIDs, inputTBIFileIDs = cur_job.rv(3), cur_job.rv(5)

            # Make sure that we don't have any alt sequences in our regions.  alt sequences
            # are inferred from the --target_regions bed file
            if alt_regions_id:
                cur_job = cur_job.addFollowOnJobFn(run_subtract_alt_regions,
                                                   context,
                                                   alt_regions_id,
                                                   regions)

                regions, alt_regions = cur_job.rv(0), cur_job.rv(1)
            else:
                alt_regions=[]
                
            if coalesce_regions_id:
                cur_job = cur_job.addFollowOnJobFn(run_read_coalesce_list,
                                                   context,
                                                   coalesce_regions_id)
                coalesce_regions = cur_job.rv()
            
            else:
                coalesce_regions=[]

            # Merge up comma-separated vcfs with bcftools merge
            cur_job = cur_job.addFollowOnJobFn(run_merge_all_vcfs, context,
                                               inputVCFFileIDs, inputVCFNames, inputTBIFileIDs)
            inputVCFFileIDs = cur_job.rv(0)
            inputVCFNames = cur_job.rv(1)
            inputTBIFileIDs = cur_job.rv(2)
                
            # Automatically make and name a bunch of vcfs
            vcf_job = cur_job.addFollowOnJobFn(run_generate_input_vcfs, context,
                                               inputVCFFileIDs, inputVCFNames, inputTBIFileIDs,
                                               regions,
                                               options.out_name,
                                               do_primary = options.primary,
                                               do_pan = options.pangenome,
                                               pos_control_sample = options.pos_control,
                                               neg_control_sample = options.neg_control,
                                               sample_graph = options.sample_graph,
                                               haplo_sample = options.haplo_sample,
                                               filter_samples = filter_samples,
                                               min_afs = options.min_af,
                                               vcf_subdir = '{}-vcfs'.format(options.out_name) if options.keep_vcfs else None)
            
            # Construct graphs
            vcf_job.addFollowOnJobFn(run_construct_all, context, inputFastaFileIDs,
                                     inputFastaNames, vcf_job.rv(),
                                     options.max_node_size, options.alt_paths or 'alt-gam' in options.indexes,
                                     options.flat_alts, options.handle_svs, regions,
                                     merge_graphs = options.merge_graphs,
                                     sort_ids = True, join_ids = True,
                                     wanted_indexes = options.indexes, 
                                     haplo_extraction_sample = haplo_extraction_sample,
                                     gbwt_prune = options.gbwt_prune,
                                     normalize = options.normalize,
                                     validate = options.validate,
                                     alt_regions_id = alt_regions_id,
                                     alt_regions = alt_regions,
                                     coalesce_regions = coalesce_regions)
                                     
            
            if inputBWAFastaID:
                # If we need to make a BWA index too, do that in parallel with everything else
                init_job.addFollowOnJobFn(run_bwa_index, context, inputBWAFastaID,
                                          copy_fasta=True,
                                          cores=context.config.bwa_index_cores, memory=context.config.bwa_index_mem,
                                          disk=context.config.bwa_index_disk)
                                     
            
            # Run the workflow
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
