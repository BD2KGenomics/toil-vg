#!/bin/bash
#################################################################################################
##
##  Script to setup the toil-vg python environment and download workflow inputs.
##
##  Inputs:
##
##  Assumptions:
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up a toil-vg python virtual environment for easier use of vg wdl and downloads
workflow inputs on the NIH Biowulf Cluster.

Inputs:
    -g PATH to the workflow input directory
    -v PATH to the toil-vg repository
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 2 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false

## Parse through arguments
while getopts "g:v:t:h" OPTION; do
    case $OPTION in
        g)
            WORKFLOW_INPUT_DIR=$OPTARG
        ;;
        v)
            TOIL_VG_DIR=$OPTARG
        ;;
        t)
            RUN_SMALL_TEST=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

module load git python/3.7

## Setup vg wdl python virtualenvironment
if [ ! -d "${TOIL_VG_DIR}" ]; then
    mkdir -p ${TOIL_VG_DIR}
    chmod 2770 ${TOIL_VG_DIR}
fi

cd ${TOIL_VG_DIR}
git clone --single-branch --branch vg_pedigree_workflow_dev https://github.com/vgteam/toil-vg.git 
git clone https://github.com/cmarkello/toil.git
python3 -m venv toilvg_venv
source toilvg_venv/bin/activate
pip install ./toil
pip install ./toil-vg
deactivate

## Setup and download workflow inputs
if [ ! -d "${WORKFLOW_INPUT_DIR}" ]; then
    mkdir -p ${WORKFLOW_INPUT_DIR}
    chmod 2770 ${WORKFLOW_INPUT_DIR}
fi

if [ $RUN_SMALL_TEST == false ]; then
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_whole_genome.txt -O ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_v1.27.0_graph_references/baseline_v1.27.0_90_ga64b70c1f_minaf_0.01.xg -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_v1.27.0_graph_references/baseline_v1.27.0_90_ga64b70c1f_minaf_0.01.gcsa -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_v1.27.0_graph_references/baseline_v1.27.0_90_ga64b70c1f_minaf_0.01.gcsa.lcp -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa.lcp
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_v1.27.0_graph_references/baseline_v1.27.0_90_ga64b70c1f_minaf_0.01.gbwt -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gbwt
else
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/path_list_21.txt -O ${WORKFLOW_INPUT_DIR}/path_list_21.txt
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.xg -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.xg
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gcsa -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gcsa.lcp -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa.lcp
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/vg_t289_graph_references/snp1kg_maf0.01_chr21_t289.gbwt -O ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gbwt
    wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/HG002.ped -O ${WORKFLOW_INPUT_DIR}/HG002.ped
fi
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.fai -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.dict -O ${WORKFLOW_INPUT_DIR}/hs37d5.dict
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/hs37d5.fa.gz -O ${WORKFLOW_INPUT_DIR}/hs37d5.fa.gz
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/snpEff_v4_3_GRCh37.75.zip -O ${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/genetic_map_GRCh37.tar -O ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar
wget https://storage.googleapis.com/cmarkell-vg-wdl-dev/eagle_data.tar -O ${WORKFLOW_INPUT_DIR}/eagle_data.tar

exit

