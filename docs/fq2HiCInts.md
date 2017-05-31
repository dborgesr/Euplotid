
### Install instructions required after first image pull
        cd /root/HiC-Pro/
        source activate py27
        R
        install.packages("ggplot2")
        install.packages("RColorBrewer")
        make configure
        make install



```python
#Take fastq in and spit out chilled normalized Hi-C matrix, use python 2.7 kernel
#Define data folders, data, and sample name
annotation="/input_dir/mm9"
tmp="/input_dir/"
input_dir="/input_dir/"
output_dir="/output_dir/"
input_fq_1="HiC_mesc_1_1M.fq.gz"
input_fq_2="HiC_mesc_2_1M.fq.gz"
sample_name="test"
bin_size="10000"
```


```python
!mkdir  $output_dir"/rawdata/mesc_test"
!cp $output_dir"$input_fq_1" $output_dir"/rawdata/mesc_test"test_1_1M.fq.gz
!cp $output_dir"$input_fq_2" $output_dir"/rawdata/mesc_test"test_2_1M.fq.gz
```


```python
%%writefile /root/HiC-Pro_2.8.1_devel/config-hicpro_mesc.txt
# %load /root/HiC-Pro_2.8.1_devel/config-hicpro.txt
# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = $tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 2
LOGFILE = hicpro.log

JOB_NAME = 
JOB_MEM = 
JOB_WALLTIME = 
JOB_QUEUE = 
JOB_MAIL = 

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _1
PAIR2_EXT = _2

#######################################################################
## Alignment options
#######################################################################

FORMAT = phred33
MIN_MAPQ = 0

BOWTIE2_IDX_PATH = 
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = mm9
GENOME_SIZE = chrom_hg19.sizes
CAPTURE_TARGET =

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = HindIII_resfrag_mm9.bed
LIGATION_SITE = AAGCTAGCTT
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1

```

    Overwriting /root/HiC-Pro_2.8.1_devel/config-hicpro_mesc.txt



```python
!/root/HiC-Pro_2.8.1_devel/bin/HiC-Pro \
    -i $output_dir -c /root/HiC-Pro_2.8.1_devel/config-hicpro_mesc.txt \
    -s mapping -s proc_hic -s quality_checks -s merge_persample -s build_contact_maps -s ice_norm \
    -o $output_dir
```

    make: *** No targets.  Stop.



```python
!/root/HiC-Pro_2.8.1_devel/bin/HiC-Pro -i $output_dir -c /root/HiC-Pro_2.8.1_devel/config-hicpro_mesc.txt -s mapping -o $output_dir
```

    make: *** No targets.  Stop.



```python
!cooler csort hg19.chrom.sizes [input_origami_contact_list]
!cooler cload tabix hg19.chrom.sizes:$bin_size  [input_origami_contact_list].sorted.txt.gz  [input_origami_contact_list].cool"
!cooler coarsegrain --no-balance  [input_origami_contact_list].cool"
```

    Usage: cooler csort [OPTIONS] PAIRS_PATH CHROMOSOMES_PATH
    
    Error: Invalid value for "pairs_path": Path "hg19.chrom.sizes" does not exist.
    /bin/sh: 1: Syntax error: Unterminated quoted string
    /bin/sh: 1: Syntax error: Unterminated quoted string



```python
!docker exec higlass-container python higlass-server/manage.py ingest_tileset \
    --filename "/data/[input_origami_contact_list].multires.cool" \
    --datatype matrix --filetype cooler --uid cooler-demo
```
