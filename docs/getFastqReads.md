
# Pull sequencing reads into fastq.gz from local storage (Tony) or SRA


```bash
#Set directories for output
#Define data folders
out_dir="/input_dir/"
```

### Find the local fastq.gz file using 
[Tony frontend](http://iona.wi.mit.edu/tonydev/dbSearch.pl)


```bash
#Subsample local file while stripping STUPID fucking header
read_num=1000000
in_fq_1 = /lab/solexa_public/Youn/solexa_public/170505_WIGTC-HISEQ2A_CALADANXX/QualityScores
zcat "input_fq.gz" | sed '1 s/.*@/@/' | head -n $read_num | \
    gzip > $out_dir/input_down_fq.gz
```

    gzip: input_fq.gz: No such file or directory


## Example below pulls 1M reads for every type of seq data supported <br> Types of seq pulled:
    * Chip-Seq and Chromatin Accessibility (ATAC-Seq and DNAse-Seq)
    * ChiA-PET
    * Hi-C
    * Hi-ChIP
    * DNAse-HiC
    * RNA-Seq
    * Gro-Seq
    * 4C
### links to pipelines for each type of data below

# Don't forget!!!
to quality-check, trim, and filter your reads using [this pipeline](fq2preppedReads.md) before running ANY of the downstream pipelines


```bash
#Chip-seq/Chromatin Accessibility
#CTCF mesc
#/root/Euplotid/src/SRA2fq SRR524848 $out_dir CTCF_mesc 1000000
#Input mesc
/root/Euplotid/src/SRA2fq SRR524849 $out_dir input_mesc 5000000
#H3k27Ac mesc
/root/Euplotid/src/SRA2fq SRR066766 $out_dir h3k27ac_mesc 5000000
#ATAC mesc
#/root/Euplotid/src/SRA2fq SRR2927023 $out_dir ATAC_mesc 5000000

#CTCF hesc
#SRA2fq SRR2056018 $out_dir CTCF_hesc 1000000
#Input hesc
#/root/Euplotid/src/SRA2fq SRR2056020 $out_dir input_hesc 1000000
#H3K27Ac hesc
#/root/Euplotid/src/SRA2fq SRR2056016 $out_dir h3k27Ac_hesc 1000000
#ATAC h7 hesc
#/root/Euplotid/src/SRA2fq SRR3689760 $out_dir ATAC_h7 1000000
```

See [this pipeline](fq2peaks.md) to call peaks using Homer and/or MACS2 as well as nucleosome positioning using nucleoatac


```bash
#chia-PET
#mesc
/root/Euplotid/src/SRA2fq SRR1296617 $out_dir ChiA_SMC1_mesc 1000000
#hesc
/root/Euplotid/src/SRA2fq SRR2054933 $out_dir ChiA_SMC1_hesc 1000000
```

    

See [this pipeline](fq2ChIAInts.md) to prep reads, align, call and normalize pairwise interactions using ChiAPet2 and/or Origami and dump into cooler format.


```bash
#Hi-c
#mesc
/root/Euplotid/src/SRA2fq SRR443883 $out_dir HiC_mesc 1000000
#hesc
/root/Euplotid/src/SRA2fq SRR400260 $out_dir HiC_hesc 1000000
```

See [this pipeline](fq2HiCInts.md) to go from fastq reads, align, normalize and dump into cooler format using HiCPro


```bash
#Hi-ChIp
#mesc
/root/Euplotid/src/SRA2fq SRR3467183 $out_dir HiChip_mesc 1000000
#GM12878
/root/Euplotid/src/SRA2fq SRR3467176 $out_dir HiChip_hesc 1000000
```

See [this pipeline](fq2HiChIPInts.md) custom pipeline to go from fastq reads through HiCPro + scripts to normalize and dump into cooler format


```bash
#Dnase Hi-c
#mesc patski cells
/root/Euplotid/src/SRA2fq SRR2033066 $out_dir dnaseHiC_patski 1000000
#hesc
/root/Euplotid/src/SRA2fq SRR1248175 $out_dir dnaseHiC_hesc 1000000
```

See [this pipeline](fq2DNAseHiCInts.md) to go from fastq reads, align, normalize and dump into cooler format using HiCPro


```bash
#RNA-Seq
#mesc 4cell
/root/Euplotid/src/SRA2fq SRR1840518 $out_dir rnaseq_mesc 1000000
#hesc mesoderm
/root/Euplotid/src/SRA2fq SRR3439456 $out_dir rnaseq_hesc 1000000
```

See [this pipeline](fq2countsFPKM.md) to take RNA-Seq reads and align and quantify/normalize expression values (FPKM) using RSEM


```bash
#Gro-seq
#mesc
/root/Euplotid/src/SRA2fq SRR935093 $out_dir groseq_mesc 1000000
#h1 hesc (our data!)
/root/Euplotid/src/SRA2fq SRR574826 $out_dir groseq_hesc 1000000
```

See [this pipeline](fq2GroRPKM.md) find nascent transcripts using FStitch and miRNA promoters using mirSTP


```bash
#4C
#mesc poised enhancers = viewpoints
/root/Euplotid/src/SRA2fq SRR4451724 $out_dir 4c_poiEnh_mesc 1000000
#hesc MT2A
/root/Euplotid/src/SRA2fq SRR1409666 $out_dir 4c_MT2A_hesc 1000000
```

See [this pipeline](fq24CInts.md) to get wiggle file from fastq reads using HiCPro and/or custom pipeline
