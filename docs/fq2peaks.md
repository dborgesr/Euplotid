
# 5. Fq2peaks
## Megatid required
Call ANY type of peak from reads stored in fq.gz files.
Can call **peaks** for:
* **factor**:
Peak finding for single contact or focal ChIP-Seq experiments or DNase-Seq.  This type of analysis is useful for transcription factors, and aims to identify the precise location of DNA-protein contact.  This type of peak finding uses a FIXED width peak size, which is automatically estimated from the Tag Autocorrelation.

* **histone**: 
Peak finding for broad regions of enrichment found in ChIP-Seq experiments for various histone marks.  This analysis finds variable-width peaks.

* **super**: 
Find Super Enhancers in your data

* **groseq**:
De novo transcript identification from strand specific GRO-Seq.  This attempts to identify transcripts from nascent RNA sequencing reads.

* **tss**
Identification of promoter/TSS from 5'RNA-Seq/CAGE or 5'GRO-Seq data.

* **dnase**
Adjusted parameters for DNase-Seq peak finding.

* **mC** 
DNA methylation analysis

## <a href="http://homer.ucsd.edu/homer/ngs/peaks.html"> Homer peak calling tutorial </a> for more info


```bash
#Define data folders, data, and sample name
#remember that fasta should have built bowtie2 indices
annotation="/input_dir/mm9"
tmp="/input_dir/"
input_dir="/input_dir/"
output_dir="/output_dir/"
input_fq_1="h3k27ac_mesc_5M.fq.gz"
input_fq_2=""
sample_name="mesc_h3k27ac_5M"
```

# Align reads to reference genome


```bash
#align and sort / filter bam using Bowtie2 and "sensitive" settings
#paired end
if [ -z "$input_dir""$input_fq_2" ]; then
    echo "Aligning using bowtie2 w/ paired end"
    bowtie2 -p 8 --sensitive -x  $annotation \
        -1 $input_dir$input_fq_1 -2 $input_dir$input_fq_2 | samtools view -bS - > $tmp$sample_name".bam"
else
    echo "Aligning using bowtie2 w/ single end"
    bowtie2 -p 8 --sensitive -x  $annotation \
        -U $input_dir$input_fq_1 | samtools view -bS - > $tmp$sample_name".bam"
fi

#sort index and filter for canonical chromosomes
samtools sort -@ 10 $tmp$sample_name".bam" -o $tmp$sample_name"_sorted.bam"
samtools index "$tmp$sample_name""_sorted.bam"
samtools view $tmp$sample_name"_sorted.bam" -hu chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX > "$output_dir$sample_name""_filt_sorted.bam"
samtools index "$output_dir$sample_name""_filt_sorted.bam"
```

    Aligning using bowtie2 w/ single end
    5000000 reads; of these:
      5000000 (100.00%) were unpaired; of these:
        704720 (14.09%) aligned 0 times
        3416197 (68.32%) aligned exactly 1 time
        879083 (17.58%) aligned >1 times
    85.91% overall alignment rate



```bash
#Take a peek at the bam
samtools view "$output_dir$sample_name""_filt_sorted.bam"| head
#get alignment stats
samtools flagstat "$output_dir$sample_name""_filt_sorted.bam"
samtools idxstats "$output_dir$sample_name""_filt_sorted.bam"
```

# Using Homer to call peaks
Homer is a pretty sweet collection of mainly perl scripts which can handle next-gen sequencing data to do things such as call peaks, below we use it to call different types of peaks, more information and a lot of knowledge is <a href="http://homer.ucsd.edu/homer/ngs/peaks.html"> Homer peak calling tutorial </a>  


```bash
#Homer run
#make tag directory file for Homer using aligned .bam file
makeTagDirectory $tmp/"$sample_name"_homer_tags "$output_dir$sample_name""_filt_sorted.bam"
```


```bash
#find Chromatin accessibility peaks 
findPeaks $tmp/"$sample_name"_homer_tags -style dnase -o auto
pos2bed.pl "$tmp""$sample_name"_homer_tags/peaks.txt > "$output_dir$sample_name""_open_peaks.bed"
```

    	Will parse file: /output_dir/mesc_atac_5M_filt_sorted.bam
    
    	Creating directory: /input_dir//mesc_atac_5M_homer_tags and removing existing *.tags.tsv
    
    	Treating /output_dir/mesc_atac_5M_filt_sorted.bam as a bam file
    	Reading alignment file /output_dir/mesc_atac_5M_filt_sorted.bam
    
    	Optimizing tag files...
    	Estimated genome size = 2638916841
    	Estimated average read density = 0.000322 per bp
    	Total Tags = 849746.0
    	Total Positions = 1534940
    	Average tag length = 124.5
    	Median tags per position = 0 (ideal: 1)
    	Average tags per position = 0.089
    	Fragment Length Estimate: 174
    	Peak Width Estimate: 607
    	Autocorrelation quality control metrics:
    		Same strand fold enrichment: 1.6
    		Diff strand fold enrichment: 4.4
    		Same / Diff fold enrichment: 0.3
    
    		Guessing sample is ChIP-Seq - uneven enrichment between same strand and
    		different strands - may have problems such as clonal amplification.
    
    	Fragment Length = 1
    	Total Tags = 849746.0
    	Tags per bp = 0.000425
    	Max tags per bp set automatically to 1.0
    	Finding peaks of size 150, no closer than 300
    		Finding peaks on chr1 (minCount=-0.9), total tags positions = 122489
    		Finding peaks on chr2 (minCount=-0.9), total tags positions = 112448
    		Finding peaks on chr3 (minCount=-0.9), total tags positions = 91766
    		Finding peaks on chr4 (minCount=-0.9), total tags positions = 93150
    		Finding peaks on chr5 (minCount=-0.9), total tags positions = 88987
    		Finding peaks on chr6 (minCount=-0.9), total tags positions = 91536
    		Finding peaks on chr7 (minCount=-0.9), total tags positions = 82812
    		Finding peaks on chr8 (minCount=-0.9), total tags positions = 74461
    		Finding peaks on chr9 (minCount=-0.9), total tags positions = 74872
    		Finding peaks on chr10 (minCount=-0.9), total tags positions = 81494
    		Finding peaks on chr11 (minCount=-0.9), total tags positions = 79499
    		Finding peaks on chr12 (minCount=-0.9), total tags positions = 67528
    		Finding peaks on chr13 (minCount=-0.9), total tags positions = 72764
    		Finding peaks on chr14 (minCount=-0.9), total tags positions = 65711
    		Finding peaks on chr15 (minCount=-0.9), total tags positions = 61693
    		Finding peaks on chr16 (minCount=-0.9), total tags positions = 60499
    		Finding peaks on chr17 (minCount=-0.9), total tags positions = 57690
    		Finding peaks on chr18 (minCount=-0.9), total tags positions = 55064
    		Finding peaks on chr19 (minCount=-0.9), total tags positions = 39903
    		Finding peaks on chrX (minCount=-0.9), total tags positions = 60574
    		Tags Used for cluster (less clonal tags) = 827104.5 / 849746.0
    	Expected tags per peak = 0.062033 (tbp = 0.000414)
    		Threshold	Peak Count	Expected Peak Count	FDR	Poisson
    		7	33.000	0.000	5.37e-07	6.64e-13
    		6	38.000	0.002	5.27e-05	7.50e-11
    		5	61.000	0.194	3.18e-03	7.27e-09
    		4	129.000	15.657	1.21e-01	5.87e-07
    		3	566.000	1012.764	1.00e+00	3.80e-05
    		2	10890.000	49234.353	1.00e+00	1.85e-03
    		1	197409.000	1603946.027	1.00e+00	6.01e-02
    		0	820539.000	26666666.667	1.00e+00	1.00e+00
    	0.10% FDR Threshold set at 6.0 (poisson pvalue ~ 7.50e-11)
    	38 peaks passed threshold
    	Local Background Filtering: 24 of 38 (63.16% passed)
    	Clonal filtering: 21 of 24 (87.50% passed)
    	Total Peaks identified = 21
    
    	Converted 21 peaks total
    



```bash
#find enhancers and super enhancers
findPeaks $tmp/"$sample_name"_homer_tags -style super -typical  $tmp/"$sample_name"_homer_tags/enh_peaks.txt -o auto
#convert regular enhancer peak file
pos2bed.pl $tmp/"$sample_name"_homer_tags/enh_peaks.txt > "$output_dir$sample_name""_enh_peaks.bed"
#save peaks as bed file in output directory
pos2bed.pl "$tmp""$sample_name"_homer_tags/peaks.txt > "$output_dir$sample_name""_se_peaks.bed"
```

## Use MACS2 to call peaks and NucleoATAC to call nucleosomes 


```bash
#Switch to python2.7 for MACS2 and NucleoATAC
source activate py27
#Call acessibility peaks, and using those peaks define area around to resolve nucleosomes
macs2 callpeak -t "$output_dir$sample_name""_filt_sorted.bam" --nomodel --shift -37 -f BAMPE \
    --extsize 73 --broad --keep-dup all -n $tmp"$sample_name""_MACS"
bedops --range 500:500 --everything $tmp"$sample_name""_MACS_peaks.broadPeak" | bedtools merge -i - > $tmp"$sample_name""_MACS_merged.bed"
nucleoatac run --bed $tmp"$sample_name""_MACS_merged.bed" --bam $output_dir"$sample_name""_filt_sorted.bam" --fasta "$annotation"".fa" --out $tmp"$sample_name""_nucATAC" --cores 8
source deactivate py27
```

    (py27) (py27) INFO  @ Tue, 30 May 2017 23:22:21: 
    # Command line: callpeak -t /output_dir/mesc_atac_5M_filt_sorted.bam --nomodel --shift -37 -f BAMPE --extsize 73 --broad --keep-dup all -n /input_dir/mesc_atac_5M_MACS
    # ARGUMENTS LIST:
    # name = /input_dir/mesc_atac_5M_MACS
    # format = BAMPE
    # ChIP-seq file = ['/output_dir/mesc_atac_5M_filt_sorted.bam']
    # control file = None
    # effective genome size = 2.70e+09
    # band width = 300
    # model fold = [5, 50]
    # qvalue cutoff for narrow/strong regions = 5.00e-02
    # qvalue cutoff for broad/weak regions = 1.00e-01
    # Larger dataset will be scaled towards smaller dataset.
    # Range for calculating regional lambda is: 10000 bps
    # Broad region calling is on
    # Paired-End mode is on
     
    INFO  @ Tue, 30 May 2017 23:22:21: #1 read fragment files... 
    INFO  @ Tue, 30 May 2017 23:22:21: #1 read treatment fragments... 
    INFO  @ Tue, 30 May 2017 23:22:35:  1000000 
    INFO  @ Tue, 30 May 2017 23:22:39: #1 mean fragment size is determined as 229 bp from treatment 
    INFO  @ Tue, 30 May 2017 23:22:39: #1 fragment size = 229 
    INFO  @ Tue, 30 May 2017 23:22:39: #1  total fragments in treatment: 1235897 
    INFO  @ Tue, 30 May 2017 23:22:39: #1 finished! 
    INFO  @ Tue, 30 May 2017 23:22:39: #2 Build Peak Model... 
    INFO  @ Tue, 30 May 2017 23:22:39: #2 Skipped... 
    INFO  @ Tue, 30 May 2017 23:22:39: #2 Use 229 as fragment length 
    INFO  @ Tue, 30 May 2017 23:22:39: #3 Call peaks... 
    INFO  @ Tue, 30 May 2017 23:22:39: #3 Call broad peaks with given level1 -log10qvalue cutoff and level2: 1.301030, 1.000000... 
    INFO  @ Tue, 30 May 2017 23:22:39: #3 Pre-compute pvalue-qvalue table... 
    INFO  @ Tue, 30 May 2017 23:22:48: #3 Call peaks for each chromosome... 
    INFO  @ Tue, 30 May 2017 23:22:49: #4 Write output xls file... /input_dir/mesc_atac_5M_MACS_peaks.xls 
    INFO  @ Tue, 30 May 2017 23:22:50: #4 Write broad peak in broadPeak format file... /input_dir/mesc_atac_5M_MACS_peaks.broadPeak 
    INFO  @ Tue, 30 May 2017 23:22:50: #4 Write broad peak in bed12/gappedPeak format file... /input_dir/mesc_atac_5M_MACS_peaks.gappedPeak 
    INFO  @ Tue, 30 May 2017 23:22:50: Done! 
    (py27) (py27) Command run:  /root/anaconda3/envs/py27/bin/nucleoatac run --bed /input_dir/mesc_atac_5M_MACS_merged.bed --bam /output_dir/mesc_atac_5M_filt_sorted.bam --fasta /input_dir/mm9.fa --out /input_dir/mesc_atac_5M_nucATAC --cores 8
    nucleoatac version 0.3.4
    start run at: 2017-05-30 23:22
    ---------Step1: Computing Occupancy and Nucleosomal Insert Distribution---------
    Making figure
    ---------Step2: Processing Vplot------------------------------------------------
    ---------Step3: Obtaining nucleosome signal and calling positions---------------
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: divide by zero encountered in log
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:120: RuntimeWarning: invalid value encountered in multiply
      nuc_lik = np.sum(np.log(nuc_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: divide by zero encountered in log
      null_lik = np.sum(np.log(null_model) * mat)
    /root/anaconda3/envs/py27/lib/python2.7/site-packages/nucleoatac/NucleosomeCalling.py:121: RuntimeWarning: invalid value encountered in multiply
      null_lik = np.sum(np.log(null_model) * mat)
    ---------Step4: Making combined nucleosome position map ------------------------
    ---------Step5: Calling NFR positions-------------------------------------------
    end run at: 2017-05-30 23:28
    (py27) 
