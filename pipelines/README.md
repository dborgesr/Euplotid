
# Pipeline Directory
## Pipelines to
1. Hello world intro to programming and Jupyter's capabilities [helloWorld](helloWorld.ipynb) O*.
2. Databases and good tools to crawl the internet for interesting datasets and hypothesis [databasesTools](databasesTools.ipynb) O*.
3. Fetch any type of sequencing data from SRA [getFastqReads](getFastqReads.ipynb) O
4. QC, trim, and filter sequencing reads [fq2preppedReads](fq2preppedReads.ipynb) O
5. Call peaks from Chip-Seq and Chromatin Accessibility reads [fq2peaks](fq2peaks.ipynb) O
6. Call normalized interactions from ChiA-PET reads [fq2ChIAInts](fq2ChIAInts.ipynb) O
7. Call normalized Interactions from HiC reads [fq2HiCInts](fq2HiCInts.ipynb) O
8. Call normalized interactions from Hi-ChIP reads [fq2HiChIPInts](fq2HiChIPInts.ipynb) O
9. Call normalized interactions from DNAse-HiC reads [fq2DNAseHiCInts](fq2DNAseHiCInts.ipynb) O
10. Call normalized expression and counts from RNA-Seq reads [fq2countsFPKM](fq2countsFPKM.ipynb) O
11. Call differentially expressed genes from RNA-seq counts [countsFPKM2DiffExp](countsFPKM2DiffExp.ipynb) O
12. Call normalized counts, miRNA promoters, and nascent transcripts from Gro-Seq reads [fq2GroRPKM](fq2GroRPKM.ipynb) O
13. Call normalized interactions from 4C [fq24CInts](fq24CInts.ipynb) O
14. Build, annotate and add INs to global graph for a given cell state using DNA-DNA interactions, Chromatin Accesibility, and FPKM [addINs](addINs.ipynb) *
15. View current built and annotated INs for all cell types [viewINs](viewINs.ipynb) O*.
16. Search for and/or manipulate annotation and other data available to euplotid [annotationManagement](annotationManagement.ipynb) O*.
17. Description of default software packages and images installed, how to get new ones, and which ones are currently installed. [packageManagement](packageManagement.ipynb) O*.

 [O] = Megatid compatible
 [*] = Euplotid compatible
 [.] = Minitid compatible

# 1. [helloWorld](helloWorld.ipynb)
Hello world intro to programming, ipython, and Euplotid

# 2. [databasesTools](databasesTools.ipynb)
Databases and good tools to crawl the internet for interesting datasets and hypothesis. Some examples include GTeX, uniprot, SRA, GEO, etc, check them out!!

# 3. [getFastqReads](getFastqReads.ipynb)
Allows you to use Tony to find local fastq.gz files OR provide an SRA number to pull from

# 4. [fq2preppedReads](fq2preppedReads.ipynb)
Take fq.gz reads and QC them using FastQC checking for over-represented sequences potentially indicating adapter contamination. Then use cutadapt and sickle to filter and remove adapters. Can also use trimmomatic for flexible trimming. 

# 5. [fq2peaks](fq2peaks.ipynb)
Take fq.gz align it using bowtie2 to the genome. Then using Homer software pick the type of peak (histone, chip-seq, dnase, etc) and chug through to get bed files of peaks. Can also use MACS2 w/ specific analysis parameters to deal with different types of peak finding problems.

# 6. [fq2ChIAInts](fq2ChIAInts.ipynb)
Take fq.gz reads, prep them by removing bridge adapters (can deal with either bridges), align, find interactions, normalize, and spit into cooler format for later viewing. Can perform analysis using either Origami or ChiA-PET2

# 7. [fq2HiCInts](fq2HiCInts.ipynb)
Take fq.gz reads and chug them through HiCPro w/ tuned relevant parameters. In the end spits out a cooler file which can be loaded for further visualization.

# 8. [fq2HiChIPInts](fq2HiChIPInts.ipynb)
Take fq.gz reads and chug them through customized Origami pipeline and customized HiCPro pipeline. In the end spits out a cooler file which can be loaded for further visualization.

# 9. [fq2DNAseHiCInts](fq2DNAseHiCInts.ipynb)
Take fq.gz reads and chug them through HiCPro pipeline. In the end spits out a cooler file which can be loaded for further visualization.

# 10. [fq2countsFPKM](fq2countsFPKM.ipynb)
Take fq.gz reads and chug them through STAR aligner and then RSEM pipeline. In the end spits out a counts vs transcripts matrix and a normalized transcript/gene FPKM matrix.

# 11. [countsFPKM2DiffExp](countsFPKM2DiffExp.ipynb)
Take RNA-seq count and FPKM matrix and run any one of many R packages (DESeq2,DESeq,EBSeq,edgeR...) to call differentially expressed genes. Plotting and interactive visualization of results included

# 12. [fq2GroRPKM](fq2GroRPKM.ipynb)
Take fq.gz reads and align them using bowtie2 then find nascent transcripts using FStitch and miRNA promoters using mirSTP

# 13. [fq24CInts](fq24CInts.ipynb)
Take fq.gz reads and align them using bowtie2. Chug them through HiCPro and/or custom pipeline to get cooler file

# 14.  [addINs](addINs.ipynb)
Build, annotate and add INs to global graph for a given cell state using DNA-DNA interactions, Chromatin Accesibility, and FPKM.

# 15. [viewINs](viewINs.ipynb)
View current built and annotated INs for all cell types

# 16. [annotationManagement](annotationManagement.ipynb)
Search for and/or manipulate annotation and other data available to euplotid

# 17. [packageManagement](packageManagement.ipynb)
Description of default image and the software packages that are installed, also how to get new packages, and how to export environment in yaml file for others to replicate analysis.
