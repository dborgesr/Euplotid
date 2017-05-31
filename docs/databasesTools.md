
## Publicly available databases are <br>expanding and improving at incredible rate
* Many publicly available databases exist, many more being added
* Aid in asking questions for many things from expression of RNA and/or protein to identification of relevant literature
* Can range from easy to use to requiring a strong programming background
* Lets see what questions and how these databases can help you with


## What can they can do for you?
* <b>Why</b> What kind of questions can this resource help you answer?
* <b>What</b>: What is the resource made up of, what data does it contain and expose to you in a digested way?
* <b>Where</b>: Where to access the resource
* <b>How</b>: How to use the resource, I’ll go through a couple of example play questions which can be answered in each resource

## Genotype-Tissue Expression (GTEx)
* <b>Why</b>:  
    * Get estimate of RNA expression of gene/isoform in a tissue of interest
    * See if any non-coding mutations associate with change in expression
* <b>What</b>:
    * Expression across ~53 tissues and 544 Donors
    * Latest eQTL analysis for all tissues
    * Intuitive and reactive browser for eQTLs in genome and gene, isoform, and exon level expression
* <b>Where</b>: 
    * <a href="http://www.gtexportal.org/">http://www.gtexportal.org/</a>
* <b>How</b>: 
    * <a href="https://www.gtexportal.org/home/documentationPage">https://www.gtexportal.org/home/documentationPage</a>

## Nextprot
* <b>Why</b>:  
    * Need to get a better understanding of a protein
    * What and where is it expressed?
    * What does it look like?
    * What might it be doing?
* <b>What</b>:
    * Amalgamation of many data dumps, such as Human Proteome
    * Interface can be taxing on computer and may slow/crash
    * Very good way to get an overview of what’s known about a protein’s function
* <b>Where</b>: 
    * <a href="http://www.nextprot.org/">http://www.nextprot.org/</a>
* <b>How</b>: 
    * <a href="https://www.nextprot.org/help/simple-search">https://www.nextprot.org/help/simple-search</a>

## Uniprot
* <b>Why</b>:  
    * Need to fetch protein sequences
    * What’s its general structure characteristics?
    * Want to find homologous proteins
    * What does the protein's 3D structure look like?
    * What Post translational modifications does it undergo?
* <b>What</b>:
    * Lots of mass spec data, phylogenetic relationships, isoform and gene level data
    * Many phylogenetic databases, i've had good luck with [InParanoid](http://inparanoid.sbc.su.se/) and [EggNog](http://eggnogdb.embl.de)
* <b>Where</b>: 
    * <a href="http://www.uniprot.org/">http://www.uniprot.org/</a>
* <b>How</b>: 
    * Just search for your protein and scroll down!

## StringDB
* <b>Why</b>:  
    * What proteins have been shown to associate with my protein?
    * What gets pulled down with it?
    * Has it been investigated in the context of another protein?
* <b>What</b>:
    * A curated database of protein to protein interactions
    * Includes natural language processing of pubmed to find co-mentioned proteins and/or concepts
    * Includes and uses knowledge across species
* <b>Where</b>: 
    * <a href="https://string-db.org/">https://string-db.org/</a>
* <b>How</b>: 
    * <a href="https://string-db.org/cgi/help.pl">https://string-db.org/cgi/help.pl</a>

## Short Read Archive (SRA)
* <b>Why</b>:  
    * Need to fetch sequencing data from a paper
* <b>What</b>:
    * Most papers that generate any sequencing data are required to deposit it here. It is stored as SRA files.
    * Petabases of information
    <img src="https://trace.ncbi.nlm.nih.gov/Traces/sra/i/g.png" style="width: 500px;">
* <b>Where</b>: 
    * <a href="https://trace.ncbi.nlm.nih.gov/Traces/sra/">https://trace.ncbi.nlm.nih.gov/Traces/sra/</a>
* <b>How</b>: 
    * [this pipeline](getFastqReads.md) shows how to pull and process reads from SRA

## Genome Expression Omnibus
* <b>Why</b>:  
    * Need to fetch quantified values such as, Wig, FPKM, BedPeaks...
* <b>What</b>:
    * Data dump of usually a more quantified version of SRA
    * Things such as FPKM, bed peaks
    * Also a big repository for microarray data
* <b>Where</b>: 
    * <a href="https://www.ncbi.nlm.nih.gov/geo/">https://www.ncbi.nlm.nih.gov/geo/</a>
* <b>How</b>: 
    * <a href="https://www.ncbi.nlm.nih.gov/geo/info/download.html">https://www.ncbi.nlm.nih.gov/geo/info/download.html</a>

## Meta
* <b>Why</b>:  
    * Literature search
* <b>What</b>:
    * AI organized pubmed
    * Can follow and combine "streams" into customized feeds
    * Can quickly find influential papers
    * A good starting point for new genes
* <b>Where</b>: 
    * <a href="http://meta.science/">http://meta.science/</a>
* <b>How</b>: 
    * <a href="http://www.editage.com/insights/using-meta-science-to-streamline-researcher-workflow-systems-a-step-by-step-guide">How to use Meta

## PRoteomics IDEntifications (PRIDE)
* <b>Why</b>:  
    * You have found a particular protein or peptide and want the raw data
    * Need to rerun MS/MS analysis to quantify a different set of peptides
* <b>What</b>:
    * Database dump of MS/MS spectra for a large volume of experiments, similar to SRA for sequencing data
* <b>Where</b>: 
    * <a href="https://www.ebi.ac.uk/pride/archive/">https://www.ebi.ac.uk/pride/archive//</a>
* <b>How</b>: 

## Harmonizome
* <b>Why</b>:  
    * Need to get a general idea of a gene
    * Curious about what medical issues have been associated with it
    * What complexes has it been found to be a part of?
    * What drugs have been shown to affect its behaviour?
* <b>What</b>:
    * "Search for genes or proteins and their functional terms extracted and organized from over a hundred publicly available resources"
    * X = gene, Y = database
    * Includes clinical information and pharmacological information related to the query gene as well
    * Awesome mobile interface!
    * http://amp.pharm.mssm.edu/Harmonizome/about
* <b>Where</b>: 
    * <a href="http://amp.pharm.mssm.edu/Harmonizome/">http://amp.pharm.mssm.edu/Harmonizome/</a>
* <b>How</b>: 
    * Simply search for your gene and click through the available databases for that given gene (LOTS of stuff!)

## Encode
* <b>Why</b>:  
    * Need to find Chip-Seq, RNA-seq, DNA-seq/ATAC-seq, etc genome track to look at a particular locus in a particular cell type
* <b>What</b>:
    * Database dump of mostly Chip-Seq, RNA-seq, and Chromatin Accessibility processed files
* <b>Where</b>: 
    * <a href="https://www.encodeproject.org/">https://www.encodeproject.org/</a>
* <b>How</b>: 
    * <a href="https://www.encodeproject.org/help/getting-started/">https://www.encodeproject.org/help/getting-started/</a>

## UCSC Xena Browser
* <b>Why</b>:  
    * Need to look at expression differences between publicly available cohorts, such as cancer patients or geographic populations
    * Dont want to process or make the raw graphs yourself


```python
from IPython.display import YouTubeVideo
YouTubeVideo("TSNc-EDjix4", start=0, autoplay=1, theme="light", color="red")
```





        <iframe
            width="400"
            height="300"
            src="https://www.youtube.com/embed/TSNc-EDjix4?color=red&start=0&theme=light&autoplay=1"
            frameborder="0"
            allowfullscreen
        ></iframe>
        



* <b>What</b>:
    * Amalgamation of many data dumps, such as Human Proteome
    * Interface can be taxing on computer and may slow/crash
    * Very good way to get an overview of what’s known about a protein’s function
* <b>Where</b>: 
    * <a href="http://www.nextprot.org/">http://www.nextprot.org/</a>
* <b>How</b>: 


```python
YouTubeVideo("go38U6iLjsw", start=0, autoplay=1, theme="light", color="red")
```





        <iframe
            width="400"
            height="300"
            src="https://www.youtube.com/embed/go38U6iLjsw?color=red&start=0&theme=light&autoplay=1"
            frameborder="0"
            allowfullscreen
        ></iframe>
        



## Ensembl
* <b>Why</b>:  
    * How far back is your gene conserved?
    * Is it lost in some species?
    * What's its synteny in humans, other species?
    * Is the protein conserved? How far?
    * What different isoforms exist?
* <b>What</b>:
    * HUGE number of different organisms and genome data
    * Much integrated knowledge surrounding phylogenetics and transcriptomics
    * Data which is accessible is enormous
    * Built in genome browser which can be used to compare a huge number of tracks
    * Can compare synteny in genome browser view across many species, being able to check for things such as conserved cis regulatory elements and/or synteny
* <b>Where</b>: 
    * <a href="http://www.ensembl.org/index.html">http://www.ensembl.org/index.html</a>
* <b>How</b>: 
    * <a href="http://www.ensembl.org/info/website/index.html">http://www.ensembl.org/info/website/index.html</a>

## MobiDB
* <b>Why</b>:  
    * Does your protein have high mobility domains (does it wiggle lots)?
* <b>What</b>:
    * A combination of curated (direct provided by <a href="http://www.disprot.org/">DisProt</a>), and indirect PDB NMR/Xray sources amalgamated into a more digestible and encompassing database
* <b>Where</b>: 
    * <a href="http://mobidb.bio.unipd.it/about">http://mobidb.bio.unipd.it/about</a>
* <b>How</b>: 
    * <a href="http://mobidb.bio.unipd.it/entries/P49711">CTCF</a>

## The Human Protein Atlas
* <b>Why</b>:  
    * Where is your protein expressed?
    * Is it upregulated in specific cancers?
    * what does the histology look like?
* <b>What</b>:
    * LOTS of staining of tissue slides. 
    * Integrates well with protein levels from mass proteomics sources
* <b>Where</b>: 
    * <a href="http://www.proteinatlas.org/">http://www.proteinatlas.org/</a>
* <b>How</b>: 
    * <a href="http://www.proteinatlas.org/ENSG00000102974-CTCF/cell">CTCF</a>
