
<h1><center>Euplotid</center></h1>
<h2><center> Defining mappings between Cis-Regulatory Elements <br> and the genes they may be regulating using DNA interactions</center></h2>
Composed of a set of constantly updating pipelines running on Docker images which enable a user to build and annotate Insulated Neighborhoods(INs)starting from raw sequencing reads of DNA-interactions, Chromatin Accessibility, and RNA-sequencing. Reads are quantified using the latest computational tools (RSEM, STAR, Origami, HiCPro) and the results are normalized, quality-checked and stored. INs are then built using a Louvain based graph partioning algorithm. CREs (as defined by Chromatin Accessibility) are then assigned to TSSs based on inclusion within the same IN and TF binding at each peak is annotated using a trained neural network architecture. The annotated output can be visualized in a variety of 1D, 2D and 3D ways.

Composed of 3 Docker images:
* [Megatid](http://localhost:8891) process sequencing data into quantified values (FPKM,peaks,etc)
* [Euplotid](http://localhost:8890) build and visualize Insulated Neighborhoods and learn/predict TFs bound at Cis-regulatory Elements
* [Minitid](http://localhost:8892) visualize and interact with built and annotated Insulated Neighborhoods

The pipelines available and their capabilities are described in [pipeline README](pipelines/README.md) which tells you which image you need. 

Each Docker image has different capabilities (packages installed in each described in [packageManagement](pipelines/packageManagement.ipynb) 
