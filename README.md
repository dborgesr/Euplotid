
<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/Title_slide.png" style="width: 1000px;"></div>


<div style="position:relative;height:0;padding-bottom:56.25%"><iframe src="https://www.youtube.com/embed/wNuoL09rqtc" width="640" height="360" frameborder="0" style="position:absolute;width:100%;height:100%;left:0" allowfullscreen></iframe></div>


# Abstract
### http://dborgesr.github.io/Euplotid/
### https://www.biorxiv.org/content/early/2018/01/08/170159

Life as we know it has continued to shock and amaze us, consistently reminding us that truth is far stranger than fiction. Euplotid is a quantized geometric model of the eukaryotic cell, a first attempt at quantifying, using planck's constant geometric shape as its base, the incredible complexity that gives rise to a living cell. By beginning from the very bottom we are able to build the pieces, which when hierarchically and combinatorially combined, produce the emergent complex behavior that even a single celled organism can show. Euplotid is composed of a set of quantized geometric 3D building blocks and constantly evolving bioinformatic pipelines encapsulated and running in Docker containers enabling a user to build and annotate the local regulatory structure of every gene starting from raw sequencing reads of DNA-interactions, chromatin accessibility, and RNA-sequencing. Reads are quantified using the latest computational tools and the results are normalized, quality-checked, and stored. The local regulatory neighborhood of each gene is built using a Louvain based graph partitioning algorithm parameterized by the chromatin extrusion model and CTCF-CTCF interactions. Cis-Regulatory Elements are defined using chromatin accessibility peaks which are then mapped to Transcription Start Sites based on inclusion within the same neighborhood. Deep Neural Networks are trained in order to provide a statistical model mimicking transcription factor binding, giving the ability to identify all TFs within a given chromatin accessibility peak. By in-silico mutating and re-applying the neural network we are able to gauge the impact of a transition mutation on the binding of any transcription factor. The annotated output can be visualized in a variety of 1D, 2D, 3D and 4D ways overlaid with existing bodies of knowledge, such as GWAS results or PDB structures. Once a particular CRE of interest has been identified by a biologist the difficulty of a Base Editor mediated transition mutation can be quantitatively observed and induced in a model organism.

<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/graphical_abstract.png" style="width: 5=800px;"></div>

## Get Docker
[**INSTALL DOCKER HERE**](https://www.docker.com/community-edition#/download")

The pipelines available and their capabilities are described in [Methods](docs/07_Methods.md)  which helps you pick the right of 3 Docker images.

Remember to define the correct directories when running the Docker image depending on your local OS machine, ex:

Cluster (Linux):
~~~
    -v "/lab/solexa_public:/input_dir" \
    -v "/home/dborgesr/work_space/tmp:/tmp_dir" \
    -v "/home/dborgesr/work_space/out_dir:/output_dir" \
    -v "/home/dborgesr/work_space/annotation:/annotation_dir" \
~~~

## Pull and run your image!

* Megatid: process sequencing data into quantified values (FPKM,peaks,etc)
~~~ 
docker run --name megatid -p 8891:8891 -tid \
	-v "/your/input/directory:/input_dir" \
	-v "/your/temporary/directory/:/tmp_dir" \
	-v "/your/output/directory/:/output_dir" \
	-v "/your/annotation/directory/:/annotation_dir" \
	dborgesr/euplotid:megatid
~~~
* Euplotid: build and visualize INs and learn/predict TFs bound at CREs
~~~
docker run --name euplotid -p 8890:8890 -tid \
	-v "/your/input/directory:/input_dir" \
	-v "/your/temporary/directory/:/tmp_dir" \
	-v "/your/output/directory/:/output_dir" \
	-v "/your/annotation/directory/:/annotation_dir" \
	dborgesr/euplotid:euplotid
~~~
* Minitid: visualize and interact with built and annotated INs
~~~
docker run --name minitid -p 8892:8892 -tid \
	-v "/your/input/directory:/input_dir" \
	-v "/your/temporary/directory/:/tmp_dir" \
	-v "/your/output/directory/:/output_dir" \
	-v "/your/annotation/directory/:/annotation_dir" \
	dborgesr/euplotid:minitid
~~~
* Nanotid: ARM architecture image to build and visualize INs and learn/predict TFs bound at CREs
~~~
docker run --name nanotid -p 8893:8893 -tid \
	-v "/your/input/directory:/input_dir" \
	-v "/your/temporary/directory/:/tmp_dir" \
	-v "/your/output/directory/:/output_dir" \
	-v "/your/annotation/directory/:/annotation_dir" \
	dborgesr/euplotid:nanotid
~~~
Your Docker image should be running on your computer (local), now you can access and use it!

## Go to your image!
* Megatid <br> local:[http://localhost:8891](http://localhost:8891) Whitehead internal:[http://airstream:8891](http://airstream:8891)
* Euplotid <br> local:[http://localhost:8890](http://localhost:8890) Whitehead internal:[http://airstream:8890](http://airstream:8890)
* Minitid <br> local:[http://localhost:8892](http://localhost:8892) Whitehead internal:[http://airstream:8892](http://airstream:8892)
* Nanotid <br> local:[http://localhost:8893](http://localhost:8893) Whitehead internal:[http://airstream:8893](http://airstream:8893)

Each Docker image has different capabilities (packages installed in each Docker image are described in [packageManagement](docs/packageManagement.md))
