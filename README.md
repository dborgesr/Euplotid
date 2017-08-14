

<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/lab_meeting_slides_Title_slide.png" style="width: 1000px;"></div>

# Abstract
### http://dborgesr.github.io/Euplotid/
### http://www.biorxiv.org/content/early/2017/08/03/170159
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Euplotid is composed of a set of constantly evolving bioinformatic pipelines encapsulated and running in Docker containers enabling a user to build and annotate the local regulatory structure of every gene starting from raw sequencing reads of DNA-interactions, chromatin accessibility, and RNA-sequencing. Reads are quantified using the latest computational tools and the results are normalized, quality-checked, and stored. The local regulatory neighborhood of each gene is built using a Louvain based graph partitioning algorithm parameterized by the chromatin extrusion model and CTCF-CTCF interactions. Cis-Regulatory Elements are defined using chromatin accessibility peaks which are then mapped to Transcription Start Sites based on inclusion within the same neighborhood. Convolutional Neural Networks are combined with Long-Short Term Memory in order to provide a statistical model mimicking transcription factor binding, one neural network for each transcription factor in the genome is trained on all available Chip-Seq and SELEX data, learning what pattern of DNA oligonucleotides the factor binds. The neural networks are then merged and trained on chromatin accessibility data, building a rationally designed neural network architecture capable of predicting chromatin accessibility. Transcription factor binding and identity at each peak is annotated using this trained neural network architecture. By in-silico mutating and re-applying the neural network we are able to gauge the impact of a transition mutation on the binding of any human transcription factor. The annotated output can be visualized in a variety of 1D, 2D and 3D ways overlaid with existing bodies of knowledge, such as GWAS results. Once a particular CRE of interest has been identified by a biologist the difficulty of a <a href="https://blog.benchling.com/base-editor/"> Base Editor 2 (BE2) </a> mediated transition mutation can be quantitatively assessed and induced in a model organism.

<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/lab_meeting_slides_graphical_abstract.png" style="width: 5=800px;"></div>

## Get Docker
[**INSTALL DOCKER HERE**](https://www.docker.com/community-edition#/download")

The pipelines available and their capabilities are described in [Methods](docs/07_Methods.md)  which helps you pick the right of 3 Docker images.

Remember to define the correct directories when running the Docker image depending on your local OS machine, ex:

Whitehead (Linux):
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
