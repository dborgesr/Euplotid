
<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/Title_slide.png" style="width: 1000px;"></div>

[![Introduction video](https://img.youtube.com/vi/VID/0.jpg)](https://www.youtube.com/embed/wNuoL09rqtc)

<div style="text-align:center"><img src="https://raw.githubusercontent.com/dborgesr/Euplotid/gh-pages/web_euplotid/graphical_abstract.png" style="width: 5=800px;"></div>

##Pull and run Euplotid Docker images

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