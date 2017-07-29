#Nanotid image to take in FPKM, normalized DNA-DNA interactions, ChiP-Seq, Chromatin Accessibility and build, annotate, and visualize INs
#x86
#FROM debian:jessie
#FROM alpine:3.6

#ARM
#FROM resin/rpi-raspbian:jessie
#FROM resin/aarch64-debian
#FROM resin/raspberry-pi-alpine
#FROM resin/aarch64-alpine

#ARM on x86
FROM resin/armv7hf-debian-qemu
#FROM amtypaldos/armv7hf-alpine-qemu
RUN [ "cross-build-start" ]

#set CUDA paths
ENV PATH /usr/local/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH /usr/lib/aarch64-linux-gnu:/usr/local/cuda-8.0/targets/aarch64-linux/lib/:/usr/lib/aarch64-linux-gnu/tegra/"

USER root

WORKDIR /root

ENV APK_PACKAGES=" \
	tk-dev \
	libpng-dev \
	gfortran \
	wget \
	git \
	vim \
	libxml2-dev \
	ca-certificates \
	build-base \
	bzip2-dev \
	libc6-compat \
	gdbm-dev \
#	hdf5-dev \
	jpeg-dev \
	lapack-dev \
	ncurses-dev \
	sqlite-dev \
	openssl-dev \
	zlib-dev \
	xz-dev \
	linux-headers \
	curl \
	python3 \
	python3-dev \
	readline-dev \
	patch \
	bash \
	cmake \
	i2c-tools"

ENV APT_PACKAGES=" \
	tk-dev \
    libpng-dev \
	gfortran \
	wget \
	git \
vim \
	libxml2-dev \
	ca-certificates \
	build-essential \
	libbz2-dev \
	libc6-dev \
	libgdbm-dev \
	libhdf5-dev \
	libjpeg-dev \
	liblapack-dev \
	libncurses-dev \
#	libncursesw5-dev \
	libsqlite3-dev \
	libssl-dev \
	zlib1g-dev \
	curl \
#	python-smbus \
	i2c-tools \
	liblzma-dev \
	xz-utils \
	libreadline-dev \
	cmake \
	python3 \
	python3-dev \
	python3-pip"

RUN apt-get update && apt-get upgrade && apt-get install -y $APT_PACKAGES \
  && if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi \
  && pip install --upgrade pip setuptools 

# Install packages necessary for compiling python using Alpine package manager, few extra commands for linking libraries
#RUN apk update  \ 
#  # replacing default repositories with edge ones
#  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/testing" > /etc/apk/repositories \
#  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories \
#  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/main" >> /etc/apk/repositories \ 
#  && apk add --no-cache $APK_PACKAGES  \
#  && python3 -m ensurepip \
#  && rm -r /usr/lib/python*/ensurepip \
#  && if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi \
#  && pip install --upgrade pip setuptools \
#  && rm -r /root/.cache \
#  && ln -s /usr/include/locale.h /usr/include/xlocale.h
  
#RUN curl -s -L https://github.com/jjhelmus/berryconda/releases/download/v1.0.0/Berryconda3-1.0.0-Linux-armv7l.sh > miniconda.sh && \
#    bash miniconda.sh -b -p /opt/conda && \
#    rm miniconda.sh && \
#    export PATH=/opt/conda/bin:$PATH  && \
#    conda config --set show_channel_urls True && \
#    conda config --add channels rpi && \
#    conda config --add channels conda-forge && \
#    conda update conda && \
#    conda install python=3.6 && \
#    conda update --all --yes && \
#    conda install conda-build && \
#    conda install anaconda-client && \
#    conda install -y pip \
#	    setuptools \
#		readline \
	#	jupyter_core \
#		pandas \
#    	scipy \
#		numpy \
#    	networkx \
    #	plotly \
	#	sense-hat \
	#	myvariant \
	#	mygene \
	#	pysam \
	#	joblib \
	#	biopython \
	#	pyliftover \
	#	nbpresent \
	#	python-louvain \
#	&& conda clean -tipy
#ENV PATH=/opt/conda/bin:$PATH

RUN pip install \
	setuptools \
	gnureadline \
	jupyter \
	plotly \
	sense-hat \
	myvariant \
	mygene \
	pysam \
    pandas \
    joblib \
    biopython \
    pyliftover \
    scipy \
    numpy \
    nbpresent \
    networkx \
    python-louvain && \
	pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U

#Get pytorch machine learning
#ENV CMAKE_PREFIX_PATH=$CONDA_DIR
#RUN cd /root && git clone https://github.com/pytorch/pytorch.git && \
#	cd pytorch && \
#	pip install -r requirements.txt && \
#	python3 setup.py install
	
#RUN conda install -c soumith magma-cuda80 # or magma-cuda75 if CUDA 7.5

#Get stuff for controlling real life
RUN cd /root && git clone git://git.drogon.net/wiringPi && cd wiringPi && ./build
RUN cd /root && git clone https://github.com/timleland/rfoutlet.git /var/www/rfoutlet
RUN cd /root && git clone https://github.com/RPi-Distro/RTIMULib.git && cd RTIMULib/Linux/python && python3 setup.py build && python3 setup.py install

#Get Euplotid  
RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git
# Set pipelines as workdir
WORKDIR /root/Euplotid

#Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.open_browser/c c.NotebookApp.open_browser = False" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.ip/c c.NotebookApp.ip = '*'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.notebook_dir/c c.NotebookApp.notebook_dir = '/root/Euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" /root/.jupyter/jupyter_notebook_config.py

EXPOSE 80
CMD ["jupyter", "notebook", "--allow-root", "--port=80", "--no-browser", "--ip=0.0.0.0"]
