# Megatid image to suck in fastq and produce DNA-DNA interactions, chromatin accessibility, and RNA expression
FROM debian:jessie
MAINTAINER Diego Borges <dborgesr@mit.edu>
USER root

# Install few tools using apt-get required to get conda package manager up and running 
RUN apt-get update && apt-get install -y --no-install-recommends \
  software-properties-common \
  wget \
  bzip2 \
  build-essential \
  vim \
  libboost-dev \
  libcairo2-dev \
#  libzmq3-dev #iperl\
#  libcairo2-dev #iperl\
  && apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

ENV CONDA_VERSION="4.2.0"  
#Download, install, and update Anaconda3
RUN cd /root && wget https://repo.continuum.io/archive/Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh && bash Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh -b && rm Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh
ENV PATH=/root/anaconda3/bin:$PATH
RUN conda update conda anaconda jupyter notebook
#RUN conda config --set channel_priority false
RUN conda config --add channels conda-forge 
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda install -y \
	nb_conda \
	fortran-magic \
	make \
	gcc \
	conda-build \
	git \
	ipyparallel \
	plotly \
	ipywidgets \
	ipywidgets \
	r-essentials \
	r-irkernel \
	r-ggplot2 \
	r-rcolorbrewer \
	autoconf \
	bamtools \
	bowtie \
	bowtie2 \
	htslib \
	cython \
	samtools \
	bedops \
	bwa \
	star \
	sra-tools \
	sickle \
	homer \
	trimmomatic \
	fastqc \
	bedtools \
	bamtools \
	perl-bioperl \
	perl-app-cpanminus \
	libmagic \
	wkhtmltopdf \
	&& conda clean -y -t -i -l -a

#Install python 2.7 kernel and environment
RUN conda create -y -n py27 python=2.7 \
	ipykernel \
	cooler \
	bx-python \
	pysam \
	nucleoatac \
	macs2 \
	&& conda clean -y -t -i -l -a

# Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py && \
	sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" /root/.jupyter/jupyter_notebook_config.py 

#install bash kernel
RUN pip install bash_kernel
RUN python -m bash_kernel.install

#Install perl package manager and perl kernel takes super long and currently broken
#RUN cpanm --notest Reply::Plugin::IPerl

#Install C kernel
RUN cd /root && wget -O - https://raw.githubusercontent.com/brendan-rius/jupyter-c-kernel/master/install.sh | sh
	
#Get Chrom3D
RUN cd /root && git clone https://github.com/CollasLab/Chrom3D.git && cd ./Chrom3D && make

#Get ChiAPet2
RUN cd /root && git clone https://github.com/GuipengLi/ChIA-PET2.git && cd ./ChIA-PET2  && chmod 755 ./bin/ChIA-PET2 && make
ENV PATH=~/ChIA-PET2/bin/:$PATH

#Get Origami compilation broken
#RUN cd /root && git clone https://github.com/danielsday/origami.git && cd origami && autoconf && ./configure --prefix=/root/origami CPPFLAGS=-I=/root/bamtools/include

#Get Euplotid
RUN cd /root && git clone https://github.com/dborgesr/Euplotid.git

#run a closing update
RUN conda update -y --all \
	&& conda clean -y -t -i -l -a

# Set root as workdir
WORKDIR /root/Euplotid/pipelines

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents
# kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

EXPOSE 8891
CMD ["jupyter", "notebook", "--port=8891", "--no-browser", "--allow-root", "--ip=0.0.0.0"]