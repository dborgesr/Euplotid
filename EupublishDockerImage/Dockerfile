# Eupublish image to write and publish ipynb into pdf
FROM debian:stretch-slim
MAINTAINER Diego Borges <dborgesrmit@gmail.com>
USER root

# Install few tools using apt-get required to get conda package manager up and running 
RUN apt-get update && apt-get install -y --no-install-recommends \
  software-properties-common \
  wget \
  bzip2 \
  vim \
  libboost-dev \
  curl \
  zip \
  gnupg2 \
  unzip \
  apt-transport-https \
  ca-certificates \
  texlive-latex-base \
  texlive-fonts-recommended \
  texlive-fonts-extra \
  texlive-latex-extra \
  texlive-generic-recommended \
  && apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
	
#Download, install, and update Miniconda3
ENV CONDA_VERSION="4.5.12"  
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh && bash Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -b && rm Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh && \
    export PATH=/root/miniconda3/bin:$PATH  && \
    conda config --set show_channel_urls True && \
	conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y pip \
	    setuptools \
		nb_conda \
		conda-build \
		ipywidgets \
		nodejs \
		ipypublish \
		latexmk \
		tornado \
	&& \
	conda clean -tipy && \
	if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi
ENV PATH=/root/miniconda3/bin:$PATH	
	
RUN pip install --ignore-installed --no-cache-dir \
		setuptools \
		ipywidgets \
		h5py \
		argparse \
		jupyter_contrib_nbextensions \
		jupyterlab \
		tornado \
		jupytext

#Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.open_browser/c c.NotebookApp.open_browser = False" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.ip/c c.NotebookApp.ip = '*'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.notebook_dir/c c.NotebookApp.notebook_dir = '/root/Euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" /root/.jupyter/jupyter_notebook_config.py 

#set bash as default shell
ENV SHELL=/bin/bash

RUN apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Get Euplotid
#RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git
RUN mkdir /root/Euplotid && mkdir /root/Euplotid/manuscript
COPY ./manuscript/ /root/Euplotid/manuscript/
WORKDIR /root/Euplotid

EXPOSE 8898
CMD ["jupyter", "lab", "--port=8898", "--no-browser", "--allow-root", "--ip=0.0.0.0"]

