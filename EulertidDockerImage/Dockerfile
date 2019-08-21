# Eulertid Anaconda3 image to get started w/ programming
FROM debian:stretch-slim
MAINTAINER Diego Borges <dborgesr@mit.edu>
USER root

# Install few tools using apt-get required to get conda package manager up and running 
RUN apt-get update && apt-get install -y --no-install-recommends \
  software-properties-common \
  wget \
  bzip2 \
#  build-essential \
  vim \
  libboost-dev \
#  libcairo2-dev \
  curl \
  zip \
  gnupg2 \
  unzip \
  apt-transport-https \
  ca-certificates \
  && apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
	
#Download, install, and update Miniconda3
ENV CONDA_VERSION="4.7.10"  
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh && bash Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh -b && rm Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh && \
    export PATH=/root/miniconda3/bin:$PATH  && \
    conda config --set show_channel_urls True && \
	conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda update conda && \
    conda update --all --yes && \
    conda install -y pip \
	    setuptools \
		nb_conda \
#		libgcc \
		conda-build \
		git \
		plotly \
		ipywidgets \
		nodejs \
	&& \
#	conda update -y --all && \
	conda clean -tipy && \
	if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi
ENV PATH=/root/miniconda3/bin:$PATH	
	
#install docker
#RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - \
#	&& add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable" \
#	&& apt-get update \
#	&& apt-get install -y docker-ce

#Get jupyter notebooks from yamaton for Euler problems
RUN cd /root && git clone https://github.com/yamaton/project-euler-jupyter.git

RUN pip install  --no-cache-dir \
		setuptools \
		ipywidgets \
		h5py \
		argparse \
		jupyter_contrib_nbextensions \
#		jupyter-c-kernel \ 
		jupyterlab 
#		google-cloud-storage
#		&& pip freeze --local | grep -v '^\-e' | cut -d = -f 1  | xargs -n1 pip install -U
	
#Enable C kernel after installed using pip above (jupyter-c-kernel)
#RUN install_c_kernel	

#Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.open_browser/c c.NotebookApp.open_browser = False" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.ip/c c.NotebookApp.ip = '*'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.notebook_dir/c c.NotebookApp.notebook_dir = '/root/project-euler-jupyter'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" /root/.jupyter/jupyter_notebook_config.py 
 #       && jupyter nbextension enable --py widgetsnbextension \
 #       && jupyter contrib nbextension install --system 
		
#Enable google drive in jupyterlab
#RUN jupyter serverextension enable --py jupyterlab --sys-prefix	&& \
#	jupyter labextension install @jupyterlab/google-drive
	
#set bash as default shell
ENV SHELL=/bin/bash

#run a closing update and clean
#RUN conda update -y --all \
#	&& conda clean -y -t -i -l -a
RUN apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Get Euplotid
#RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git
RUN mkdir /root/Euplotid

# Set Euler problems as working directory
WORKDIR /root/project-euler-jupyter

EXPOSE 8896
CMD ["jupyter", "lab", "--port=8896", "--no-browser", "--allow-root", "--ip=0.0.0.0"]

