FROM ubuntu:latest
USER root

# Install few tools using apt-get required to get conda package manager up and running 
RUN apt-get update && apt-get install -y --no-install-recommends \
  software-properties-common \
  wget \
  bzip2 \
  vim \
  libboost-dev \
  libcairo2-dev \
  curl \
  zip \
  unzip \
  apt-transport-https \
  ca-certificates \
  curl \
  && apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
	
#Download, install, and update Anaconda3
ENV CONDA_VERSION="5.1.0"  
RUN wget https://repo.continuum.io/archive/Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh && bash Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh -b && rm Anaconda3-${CONDA_VERSION}-Linux-x86_64.sh && \
    export PATH=/root/anaconda3/bin:$PATH  && \
    conda config --set show_channel_urls True && \
    conda update conda && \
    conda update --all --yes && \
    conda install -y pip \
	    setuptools \
		nb_conda \
		libgcc \
		conda-build \
		git \
		plotly \
		ipywidgets \
		nodejs \
	&& \
	conda clean -tipy && \
	if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi
ENV PATH=/root/anaconda3/bin:$PATH

RUN pip install --ignore-installed --no-cache-dir \
		setuptools \
		ipywidgets \
		h5py \
		argparse \
		jupyter_contrib_nbextensions \
		jupyter-c-kernel \ 
		jupyterlab 	

#install bash kernel
RUN pip install bash_kernel && \
	python -m bash_kernel.install

#Enable C kernel after installed using pip above (jupyter-c-kernel)
RUN install_c_kernel	
	
#Configure Jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.open_browser/c c.NotebookApp.open_browser = False" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.ip/c c.NotebookApp.ip = '*'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.notebook_dir/c c.NotebookApp.notebook_dir = '/root/Euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" /root/.jupyter/jupyter_notebook_config.py \
        && jupyter nbextension enable --py widgetsnbextension \
        && jupyter contrib nbextension install --system \
        && jupyter nbextension install nbpresent --py --overwrite --system \
		&& jupyter nbextension enable nbpresent --py --system \
		&& jupyter serverextension enable nbpresent --py --system
		
#Enable google drive in jupyterlab
RUN jupyter serverextension enable --py jupyterlab --sys-prefix	&& \
	jupyter labextension install @jupyterlab/google-drive && \
	jupyter labextension install @jupyterlab/plotly-extension

#set bash as default shell
ENV SHELL=/bin/bash


#Get Euplotid
RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git

# Set Euplotid as working directory
WORKDIR /root/Euplotid

EXPOSE 8889
CMD ["jupyter", "lab", "--port=8889", "--no-browser", "--allow-root", "--ip=0.0.0.0"]