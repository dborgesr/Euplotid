# Euplotid image to take in FPKM, DNA-DNA interactions, ChiP-Seq, Chromatin Accessibility and calls and annotate/visualize INs
FROM rocm/pytorch:rocm3.0_ubuntu16.04_py3.6_pytorch
MAINTAINER Diego Borges <dborgesrmit@gmail.com>
USER root

ENV APT_PACKAGES=" \
		git \
		software-properties-common \
		apt-transport-https \
		ca-certificates \
		apt-utils \
#		libssl-dev \
#		libzmq3-dev \
#		libhdf5-7 \
#		hdf5-tools \
#		libhdf5-dev \
		libreadline-dev \
		zlib1g-dev \
		libgsl0-dev \
		wget \
#		vim \
#		mercurial \
#		wkhtmltopdf \
#		libtbb-dev \
#		libboost-dev \
#		cmake \
#		libmagickwand-dev \
		build-essential \
		curl \
		zip \
		"

ENV DEBIAN_FRONTEND=noninteractive
#Add ROCm 3.0 GPG key and fix broken apt-get
RUN wget -q -O - http://repo.radeon.com/rocm/apt/3.0/rocm.gpg.key | sudo apt-key add - && \
	add-apt-repository 'deb [arch=amd64] http://repo.radeon.com/rocm/apt/3.0/ xenial main' && \
	apt-get install -y apt-utils apt-transport-https

#Get apt packages
RUN apt-get update && \
	apt-get install -y --no-install-recommends $APT_PACKAGES && \
 	rm -rf "/var/lib/apt/lists/*" && \
	apt-get clean && \
	rm -rf /var/cache/apt

#Get conda dependencies
COPY environment.yml /tmp/environment.yml
#Get Anaconda and install needed packages
ENV CONDA_VERSION="4.7.12.1"
RUN curl -s -L https://repo.anaconda.com/miniconda/Miniconda3-$CONDA_VERSION-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -u -b -p /opt/conda && \
    rm miniconda.sh && \
    export PATH=/opt/conda/bin:$PATH && \
    conda env create -f /tmp/environment.yml && \
#    conda update conda && \
#    conda update --all --yes && \
	conda clean -tipy
ENV PATH=/opt/conda/envs/env/bin/:/opt/conda/bin:$PATH

#Move over built pytorch into conda envs and install selene-sdk using pip from installed environment
RUN cp -r /root/.local/lib/python3.6/site-packages /opt/conda/envs/env/lib/python3.6/ && \
	/opt/conda/envs/env/bin/pip install --no-cache-dir selene-sdk && \
	rm -rf /root/.local/lib/python3.6

# Set Euplotid as working directory
USER root
RUN mkdir /data && \
 		mkdir /data/Euplotid && \
		chmod -R 0777 /data/Euplotid
WORKDIR /data/Euplotid
COPY ./applotid.py ./start.sh /app/
RUN chmod 755 /app/applotid.py /app/start.sh

#set bash as default shell and add activation of environment to bash start
ENV SHELL=/bin/bash
#Set bash up to start and activate env environment
RUN echo "export PATH=/opt/conda/envs/env/bin/:/opt/conda/bin/:$PATH" > ~/.bashrc && \
	echo "source activate base" >> ~/.bashrc

#Configure jupyter
ENV JUPYTER_PREFIX="/root/"
RUN jupyter notebook --generate-config && \
	sed -i "/c.NotebookApp.open_browser/c c.NotebookApp.open_browser = False" "$JUPYTER_PREFIX".jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.ip/c c.NotebookApp.ip = '*'" "$JUPYTER_PREFIX".jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.notebook_dir/c c.NotebookApp.notebook_dir = '/data/Euplotid'" "$JUPYTER_PREFIX".jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" "$JUPYTER_PREFIX".jupyter/jupyter_notebook_config.py \
        && sed -i "/c.NotebookApp.iopub_data_rate_limit/c c.NotebookApp.iopub_data_rate_limit = 10000000" "$JUPYTER_PREFIX".jupyter/jupyter_notebook_config.py

#Define default startup behavior
ENV EUPLO_IMAGE=euplotid
ENV EUPLO_DEPLOY=false
ENV EUPLO_ARCH=x86
ENV EUPLO_OS=debian
ENV JUPYTER_PORT=8894

EXPOSE $JUPYTER_PORT
EXPOSE 80
CMD ["bash", "/app/start.sh"]
