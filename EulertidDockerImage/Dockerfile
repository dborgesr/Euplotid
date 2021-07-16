# Eulertid Anaconda3 image to get started w/ programming
FROM debian:stretch-slim
USER root

#Set up shell for install
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash","-c"]
ENV BASH_EN=~/.bashrc
ENV SHELL=/bin/bash

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
  git \
  && apt-get update \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 
	
#Get micromamba
COPY environment.yml /tmp/environment.yml
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
RUN  wget -qO- https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj bin/micromamba --strip-components=1 && \
	chmod 755 ./micromamba && \
	mkdir -p $(dirname $MAMBA_ROOT_PREFIX) && \
	./micromamba shell init -s bash -p $MAMBA_ROOT_PREFIX && \
	echo "micromamba activate base" >> /root/.bashrc && \
	source ~/.bashrc && \
	./micromamba install -y -n base -f /tmp/environment.yml && \
	rm ${MAMBA_ROOT_PREFIX}/lib/*.a && \
	rm -rf ${MAMBA_ROOT_PREFIX}/pkgs

#Get jupyter notebooks from yamaton for Euler problems
RUN cd /root && git clone https://github.com/yamaton/project-euler-jupyter.git

#Enable C kernel after installed using pip above (jupyter-c-kernel)
#RUN install_c_kernel	

#Configure jupyter
ENV JUPYTER_PREFIX="/root/"
RUN jupyter lab --generate-config && \
	sed -i "/c.ServerApp.open_browser/c c.ServerApp.open_browser = False" "$JUPYTER_PREFIX".jupyter/jupyter_lab_config.py \
        && sed -i "/c.ServerApp.ip/c c.ServerApp.ip = '*'" "$JUPYTER_PREFIX".jupyter/jupyter_lab_config.py \
        && sed -i "/c.ServerApp.root_dir/c c.ServerApp.root_dir = '$WORK_DIR'" "$JUPYTER_PREFIX".jupyter/jupyter_lab_config.py \
        && sed -i "/c.ServerApp.token/c c.ServerApp.token = 'eulertid'" "$JUPYTER_PREFIX".jupyter/jupyter_lab_config.py \
        && sed -i "/c.ServerApp.iopub_data_rate_limit/c c.ServerApp.iopub_data_rate_limit = 10000000" "$JUPYTER_PREFIX".jupyter/jupyter_lab_config.py

# Set Euler problems as working directory
WORKDIR /root/project-euler-jupyter

EXPOSE 8896
CMD ["jupyter", "lab", "--port=8896", "--no-browser", "--allow-root", "--ip=0.0.0.0"]

