# Euplotid image to take in FPKM, DNA-DNA interactions, ChiP-Seq, Chromatin Accessibility and calls and annotate/visualize INs
FROM debian:jessie-slim
MAINTAINER Diego Borges <dborgesrmit@gmail.com>
USER root

ENV APT_PACKAGES=" \
		git \
		software-properties-common \
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
RUN	apt-get update && \
	apt-get install -y --no-install-recommends $APT_PACKAGES && \
 	rm -rf "/var/lib/apt/lists/*" && \
	apt-get clean && \
	rm -rf /var/cache/apt
	
#Get conda dependencies 
COPY environment.yml /tmp/environment.yml
#Get Anaconda and install needed packages
ENV CONDA_VERSION="4.7.12.1"
RUN curl -s -L https://repo.anaconda.com/miniconda/Miniconda3-$CONDA_VERSION-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh && \
    export PATH=/opt/conda/bin:$PATH && \
	conda env create -f /tmp/environment.yml && \
#    conda update conda && \
#    conda update --all --yes && \
	conda clean -tipy 
ENV PATH=/opt/conda/envs/env/bin/:$PATH	

# Run Torch7 installation scripts
#RUN apt-get install sudo
#RUN cd /root && git clone https://github.com/torch/distro.git /root/torch --recursive && cd /root/torch && \
#  bash ./install-deps && \
#  ./install.sh
  
# Export environment variables manually
ARG LUA_PATH='/root/.luarocks/share/lua/5.1/?.lua;/root/.luarocks/share/lua/5.1/?/init.lua;/root/torch/install/share/lua/5.1/?.lua;/root/torch/install/share/lua/5.1/?/init.lua;./?.lua;/root/torch/install/share/luajit-2.1.0-beta1/?.lua;/usr/local/share/lua/5.1/?.lua;/usr/local/share/lua/5.1/?/init.lua'
ENV LUA_PATH $LUA_PATH
ENV LUA_CPATH "/root/.luarocks/lib/lua/5.1/?.so;/root/torch/install/lib/lua/5.1/?.so;./?.so;/usr/local/lib/lua/5.1/?.so;/usr/local/lib/lua/5.1/loadall.so"
ENV PATH /root/torch/install/bin:$PATH
ENV LD_LIBRARY_PATH /root/torch/install/lib:$LD_LIBRARY_PATH
ENV DYLD_LIBRARY_PATH /root/torch/install/lib:$DYLD_LIBRARY_PATH
ENV LUA_CPATH /root/torch/install/lib/?.so:$LUA_CPATH
ENV JOBLIB_TEMP_FOLDER /output_data


# Install lua packages
#RUN luarocks install luafilesystem
#RUN luarocks install dpnn
#RUN luarocks install dp
#RUN luarocks install torch
#RUN luarocks install nn
#RUN luarocks install inn

# clone Basset
RUN cd /root && git clone https://github.com/davek44/Basset.git /root/Basset --recursive
# Run torch-hdf5 installation scripts
#RUN cd /root/Basset/src && git clone https://github.com/davek44/torch-hdf5.git && cd torch-hdf5 && \
#  luarocks make 

#Get Deepbind
#RUN cd /root && git clone https://github.com/kundajelab/deepbind.git 

#Get Euplotid
#RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git

# Set Euplotid as working directory
WORKDIR /root/Euplotid

#Define default startup behavior
ENV EUPLO_IMAGE=euplotid
ENV EUPLO_DEPLOY=false
ENV EUPLO_ARCH=x86
ENV EUPLO_OS=debian
ENV JUPYTER_PORT=8894

COPY ./applotid.py ./start.sh /app/
RUN chmod 755 /app/applotid.py /app/start.sh
EXPOSE $JUPYTER_PORT
EXPOSE 80
CMD ["bash", "/app/start.sh"]

