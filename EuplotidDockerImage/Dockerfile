# Euplotid Anaconda2 image to take in FPKM, DNA-DNA interactions, ChiP-Seq, Chromatin Accessibility and calls and annotate/visualize INs
FROM debian:jessie
MAINTAINER Diego Borges <dborgesr@mit.edu>
USER root

# Install git, apt-add-repository and dependencies for iTorch
RUN apt-get update && apt-get install -y \
  git \
  software-properties-common \
  libssl-dev \
  libzmq3-dev \
  python-zmq \
  python-pip \
  bedtools \
  samtools \
#  libhdf5-7 \
  hdf5-tools \
  libhdf5-dev \
  wget \
  vim \
  mercurial \
  wkhtmltopdf \
  libtbb-dev \
  libboost-dev \
  cmake \
  libmagickwand-dev \
  && apt-get clean

#Get Anaconda
RUN cd /root && wget https://repo.continuum.io/archive/Anaconda2-4.2.0-Linux-x86_64.sh && bash Anaconda2-4.2.0-Linux-x86_64.sh -b 
ENV PATH=/root/anaconda2/bin:$PATH

# Install Jupyter Notebook for iTorch and stray libraries
RUN pip install --upgrade pip
RUN pip install notebook ipywidgets seaborn pysam myvariant mygene weblogo h5py wand mechanize rst2pdf argparse python-louvain uuid plotly joblib biopython pyliftover nbbrowserpdf bx-python MACS2 cython

#Install bash kernel
RUN pip install bash_kernel
RUN python -m bash_kernel.install

# Run Torch7 installation scripts
RUN apt-get install sudo
RUN cd /root && git clone https://github.com/torch/distro.git /root/torch --recursive && cd /root/torch && \
  bash ./install-deps && \
  ./install.sh
  
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
RUN luarocks install luafilesystem
RUN luarocks install dpnn
RUN luarocks install dp
RUN luarocks install torch
RUN luarocks install nn
#RUN luarocks install inn

# clone Basset
RUN cd /root && git clone https://github.com/davek44/Basset.git /root/Basset --recursive
# Run torch-hdf5 installation scripts
RUN cd /root/Basset/src && git clone https://github.com/davek44/torch-hdf5.git && cd torch-hdf5 && \
  luarocks make 
  
#Get bedops
RUN cd /root && git clone https://github.com/bedops/bedops.git && cd ./bedops && make && make install
ENV PATH="/root/bedops/bin:$PATH"

#Get Deepbind
RUN cd /root && git clone https://github.com/kundajelab/deepbind.git 

#Get Euplotid
RUN cd /root && git clone -b master --single-branch https://github.com/dborgesr/Euplotid.git

# Set Euplotid as working directory
WORKDIR /root/Euplotid

RUN apt-get clean \
 && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Add Tini. Tini operates as a process subreaper for jupyter. This prevents
# kernel crashes.
ENV TINI_VERSION v0.6.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT ["/usr/bin/tini", "--"]

EXPOSE 8890
CMD ["jupyter", "notebook", "--port=8890", "--no-browser", "--ip=0.0.0.0"]
