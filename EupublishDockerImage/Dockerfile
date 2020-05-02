# Eupublish image to write and publish ipynb into pdf
FROM debian:stretch-slim
MAINTAINER Diego Borges <dborgesrmit@gmail.com>
USER root

ENV APT_PACKAGES=" \
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
  build-essential \
  "

ENV DEBIAN_FRONTEND=noninteractive
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

#get jekyll things for WEBSITE_PORT
RUN gem install jekyll \
                jekyll-include-cache \
                jekyll-feed \
                jekyll-gist \
                jekyll-sitemap \
                jekyll-paginate

# Set Euplotid as working directory
USER root
RUN mkdir /data && \
 		mkdir /data/Euplotid && \
		chmod -R 0777 /data/Euplotid
WORKDIR /data/Euplotid
COPY ./start.sh /app/
RUN chmod 755 /app/start.sh

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
ENV JUPYTER_PORT=8895
ENV WEBSITE_PORT=4000

EXPOSE $JUPYTER_PORT
EXPOSE $WEBSITE_PORT
EXPOSE 80
CMD ["bash", "/app/start.sh"]
