# Minitid image to take in FPKM, normalized DNA-DNA interactions, ChiP-Seq, Chromatin Accessibility and build, annotate, and visualize INs
FROM frolvlad/alpine-glibc
MAINTAINER Diego Borges <dborgesr@mit.edu>
USER root

ENV CONDA_DIR="/opt/conda"
ENV PATH="$CONDA_DIR/bin:$PATH"

# Install conda
RUN CONDA_VERSION="4.0.5" && \
    CONDA_MD5_CHECKSUM="b1b15a3436bb7de1da3ccc6e08c7a5df" && \
    \
    apk add --no-cache --virtual=.build-dependencies wget ca-certificates bash && \
    \
    mkdir -p "$CONDA_DIR" && \
    wget "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh" -O miniconda.sh && \
    echo "$CONDA_MD5_CHECKSUM  miniconda.sh" | md5sum -c && \
    bash miniconda.sh -f -b -p "$CONDA_DIR" && \
    echo "export PATH=$CONDA_DIR/bin:\$PATH" > /etc/profile.d/conda.sh && \
    rm miniconda.sh && \
    \
    conda install -y -c soumith \
    	pytorch \
    	torchvision \
    	jupyter \
    	notebook \
    	nbpresent \
    	git && \
    conda install -y -c edurand \
    	myvariant \
    	mygene && \
    conda install -y -c bioconda \
    	pysam \
    	networkx \
    	pandas \
    	joblib \
    	biopython \
    	pyliftover \
    	scipy \
    	bedops \
    	plotly && \
    conda update --all --yes && \
    conda config --set auto_update_conda False && \
    conda clean --all --yes && \
    echo "Removing older versions of Python due to the bug #3177 in 'conda clean'" && \
    rm -rf $(find "$CONDA_DIR/pkgs/" -maxdepth 1 -name "python-3*" | sort | head -n -1) && \
    \
    apk del --purge .build-dependencies && \
    \
    mkdir -p "$CONDA_DIR/locks" && \
    chmod 777 "$CONDA_DIR/locks"

RUN cd /root && git clone https://github.com/taynaud/python-louvain.git && cd python-louvain && python setup.py install
RUN cd /root && git clone https://github.com/dborgesr/Euplotid.git
# Set pipelines as workdir
WORKDIR /root/Euplotid

# Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py    
ENV TINI_VERSION v0.14.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

EXPOSE 8892
CMD ["jupyter", "notebook", "--allow-root", "--port=8892", "--no-browser", "--ip=0.0.0.0"]