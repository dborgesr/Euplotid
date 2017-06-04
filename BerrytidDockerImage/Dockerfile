# This file creates a container that runs a jupyter notebook server on Raspberry Pi

FROM easypi/alpine-arm

RUN apk add --no-cache \
	bash \
	wget \
	git \
	alpine-sdk \
	python3 \
	python3-dev 

RUN pip3 install --upgrade pip  && \
	pip3 install \ 
		jupyter \
		plotly  
#Get Louvain
RUN cd /root && git clone https://github.com/taynaud/python-louvain.git && cd python-louvain && python3 setup.py install

#Get Euplotid  
RUN cd /root && git clone https://github.com/dborgesr/Euplotid.git
# Set pipelines as workdir
WORKDIR /root/Euplotid

#Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py    
ENV TINI_VERSION v0.14.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

EXPOSE 8893
CMD ["jupyter", "notebook", "--allow-root", "--port=8893", "--no-browser", "--ip=0.0.0.0"]

