
FROM  maxjiang/rpi-jupyter:datascience


RUN apt-get update && apt-get install \
	libjpeg-dev

RUN pip3 install  \
	plotly \
	sense-hat

RUN apt-get update && apt-get install -y git
#Get Louvain
RUN cd /root && git clone https://github.com/taynaud/python-louvain.git && cd python-louvain && python3 setup.py install

#Get Euplotid  
RUN cd /root && git clone https://github.com/dborgesr/Euplotid.git
# Set pipelines as workdir
WORKDIR /root/Euplotid

#Configure jupyter
RUN sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py


EXPOSE 8893
CMD ["jupyter", "notebook", "--port=8893", "--no-browser", "--ip=0.0.0.0"]

