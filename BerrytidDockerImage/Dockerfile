# This file creates a container that runs a jupyter notebook server on Raspberry Pi

#FROM debian/jessie-slim
FROM easypi/alpine-arm

ENV ALPINE_VERSION=3.6

ENV PACKAGES="\
  dumb-init \
  musl \
  linux-headers \
  build-base \
  bash \
  git \
  ca-certificates \
  python3.4 \
  python3.4-dev \
  wget \
  libxml2-dev \
  libxslt-dev \
  libffi-dev \
  gcc \
  musl-dev \
  libgcc \
  openssl-dev \
  curl \
  jpeg-dev \
  zlib-dev \
  freetype-dev \
  lcms2-dev \
  openjpeg-dev \
  tiff-dev \
  tk-dev \
  tcl-dev

RUN echo \
  # replacing default repositories with edge ones
  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/testing" > /etc/apk/repositories \
  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories \
  && echo "http://dl-cdn.alpinelinux.org/alpine/edge/main" >> /etc/apk/repositories \

  # Add the packages, with a CDN-breakage fallback if needed
  && apk add --no-cache $PACKAGES || \
    (sed -i -e 's/dl-cdn/dl-4/g' /etc/apk/repositories && apk add --no-cache $PACKAGES) \

  # turn back the clock -- so hacky!
  && echo "http://dl-cdn.alpinelinux.org/alpine/v$ALPINE_VERSION/main/" > /etc/apk/repositories \
  # && echo "@edge-testing http://dl-cdn.alpinelinux.org/alpine/edge/testing" >> /etc/apk/repositories \
  # && echo "@edge-community http://dl-cdn.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories \
  # && echo "@edge-main http://dl-cdn.alpinelinux.org/alpine/edge/main" >> /etc/apk/repositories \

  # make some useful symlinks that are expected to exist
  && if [[ ! -e /usr/bin/python ]];        then ln -sf /usr/bin/python3.4 /usr/bin/python; fi \
  && if [[ ! -e /usr/bin/python-config ]]; then ln -sf /usr/bin/python3.4-config /usr/bin/python-config; fi \
  && if [[ ! -e /usr/bin/idle ]];          then ln -sf /usr/bin/idle3.4 /usr/bin/idle; fi \
  && if [[ ! -e /usr/bin/pydoc ]];         then ln -sf /usr/bin/pydoc3.4 /usr/bin/pydoc; fi \
  && if [[ ! -e /usr/bin/easy_install ]];  then ln -sf /usr/bin/easy_install-3.4 /usr/bin/easy_install; fi \

  # Install and upgrade Pip
  && easy_install pip \
  && pip install --upgrade pip \
  && if [[ ! -e /usr/bin/pip ]]; then ln -sf /usr/bin/pip3.4 /usr/bin/pip; fi

RUN ln -s /usr/include/locale.h /usr/include/xlocale.h
RUN pip install  \
	jupyter \
	plotly \
	sense-hat

#Get Louvain
RUN cd /root && git clone https://github.com/taynaud/python-louvain.git && cd python-louvain && python setup.py install

#Get RTIMULIB, not working yet :(
RUN cd /root && git clone https://github.com/RPi-Distro/RTIMULib.git

#Get Euplotid  
RUN cd /root && git clone https://github.com/dborgesr/Euplotid.git
# Set pipelines as workdir
WORKDIR /root/Euplotid

#Configure jupyter
RUN jupyter notebook --generate-config --allow-root && \
	sed -i "/c.NotebookApp.token/c c.NotebookApp.token = 'euplotid'" /root/.jupyter/jupyter_notebook_config.py


EXPOSE 8893
CMD ["jupyter", "notebook", "--allow-root", "--port=8893", "--no-browser", "--ip=0.0.0.0"]

