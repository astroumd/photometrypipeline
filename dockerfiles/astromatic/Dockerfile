FROM continuumio/anaconda3

WORKDIR /work
VOLUME /myvolume

RUN conda install numpy scipy matplotlib astropy jsonschema ipython
# && conda install -c cefca pyfits

RUN apt-get update -y \
	&& apt-get install -y apt-utils \
	git \
	sextractor \
	swarp \
	scamp \
	missfits \
	gfortran \
	libatlas-base-dev \
	libpng-dev \
	libfreetype6-dev \
	libxft-dev \
	python-tk \
	libsm-dev \
	x11-apps \
	saods9 \
	make \
	wget

RUN cp /usr/bin/sextractor /usr/bin/sex && cp /usr/bin/SWarp /usr/bin/swarp

RUN mkdir dependencies && cd dependencies

RUN cd dependencies && wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz \
    && tar -zxvf cdsclient.tar.gz \
    && cd cdsclient-3.84/ \
    && ./configure \
    && make \
    && make install

RUN cd dependencies && wget http://ds9.si.edu/download/debian9/ds9.debian9.8.1.tar.gz  \
    && tar xzvf ds9.debian9.8.1.tar.gz \
	&& mv ds9 /usr/local/bin \
	&& rm /usr/bin/ds9
