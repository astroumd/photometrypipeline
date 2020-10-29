FROM stelais/astromatic

WORKDIR /work

RUN apt-get install -y vim nano

#RUN git clone --branch py3 https://c2499497b5dd9fc93f9d67000bc15f1c514b850f@github.com/RIMAS-RATIR-DCT/photometrypipeline.git \
RUN git clone --branch py3_astrometry https://github.com/astroumd/photometrypipeline.git \
	&& cd photometrypipeline \
	&& python setup.py build

#RUN export LANG=C.UTF-8 \
#    && conda update conda \
#    && conda install conda-build \
#    && conda develop /work/photometrypipeline

