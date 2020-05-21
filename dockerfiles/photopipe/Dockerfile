FROM stelais/astromatic

WORKDIR /work

RUN apt-get install -y vim nano emacs gedit

RUN git clone https://github.com/RIMAS-RATIR-DCT/photometrypipeline.git \
	&& cd photometrypipeline \
	&& python setup.py build \
	&& rm /usr/local/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc \
	&& cp /work/photometrypipeline/matplotlibrc /usr/local/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda develop /work/photometrypipeline
	