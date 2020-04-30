FROM joedurbak/astromatic

WORKDIR /work

RUN apt-get install -y vim nano emacs gedit

RUN git clone  https://github.com/astroumd/photometrypipeline.git \
	&& cd photometrypipeline \
	&& python setup.py build \
	&& rm /usr/local/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc \
	&& cp /work/photometrypipeline/matplotlibrc /usr/local/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc

RUN conda install conda-build \
    && conda develop /work/photometrypipeline
	