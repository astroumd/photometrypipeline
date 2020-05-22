FROM stelais/astromatic

WORKDIR /work

RUN apt-get install -y vim nano emacs gedit

RUN git clone --branch py3 https://c2499497b5dd9fc93f9d67000bc15f1c514b850f@github.com/RIMAS-RATIR-DCT/photometrypipeline.git \
	# && git checkout py3 \
	&& cd photometrypipeline \
	&& python setup.py build \
	&& rm /usr/local/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc \
	&& cp /work/photometrypipeline/matplotlibrc /usr/local/lib/python3.7/site-packages/matplotlib/mpl-data/matplotlibrc

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda develop /work/photometrypipeline
