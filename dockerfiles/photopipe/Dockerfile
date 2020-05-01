FROM joedurbak/astromatic

WORKDIR /work

RUN apt-get install -y vim nano emacs gedit

RUN git clone  https://github.com/astroumd/photometrypipeline.git \
	&& cd photometrypipeline \
	&& python setup.py build \
	&& rm /usr/local/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc \
	&& cp /work/photometrypipeline/matplotlibrc /usr/local/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda develop /work/photometrypipeline

RUN cd photometrypipeline && git checkout docker_development

RUN ln -s /work/photometrypipeline/reduction.py /usr/local/bin/photo-reduce \
    && ln -s /work/photometrypipeline/test.py /usr/local/bin/photo-test

CMD ['bash', 'photo-reduce', 'photo-test']
