FROM stelais/astromatic

WORKDIR /work

RUN pip install astride
RUN conda install -c conda-forge astroquery

RUN git clone https://github.com/joedurbak/photometrypipeline.git \
    && cd photometrypipeline \
    && python setup.py build \
    && chmod a+x photopipe_cli.py

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda develop /work/photometrypipeline

ENTRYPOINT ["/work/photometrypipeline/photopipe_cli.py"]
CMD ["-h"]
