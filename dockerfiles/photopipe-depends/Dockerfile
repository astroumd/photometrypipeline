FROM stelais/astromatic

WORKDIR /work

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda install -y -c conda-forge astroquery \
    && pip install astride astroscrappy
