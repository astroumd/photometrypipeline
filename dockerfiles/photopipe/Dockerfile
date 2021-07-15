FROM stelais/astromatic

WORKDIR /work

RUN pip install astride

#RUN apt-get install -y vim nano emacs gedit

#RUN git clone --branch py3 https://c2499497b5dd9fc93f9d67000bc15f1c514b850f@github.com/RIMAS-RATIR-DCT/photometrypipeline.git \

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
