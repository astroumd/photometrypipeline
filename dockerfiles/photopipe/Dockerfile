FROM stelais/astromatic

WORKDIR /work

RUN pip install astride
RUN conda install -c conda-forge astroquery

#RUN apt-get install -y vim nano emacs gedit

RUN git clone --branch rimas_py3 https://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@github.com/RIMAS-RATIR-DCT/photometrypipeline.git \

#RUN git clone https://github.com/joedurbak/photometrypipeline.git \
    && cd photometrypipeline \
    && python setup.py build \
    && chmod a+x photopipe_cli.py

RUN export LANG=C.UTF-8 \
    && conda update conda \
    && conda install conda-build \
    && conda develop /work/photometrypipeline

CMD ["git -C /work/photometrypipeline pull"]
ENTRYPOINT ["/work/photometrypipeline/photopipe_cli.py"]
CMD ["-h"]
