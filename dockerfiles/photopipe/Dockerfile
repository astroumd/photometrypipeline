FROM joedurbak/photopipe-depends

WORKDIR /work

RUN git clone --branch rimas_py3 https://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@github.com/RIMAS-RATIR-DCT/photometrypipeline.git \
    && cd photometrypipeline \
    && python setup.py build \
    && chmod a+x photopipe_cli.py

RUN conda develop /work/photometrypipeline

CMD ["git -C /work/photometrypipeline pull"]
ENTRYPOINT ["/work/photometrypipeline/photopipe_cli.py"]
CMD ["-h"]
