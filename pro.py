from photopipe.reduction.auto.autoproc import autoproc

autoproc(
    datadir='/mnt/data/selected/',
    imdir='/mnt/data/reduced/',
    nomastersky=True, redo=1, step='astrometry'
)
