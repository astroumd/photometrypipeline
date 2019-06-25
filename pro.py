from photopipe.reduction.auto.autoproc import autoproc

autoproc(
    datadir='/mnt/data/science/selected/',
    imdir='/mnt/data/science/selected/reduced/',
    nomastersky=True, redo=1
)
