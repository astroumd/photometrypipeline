from photopipe.reduction.auto.autoproc import autoproc

autoproc(
    datadir='/mnt/data/science/',
    imdir='/mnt/data/reduced/',
    nomastersky=True, redo=1
)
