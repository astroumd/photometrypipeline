from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import copyfile, copytree


bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir='/mnt/data/july/bias/',
    cams=[0],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=False
)

flat_calib = preproc.choose_calib(
    'lmi',
    'flat',
    workdir='/mnt/data/july/dome_flats/',
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=False
)

science_dict = preproc.choose_science(
    'lmi',
    workdir='/mnt/data/july/GRB161004A/r/',
    targetdir='/mnt/data/july/GRB161004A/r/selected/',
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=False
)

preproc.mkmaster('lmi', bias_calib, 'bias')
preproc.mkmaster('lmi', flat_calib, 'flat')

copyfile('flat_SDSS-R.fits', '/mnt/data/july/GRB161004A/r/selected/flat_SDSS-R.fits')
copyfile('bias_C0.fits', '/mnt/data/july/GRB161004A/r/selected/bias_C0.fits')

autoproc(datadir='/mnt/data/july/GRB161004A/r/selected/',
         imdir='/work/july/GRB161004A/r/reduced/',
         redo=1, nomastersky=True)

copytree('/work/july/GRB161004A/r/reduced/', '/mnt/data/output/july/GRB161004A/r/july/GRB161004A/r/reduced/')

