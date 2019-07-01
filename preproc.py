from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import copyfile


bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir='/mnt/data/i/bias/',
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
    workdir='/mnt/data/i/flat/',
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=False
)

science_dict = preproc.choose_science(
    'lmi',
    workdir='/mnt/data/i/science/',
    targetdir='/mnt/data/i/science_selected/',
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=False
)

preproc.mkmaster('lmi', bias_calib, 'bias')
preproc.mkmaster('lmi', flat_calib, 'flat')

copyfile('flat_SDSS-I.fits', '/mnt/data/science_selected/flat_SDSS-I.fits')
copyfile('bias_C0.fits', '/mnt/data/science_selected/bias_C0.fits')

autoproc(datadir='/mnt/data/i/science_selected/',
         imdir='/mnt/data/i/reduced/',
         redo=1, nomastersky=True)
