from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import copyfile, copytree


bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir='/mnt/data/GW_TOO_18Aug2019/biases/',
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
    workdir='/mnt/data/GW_TOO_18Aug2019/dome_flats/',
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=False
)

science_dict = preproc.choose_science(
    'lmi',
    workdir='/mnt/data/GW_TOO_18Aug2019/DG19yjnjc/',
    targetdir='/work/data/GW_TOO_18Aug2019/DG19yjnjc/selected/',
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=False
)

preproc.mkmaster('lmi', bias_calib, 'bias')
preproc.mkmaster('lmi', flat_calib, 'flat')

copyfile('flat_SDSS-Z.fits', '/work/data/GW_TOO_18Aug2019/DG19yjnjc/selected/flat_SDSS-Z.fits')
copyfile('bias_C0.fits', '/work/data/GW_TOO_18Aug2019/DG19yjnjc/selected/bias_C0.fits')

autoproc(datadir='/work/data/GW_TOO_18Aug2019/DG19yjnjc/selected/',
         imdir='/work/data/GW_TOO_18Aug2019/DG19yjnjc/reduced/',
         redo=1, nomastersky=True)

copytree('/work/july/GRB161004A/r/reduced/', '/mnt/data/output/july/GRB161004A/r/july/GRB161004A/r/reduced/')

