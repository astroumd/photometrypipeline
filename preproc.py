from photopipe.reduction import preproc

bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir='/mnt/data/bias/',
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
    workdir='/mnt/data/flat/',
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=False
)

preproc.mkmaster('lmi', bias_calib, 'bias')
preproc.mkmaster('lmi', flat_calib, 'flat')
