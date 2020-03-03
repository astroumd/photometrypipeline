from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import copyfile, copytree
import os
import glob
from zipfile import ZipFile

base_path = os.path.dirname(os.path.abspath(__file__))
test_path = os.path.join(base_path, 'test')
copy_path = os.path.join(base_path, 'copy')
selected_path = os.path.join(copy_path, 'selected')
reduced_path = os.path.join(copy_path, 'reduced')
zip_file = os.path.join('test_path', 'test.zip')
zf = ZipFile(zip_file, 'r')
zf.extractall(copy_path)

bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir=copy_path,
    cams=[0],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=True
)

flat_calib = preproc.choose_calib(
    'lmi',
    'flat',
    workdir=copy_path,
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=True
)

science_dict = preproc.choose_science(
    'lmi',
    workdir=copy_path,
    targetdir=selected_path,
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=True
)

preproc.mkmaster('lmi', bias_calib, 'bias')
preproc.mkmaster('lmi', flat_calib, 'flat')

for f in glob.glob(os.path.join(base_path, 'bias*.fits')):
    copyfile(os.path.join(base_path, f), os.path.join(selected_path, f))

for f in glob.glob(os.path.join(base_path, 'flat*.fits')):
    copyfile(os.path.join(base_path, f), os.path.join(selected_path, f))

autoproc(datadir=selected_path,
         imdir=reduced_path,
         redo=1, nomastersky=True)
