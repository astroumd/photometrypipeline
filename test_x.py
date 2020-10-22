from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import move
import os
import glob


# testing git commit and push
print("creating paths")
base_path = os.path.dirname(os.path.abspath(__file__))
print('base_path', base_path)
test_path = os.path.join(base_path, 'test')
print('test_path', test_path)
copy_path = os.path.join(test_path, 'IR')
print('copy_path', copy_path)
selected_path = os.path.join(copy_path, 'selected')
print('selected_path', selected_path)
reduced_path = os.path.join(copy_path, 'reduced')
print('reduced_path', reduced_path)

bias_calib = preproc.choose_calib(
    'ratir',
    'bias',
    workdir=copy_path+os.path.sep,
    cams=[0,1,2,3],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=True
)

print('start dark calibration')
dark_calib = preproc.choose_calib(
    'ratir',
    'dark',
    workdir=copy_path+os.path.sep,
    cams=[0,1,2,3],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=True
)

print('start mkmaster bias')
preproc.mkmaster('ratir', bias_calib, 'bias')
print('start mkmaster dark')
preproc.mkmaster('ratir', dark_calib, 'dark')

