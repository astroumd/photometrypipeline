import argparse
import os
import glob
from shutil import move

from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc

parser = argparse.ArgumentParser()
parser.add_argument('data_directory', help='directory that contains your data, e.g. /mnt/data')
parser.add_argument('--noplot', help='turns off plots during the frame selection process', action='store_true')
args = parser.parse_args()

noplot = args.noplot

print("creating paths")
# base_path = os.path.abspath(os.path.dirname(__file__))  # /path/to/data
# print('base_path', base_path)
# test_path = os.path.join(base_path, 'test')  # /path/to/data/test
# print('test_path', test_path)
input_data_path = args.data_directory  # /path/to/data/
print('input_data_path', input_data_path)
selected_path = os.path.join(input_data_path, 'selected')  # /path/to/data/test/extracted_data/selected
print('selected_path', selected_path)
reduced_path = os.path.join(input_data_path, 'reduced')  # /path/to/data/test/extracted_data/reduced
print('reduced_path', reduced_path)

print('start bias calibration selection')
bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir=input_data_path + os.path.sep,
    cams=[0],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=noplot
)

print('start flat calibration selection')
flat_calib = preproc.choose_calib(
    'lmi',
    'flat',
    workdir=input_data_path + os.path.sep,
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=True,
    save_select=True,
    noplot=noplot
)

print('start science frame selection')
science_dict = preproc.choose_science(
    'lmi',
    workdir=input_data_path+os.path.sep,
    targetdir=selected_path+os.path.sep,
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=noplot
)

print('start mkmaster bias')
preproc.mkmaster('lmi', bias_calib, 'bias')
print('start mkmaster flat')
preproc.mkmaster('lmi', flat_calib, 'flat')

print('start move master bias files to "selected" folder')
for f in glob.glob('bias*.fits'):
    filename = os.path.basename(f)  # e
    print('selected_path', selected_path)
    f_selected_path = os.path.join(selected_path, filename)
    print('f_selected_path', f_selected_path)
    print('moving {0} to {1}'.format(filename, f_selected_path))
    move(filename, f_selected_path)

print('start move files master flats to selected folder')
for f in glob.glob('flat*.fits'):
    filename = os.path.basename(f)
    print('selected_path', selected_path)
    f_selected_path = os.path.join(selected_path, filename)
    print('f_selected_path', f_selected_path)
    print('moving {0} to {1}'.format(filename, f_selected_path))
    move(filename, f_selected_path)

print('start reduction')
autoproc(datadir=selected_path+os.path.sep,
         imdir=reduced_path+os.path.sep,
         redo=1, nomastersky=True)
