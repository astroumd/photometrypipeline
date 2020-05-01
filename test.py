from photopipe.reduction import preproc
from photopipe.reduction.auto.autoproc import autoproc
from shutil import move
import os
import glob
from zipfile import ZipFile

# testing git commit and push
print("creating paths")
base_path = os.path.abspath(os.path.dirname(__file__))  # /path/to/photometrypipeline
print('base_path', base_path)
test_path = os.path.join(base_path, 'test')  # /path/to/photometrypipeline/test
print('test_path', test_path)
input_data_path = os.path.join(test_path, 'extracted_data')  # /path/to/photometrypipeline/test/extracted_data
print('copy_path', input_data_path)
selected_path = os.path.join(input_data_path, 'selected')  # /path/to/photometrypipeline/test/extracted_data/selected
print('selected_path', selected_path)
reduced_path = os.path.join(input_data_path, 'reduced')  # /path/to/photometrypipeline/test/extracted_data/reduced
print('reduced_path', reduced_path)
zip_file = os.path.join(test_path, 'test.zip')  # /path/to/photometrypipeline/test/test.zip
print(zip_file, zip_file)
zf = ZipFile(zip_file, 'r')
print('extracting files')
zf.extractall(input_data_path)
print('extraction complete')

print('start bias calibration')
bias_calib = preproc.choose_calib(
    'lmi',
    'bias',
    workdir=input_data_path + os.path.sep,
    cams=[0],
    auto=True,
    amin=0.0, amax=1.0,
    reject_sat=False,
    save_select=True,
    noplot=True
)

print('start flat calibration')
flat_calib = preproc.choose_calib(
    'lmi',
    'flat',
    workdir=input_data_path + os.path.sep,
    cams=[0],
    auto=True,
    amin=0.2, amax=0.8,
    reject_sat=False,
    save_select=True,
    noplot=True
)

print('start choose science')
science_dict = preproc.choose_science(
    'lmi',
    workdir=input_data_path+os.path.sep,
    targetdir=selected_path+os.path.sep,
    cams=[0],
    auto=True,
    save_select=True,
    calibrate=False,
    noplot=True
)

print('start mkmaster bias')
preproc.mkmaster('lmi', bias_calib, 'bias')
print('start mkmaster flat')
preproc.mkmaster('lmi', flat_calib, 'flat')

print('start move files master biases to selected folder')
for f in glob.glob('bias*.fits'):
    filename = os.path.basename(f)  # e
    print('selected_path', selected_path)
    # f_base_path = os.path.join(base_path, filename)
    # print('f_base_path', f_base_path)
    f_selected_path = os.path.join(selected_path, filename)
    print('f_selected_path', f_selected_path)
    print('moving {0} to {1}'.format(filename, f_selected_path))
    move(filename, f_selected_path)

print('start move files master flats to selected folder')
for f in glob.glob('flat*.fits'):
    filename = os.path.basename(f)
    print('selected_path', selected_path)
    # f_base_path = os.path.join(base_path, filename)
    # print('f_base_path', f_base_path)
    f_selected_path = os.path.join(selected_path, filename)
    print('f_selected_path', f_selected_path)
    print('moving {0} to {1}'.format(filename, f_selected_path))
    move(filename, f_selected_path)

autoproc(datadir=selected_path+os.path.sep,
         imdir=reduced_path+os.path.sep,
         redo=1, nomastersky=True)
