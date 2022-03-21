#!/usr/bin/env python
# coding=utf-8
import argparse
import glob
import os
import shutil
from zipfile import ZipFile

from photopipe.reduction.preprocess import choose, master
from photopipe.reduction.auto.autoproc import autoproc
from photopipe.photometry.autoredux import autoredux


def list_to_options_str(options):
    _options = ["'{}'".format(i) for i in options]
    _str = ', '.join(_options)
    return ', options are: ' + _str


parser = argparse.ArgumentParser()

cmd_options = ['pipeline_all', 'preprocess', 'autoproc', 'photometry', 'ipython', 'bash', 'test']
cmd_str = list_to_options_str(cmd_options)
instrument_options = ['lmi', 'ratir', 'rimas']
instrument_str = list_to_options_str(instrument_options)
autoproc_step_options = ['prepare', 'flatten', 'makesky', 'skysub', 'crclean', 'astrometry', 'zpoint', 'stack']
autoproc_step_str = list_to_options_str(autoproc_step_options)

# overall arguments
parser.add_argument('cmd', type=str, help="'command to run"+cmd_str)
parser.add_argument(
    '-d', '--datadir', type=str, default='/data', help="""
raw data directory, default: /data

pipeline will create nested directories for the following steps, for example:

data  
│
|
│   20160628T032914C0b.fits
│   20160628T032914C1b.fits
│   ...
|
└───calibration
|
└───selected
│
└───reduced
│
└───coadd

*NOTE* If you are using docker, you need to have the data directory mounted as a volume, e.g. -v /path/to/data:/data
"""
)
parser.add_argument('--nopreprocess', action="store_true", default=False, help="""
skips all preprocessing steps: 
    - master calibration frame creation
        - flats, darks, biases
    - science data formatting
""")
parser.add_argument('--noautoproc', action="store_true", default=False, help="""
skips all of the reduction steps:
    - prepare, flatten, makesky, skysub, crclean, astrometry, zpoint, stack
""")
parser.add_argument('--nophotometry', action="store_true", default=False, help="""
skips the photometry step
""")

# preprocess arguments
parser.add_argument(
    '--instrument', type=str, default='lmi', help='instrument data being processed, default: "lmi"'+instrument_str
)
parser.add_argument('--cams', type=str, default=None, help="camera numbers, example: --cams '0,1,2,3', default: all")
parser.add_argument('--noautoselect', action="store_true", default=False, help="""
- no automated frame selection. If 'bias', will select all, if 'flat' will select
    non-saturated frames with sufficient counts
""")
parser.add_argument('--norejectsat', default=False, action="store_true",
                    help="don't reject frames with saturated pixels")
parser.add_argument(
    '--amin', type=float, default=0.2,
    help='minimum fraction of saturation value for median during automated flat selection, default: 0.2'
)
parser.add_argument(
    '--amax', type=float, default=0.8,
    help='maximum fraction of saturation value for median during automated flat selection, default: 0.8'
)
parser.add_argument(
    '--nosaveselect', default=False, action="store_true",
    help="don't save dictionary of selected frames to python pickle file during preprocessing"
)
parser.add_argument(
    '--fmin', default=5, type=int, help='Minimum number of files needed to make a master frame, default: 5'
)
parser.add_argument(
    '--noplot', default=False, action="store_true",
    help="don't show selected frames"
)
parser.add_argument(
    '--nobias', action="store_true", default=False,
    help="skip master bias creation, useful if you already have a master bias in the datadir or have no bias data"
)
parser.add_argument(
    '--noflat', action="store_true", default=False,
    help="skip master flat creation, useful if you already have a master flat in the datadir or have no flat data"
)
parser.add_argument(
    '--nodark', action="store_true", default=False,
    help="skip master dark creation, useful if you already have a master dark in the datadir or have no dark data"
)
parser.add_argument(
    '--noscienceselect', action="store_true", default=False,
    help="skip science data selection, useful if you already did this"
)

# autoproc arguments
parser.add_argument('--start', type=str, help="autoproc starting step"+autoproc_step_str, default=None)
parser.add_argument('--stop', type=str, help="autoproc stopping step"+autoproc_step_str, default=None)
parser.add_argument('--step', type=str, help="autoproc only do this step"+autoproc_step_str, default=None)
parser.add_argument('--nocrclean', action="store_true", help='skip cosmic ray removal', default=False)
parser.add_argument(
    '--nomastersky', action="store_true", help='skip master sky creation and use median value instead', default=False
)
parser.add_argument(
    '--noskyflattarg', action="store_true",
    help='skip creating master sky by target, make master sky using entire observation run instead', default=False
)
parser.add_argument('--redo', action="store_true", help='redo steps that have already been completed', default=False)
parser.add_argument('--quiet', action="store_true", help='reduces the number of print statements', default=False)
parser.add_argument('--rmifiles', action="store_true", help='removes intermediate files', default=False)
parser.add_argument(
    '--customcat', type=str, default=None,
    help='''
    Custom catalog (txt file) to determine instrumental zeropoint
    corrections
    Must be in same format as what get_SEDs.py produces
    (i.e. ra(deg)	dec(deg)	u	g	r	i	z	y	B	V	R
    I	J	H	K	u_err	g_err	r_err	i_err	y_err	B_err
    V_err	R_err	I_err	J_err	H_err	K_err	Mode)
    First line will be skipped (use for headers)
    everything but JHK are expected to be in AB magnitudes, JHK should
    be in same units as 2MASS (Vega mags)
    
    File is expected to be in datadir location
    '''
)
parser.add_argument(
    '--customcatfilt', type=str, default=None, help='''
    Filters relevant to custom catalog file (all other filters will
    use get_SEDs.py to calculate catalog from 2MASS + (SDSS or APASS
    or USNOB1) in that order
    
    Usage: --customcatfilt "u,g,z"
    '''

)
parser.add_argument(
    '--nogaia', action="store_true", default=False,
    help="Do not use Gaia catalog for astrometric calibration"
)
parser.add_argument('--debug', help='turn on debug outputs', default=False, action='store_true')
args = parser.parse_args()
data_dir = args.datadir + '/'
cal_dir = os.path.join(args.datadir, 'calibration') + os.path.sep
selected_dir = os.path.join(args.datadir, 'selected') + os.path.sep
reduced_dir = os.path.join(args.datadir, 'reduced') + os.path.sep
coadd_dir = os.path.join(args.datadir, 'coadd') + os.path.sep

for d in (data_dir, cal_dir, selected_dir, reduced_dir, coadd_dir):
    if not os.path.isdir(d):
        os.makedirs(d)

if not args.quiet:
    print('Data directories generated')
    print('cal_dir', cal_dir)
    print('selected_dir', selected_dir)
    print('reduced_dir', reduced_dir)
    print('coadd_dir', coadd_dir)


def str_list_to_list(str_list, delimiter=','):
    if str_list is None:
        return None
    assert isinstance(str_list, str)
    return str_list.split(delimiter)


def str_list_to_int_list(str_list, delimiter=','):
    return [int(i) for i in str_list_to_list(str_list, delimiter=delimiter)]


def preprocess_cli():
    if args.cams is None:
        if args.instrument == 'lmi':
            cams = [0]
        elif args.instrument == 'ratir':
            cams = [0, 1, 2, 3]
        elif args.instrument == 'rimas':
            cams = [0, 1]
        else:
            cams = [0, 1, 2, 3]
    else:
        cams = str_list_to_int_list(args.cams)
    preprocess_kwargs = {
        'instrument': args.instrument, 'workdir': data_dir, 'cams': cams,
        'auto': not args.noautoselect, 'reject_sat': not args.norejectsat, 'amin': args.amin, 'amax': args.amax,
        'save_select': not args.nosaveselect, 'noplot': args.noplot,
    }

    if not args.nobias:
        preprocess_kwargs['ftype'] = 'bias'
        preprocess_dict = choose.choose_calib(**preprocess_kwargs)
        master.mkmaster(args.instrument, preprocess_dict, preprocess_kwargs['ftype'], args.fmin, master_dir=cal_dir)

    if not args.noflat:
        preprocess_kwargs['ftype'] = 'flat'
        preprocess_dict = choose.choose_calib(**preprocess_kwargs)
        master.mkmaster(args.instrument, preprocess_dict, preprocess_kwargs['ftype'], args.fmin, master_dir=cal_dir)

    if not args.nodark:
        preprocess_kwargs['ftype'] = 'dark'
        preprocess_dict = choose.choose_calib(**preprocess_kwargs)
        master.mkmaster(args.instrument, preprocess_dict, preprocess_kwargs['ftype'], args.fmin, master_dir=cal_dir)

    if not args.noscienceselect:
        science_kwargs = preprocess_kwargs.copy()
        del science_kwargs['reject_sat']
        del science_kwargs['amin']
        del science_kwargs['amax']
        del science_kwargs['ftype']
        science_kwargs['targetdir'] = selected_dir
        choose.choose_science(**science_kwargs)


def autoproc_cli():
    autoproc_kwargs = {
        'datadir': selected_dir, 'imdir': reduced_dir, 'caldir': cal_dir, 'start': args.start, 'stop': args.stop,
        'only': None, 'step': args.step, 'mastersky': not args.nomastersky, 'skyflattarg': not args.noskyflattarg,
        'redo': args.redo, 'quiet': args.quiet, 'rmifiles': args.rmifiles,
        'customcat': args.customcat, 'customcatfilt': str_list_to_list(args.customcatfilt),
        'nogaia': args.nogaia,
    }
    autoproc(**autoproc_kwargs)

    for f in glob.glob(os.path.join(reduced_dir, 'coadd*.fits')):
        filename = os.path.basename(f)
        f_reduced_dir = os.path.join(reduced_dir, filename)
        f_coadd_dir = os.path.join(coadd_dir, filename)
        if not args.quiet:
            print('moving {0} to {1}'.format(f_reduced_dir, f_coadd_dir))
        shutil.move(f_reduced_dir, f_coadd_dir)


def photometry_cli():
    autoredux(coadd_dir)


def photopipe_all_cli():
    if not args.nopreprocess:
        preprocess_cli()
    if not args.noautoproc:
        autoproc_cli()
    if not args.nophotometry:
        photometry_cli()


def test_cli():
    base_path = os.path.dirname(os.path.abspath(__file__))
    test_path = os.path.join(base_path, 'test')
    zip_file = os.path.join(test_path, 'test.zip')
    copy_path = os.path.join(test_path, 'copy')
    zf = ZipFile(zip_file, 'r')
    zf.extractall(copy_path)
    test_files = os.listdir(copy_path)

    for f in test_files:
        original_path = os.path.join(copy_path, f)
        new_path = os.path.join(data_dir, f)
        # print('moving {} to {}'.format(original_path, new_path))
        shutil.move(original_path, new_path)

    photopipe_all_cli()


if args.cmd == 'bash':
    os.system('bash')
elif args.cmd == 'ipython':
    os.system('ipython')
elif args.cmd == 'pipeline_all':
    photopipe_all_cli()
elif args.cmd == 'preprocess':
    preprocess_cli()
elif args.cmd == 'autoproc':
    autoproc_cli()
elif args.cmd == 'photometry':
    photometry_cli()
elif args.cmd == 'test':
    test_cli()
else:
    print('Received command {}, which is not available'.format(args.cmd))
    parser.print_help()
