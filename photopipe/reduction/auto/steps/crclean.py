import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.autoproc_depend as apd
from astropy import wcs
import re
import datetime
from astropy.time import Time
import sys
from scipy import interpolate
from photopipe.reduction.astrom import autoastrometry3
from photopipe.photometry.dependencies import get_SEDs


inpipevar = {
    'autoastrocommand': 'autoastrometry', 'getsedcommand': 'get_SEDs', 'sexcommand': 'sex', 'swarpcommand': 'swarp',
    'rmifiles': 0, 'prefix': '', 'datadir': '', 'imworkingdir': '', 'overwrite': 0, 'verbose': 1, 'flatfail': '',
    'fullastrofail': '',	'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}


def autopipecrcleanim(pipevar=None):
    """
    NAME:
        autopipecrcleanim
    PURPOSE:
        Removes cosmic rays
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipecrcleanim(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.cosmiczap
    FUTURE IMPROVEMENTS:
        Slow, alter cosmics.py?
        Get readnoise from header
    """

    print('CRCLEAN')
    if pipevar is None:
        pipevar = inpipevar
    # Find data that needs to be cosmic ray zapped
    files = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    # For each file check that objects meet count limits and exposure time
    # (i.e. short exposure time with lot of counts will be ignored), also targets that are
    # calibration files will be ignored.
    # Run cosmiczap on the files and have output files be 'z'+file plus weight files

    for f in files:

        head = pf.getheader(f)

        try:
            target = head['TARGNAME']
        except:
            print('Requires header keywords: TARGNAME. Check file.')
            continue

        if 'flat' in target.lower():
            continue
        if 'twilight' in target.lower():
            continue

        fileroot = os.path.basename(f)
        outfile = pipevar['imworkingdir'] + 'z' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping crzap for ' + f + '. File already exists')
            continue

        if pipevar['verbose'] > 0:
            print('Cleaning cosmic rays from', f)

        # Runs cosmics.py
        apd.cosmiczap(f, outfile, sigclip=6.0, maxiter=1, verbose=pipevar['verbose'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')