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
    'fullastrofail': '', 'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}


def autopipedefaults(pipevar=None):
    """
    NAME:
        autopipedefaults
    PURPOSE:
        Sets commonly used variables for pipeautoproc to use throughout each step
        Uses pipeautoproc.par to set variables, otherwise set to default values
        Saves in a dictionary
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                    but can be set to default)
    EXAMPLE:
        autopipedefaults(pipevar=inpipevar)
    """
    print('Setting pipeline parameters (DEFAULTS)')
    if pipevar is None:
        pipevar = inpipevar
    path = os.path.dirname(os.path.abspath(__file__))

    pipevar['pipeautopath'] = path
    sfile = path + '/pipeautoproc.par'

    if os.path.isfile(sfile):
        f = open(sfile, 'r')
        for line in f.readlines():
            line = ''.join(line.split())
            colpos = line.find(':')
            keyword = line[0:colpos]
            value = line[colpos + 1:]
            pipevar[keyword] = value
        f.close()

    if pipevar['refdatapath'] == '':
        pipevar['refdatapath'] = pipevar['pipeautopath'] + '/refdata'

    if pipevar['defaultspath'] == '':
        pipevar['defaultspath'] = pipevar['pipeautopath'] + '/defaults'

    if pipevar['imworkingdir'] != '' and not (os.path.exists(pipevar['imworkingdir'])):
        print('Creating imaging working directory: ', pipevar['imworkingdir'])
        os.makedirs(pipevar['imworkingdir'])

    pipevar['autoastrocommand'] = os.path.abspath(autoastrometry3.__file__)
    pipevar['getsedcommand'] = os.path.abspath(get_SEDs.__file__)


def autopipeprepare(pipevar=None):
    """
    NAME:
        autopipeprepare
    PURPOSE:
        Runs pipeprepare on every valid file and saves files with prefix 'p'.  Changes
        header with more manageable keywords and does bias/dark subtraction if bias/dark
        master exists (compares header keywords in files and bias/dark master)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeprepare(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    """

    print('PREPARE')
    if pipevar is None:
        pipevar = inpipevar
    # Looks for existing files in given data directory using prefix
    files = glob.glob(pipevar['datadir'] + pipevar['prefix'] + '*.fits')
    pfiles = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    if pipevar['verbose']:
        print('Found', len(files), 'files')

    # Finds any master bias files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    biasfiles = glob.glob(pipevar['imworkingdir'] + 'bias*')
    biascamera = []
    if len(biasfiles) > 0:
        for bfile in biasfiles:
            head = pf.getheader(bfile)
            camera = int(head['CAMERA'])
            biascamera += [camera]

    # Finds any master dark files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    darkfiles = glob.glob(pipevar['imworkingdir'] + 'dark*')
    darkcamera = []
    if len(darkfiles) > 0:
        for dfile in darkfiles:
            head = pf.getheader(dfile)
            camera = int(head['CAMERA'])
            darkcamera += [camera]

            # For each file (that doesn't have an existing p file or can be overwritten),
    # run pipeprepare on it with output file being saved into the imworkingdir,
    # will run bias subtraction if bias master available (checks based on how bias
    # file and data file are named
    for f in files:

        fileroot = os.path.basename(f)
        outnameim = pipevar['imworkingdir'] + 'p' + fileroot

        head = pf.getheader(f)
        camera = int(head['CAMERA'])

        try:
            bcamloc = biascamera.index(camera)
            biasfile = biasfiles[bcamloc]
        except:
            biasfile = None

        try:
            dcamloc = darkcamera.index(camera)
            darkfile = darkfiles[dcamloc]
        except:
            darkfile = None

        if (outnameim not in pfiles) or (pipevar['overwrite'] != 0):
            apd.pipeprepare(
                f, outname=outnameim, biasfile=biasfile, darkfile=darkfile, verbose=pipevar['verbose']
            )
        else:
            print('Skipping prepare. File already exists')


def autopipeimflatten(pipevar=None):
    """
    NAME:
        autopipeflatten
    PURPOSE:
        Flatten data using flat with matching filter name
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeflatten(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.flatpipeproc()
    """

    print('FLATTEN')
    if pipevar is None:
        pipevar = inpipevar
    # Finds prepared files and checks to see if there are any existing flattened files
    # Find flats in imworkingdir with name flat somewhere in a fits file name
    print(pipevar['imworkingdir'])
    files = glob.glob(pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
    ffiles = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    flats = glob.glob(pipevar['imworkingdir'] + '*flat*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    # If there are flats, then grab the filter from each of them,
    # otherwise end program
    flatfilts = []
    if len(flats) > 0:
        for flat in flats:
            head = pf.getheader(flat)
            head_filter = head['FILTER']
            flatfilts += [head_filter]
    else:
        print('No flats found for any filter!')
        return

    # Create outfile name and check to see if outfile already exists.  If it doesn't or
    # overwrite enabled then take filter from file and find where the flat filter matches
    # If no flats match filter, store in pipevar.flatfail, otherwise run flatproc on file
    for f in files:
        print(f)
        fileroot = os.path.basename(f)
        outnameim = pipevar['imworkingdir'] + 'f' + fileroot

        if (outnameim not in ffiles) or (pipevar['overwrite'] != 0):
            head = pf.getheader(f)
            head_filter = head['FILTER']

            try:
                flatfileno = flatfilts.index(head_filter)
            except:
                print('Flat field not found for ' + f + ' (filter=' + head_filter + ')')
                pipevar['flatfail'] += ' ' + f
                continue

            flatfile = flats[flatfileno]

            if pipevar['verbose']:
                print('Flattening', f, 'using', flatfile)

            apd.flatpipeproc(f, flatfile, flatminval=0.3)

        else:
            print('Skipping flatten. File already exists')

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')