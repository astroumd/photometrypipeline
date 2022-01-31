import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.steps.autoproc_depend as apd
from photopipe.reduction.astrom import vlt_autoastrometry as autoastro
from photopipe.SEDs import get_SEDs

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

    pipevar['pipeautopath'] = path + "/../"
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

    pipevar['autoastrocommand'] = os.path.abspath(autoastro.__file__)
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
    biasfiles = apd.findcals(pipevar, 'bias*.fits')

    biascamera = []
    if len(biasfiles) > 0:
        for bfile in biasfiles:
            head = pf.getheader(bfile)
            camera = int(head['CAMERA'])
            biascamera += [camera]
    else:
        print('Did not find any BIAS files! Check your data directory path!')


    # Finds any master dark files and filter name from header keyword
    # Assumes camera name is in header under CAMERA
    darkfiles = apd.findcals(pipevar, 'dark*.fits')

    darkcamera = []
    if len(darkfiles) > 0:
        for dfile in darkfiles:
            head = pf.getheader(dfile)
            camera = int(head['CAMERA'])
            darkcamera += [camera]
    else:
        print('Did not find any DARK files! Check your data directory path!')

    # Flat files are already bias/dark subtracted
    flatfiles = apd.findcals(pipevar, 'flat*.fits')

    # For each file (that doesn't have an existing p file or can be overwritten),
    # run pipeprepare on it with output file being saved into the imworkingdir,
    # will run bias subtraction if bias master available (checks based on how bias
    # file and data file are named
    for f in files:
        print(f)

        if f in biasfiles:
            continue
        elif f in darkfiles:
            continue
        elif f in flatfiles:
            continue

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
            pipeprepare(f, outname=outnameim, biasfile=biasfile, darkfile=darkfile, verbose=pipevar['verbose'])
        else:
            print('Skipping prepare. File already exists')

def pipeprepare(filename, outname=None, biasfile=None, darkfile=None, verbose=1):

    """
    NAME:
        pipeprepare
    PURPOSE:
        Adds additional header keywords needed later on in the pipeline and removes
        unnecessary header keywords by looking through a list of mandatory keywords.
        Also runs bias and dark subtraction for filters with an existing master bias/dark
        (CCDs).  The prepared images are written to disk with outname
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames
    OPTIONAL KEYWORDS:
        outname  - specify output file to write to disk
        biasfile - full name (including path) of master bias
        darkfile - full name (including path) of master dark
        verbose  - print out comments
    EXAMPLE:
        pipeprepare(filename, outname=outname, biasfile=biasfile, darkfile=darkfile, verbose=1)
    DEPENDENCIES:
        autoproc_depend.pipeprepare()
    """

    # ------ Process input filenames(s) ------

    # Check for empty filename
    if len(filename) == 0:
        print('No filename specified')
        return

    # If string, check if a file of items or wildcards
    # otherwise store all files
    if isinstance(filename, str):
        fileext = os.path.splitext(filename)[1][1:]

        files = [filename]

        if fileext in ['cat', 'lis', 'list', 'txt']:
            f = open(filename, 'r')
            files = f.read().splitlines()
            f.close()

        if '?' in filename or '*' in filename:
            files = glob.glob(filename)
            if len(files) == 0:
                print('Cannot find any files matching ', filename)
                return
    else:
        files = [filename]

    # ------ Read data and process header information ------
    for pipe_file in files:
        f = pf.open(pipe_file)
        head = f[0].header
        data = f[0].data
        f.close()

        # If these keys exist keep, otherwise delete all extraneous keywords
        mandatorykey = [
            'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
            'HISTORY', 'DATE-OBS', 'EXPTIME', 'INSTRUME',
            'LATITUDE', 'LONGITUD', 'BINNING', 'BINY', 'BINX',
            'CAMERA', 'TARGNAME', 'UTC', 'OBJECT', 'OBJNAME', 'AIRMASS',
            'GAIN', 'SATURATE', 'PIXSCALE', 'FILTER', 'WAVELENG',
            'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
            'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CTYPE1', 'CTYPE2',
            'PV1_1', 'PV2_1', 'PV1_17', 'PV2_17', 'PV1_19', 'PV2_19', 'PV1_21', 'PV2_21',
            'PV1_31', 'PV2_31', 'PV1_33', 'PV2_33', 'PV1_35', 'PV2_35', 'PV1_37', 'PV2_37', 'OBSRA', 'OBSDEC'
        ]

        # Finds list of unnecessary keywords, then deletes extraneous
        newhead = head.copy()
        for oldkey in head.keys():
            if oldkey not in mandatorykey:
                try:
                    del newhead[oldkey]
                except KeyError:
                    pass

        # Create binary array of saturated pixels, save for later use
        find_sats(outname, data, newhead)

        # If biasfile keyword set subtract master bias from current file with given master bias file
        # If they are not the same size, quit program without saving with preparation prefix (will not move
        # on in following processing steps)
        if biasfile is not None:
            bias = pf.getdata(biasfile)

            if np.shape(data) != np.shape(bias):
                print(pipe_file + ' could not be bias subtracted because it is not the same' +
                      ' size as the master bias, remove file to avoid confusion')
                return

            if verbose > 0:
                print('    bias subtracting')

            newdata = data - bias

            # If darkfile keyword set subtract master dark from current file with given master dark file
            # If they are not the same size, quit program without saving with preparation prefix (will not move
            # on in following processing steps)
            if darkfile is not None:
                dark = pf.getdata(darkfile) * newhead['EXPTIME']

                if np.shape(data) != np.shape(dark):
                    print(' ')
                    print(pipe_file + ' could not be dark subtracted because it is not the same' +
                          ' size as the master dark, remove file to avoid confusion')
                    return

                if verbose > 0:
                    print('    dark subtracting')

                newdata = newdata - dark
            else:
                print(pipe_file, 'could not be dark subtracted because the master dark file was not provided')
        else:
            newdata = data

        # Write changes to disk
        apd.write_fits(outname, newdata, newhead)

        if verbose > 0:
            print(pipe_file, '-> ', outname)

def find_sats(fname, data, header):
    sat = header['SATURATE']
    saturated = np.where(data > sat, 0, 1)
    print("# of Saturated Pixels: {}".format(np.sum(saturated)))
    fileroot = os.path.basename(fname)
    filedir = os.path.dirname(fname)
    outname = filedir + "/" + fileroot.replace("p", "SAT_", 1)
    apd.write_fits(outname, saturated, header)