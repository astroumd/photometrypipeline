import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.steps.autoproc_depend as apd
import datetime

inpipevar = {
    'autoastrocommand': 'autoastrometry', 'getsedcommand': 'get_SEDs', 'sexcommand': 'sex', 'swarpcommand': 'swarp',
    'rmifiles': 0, 'prefix': '', 'datadir': '', 'imworkingdir': '', 'overwrite': 0, 'verbose': 1, 'flatfail': '',
    'fullastrofail': '', 'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}


def autopipeimflatten(pipevar=None):
    """
    NAME:
        autopipeflatten
    PURPOSE:
        Flatten data using flat with matching band_filter name
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
    flats = apd.findcals(pipevar, 'flat*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    # If there are flats, then grab the band_filter from each of them,
    # otherwise end program
    flatfilts = []
    if len(flats) > 0:
        for flat in flats:
            head = pf.getheader(flat)
            head_filter = head['FILTER']
            flatfilts += [head_filter]
    else:
        print('No flats found for any band_filter!')
        return

    # Create outfile name and check to see if outfile already exists.  If it doesn't or
    # overwrite enabled then take band_filter from file and find where the flat band_filter matches
    # If no flats match band_filter, store in pipevar.flatfail, otherwise run flatproc on file
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
                print('Flat field not found for ' + f + ' (band_filter=' + head_filter + ')')
                pipevar['flatfail'] += ' ' + f
                continue

            flatfile = flats[flatfileno]

            if pipevar['verbose']:
                print('Flattening', f, 'using', flatfile)

            flatpipeproc(f, flatfile, flatminval=0.3)

        else:
            print('Skipping flatten. File already exists')

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')


def flatpipeproc(filename, flatname, flatminval=0, flatmaxval=0):
    """
    NAME:
        flatpipeproc
    PURPOSE:
        Checks if flat is same size as data, then divides for correct band_filter
    INPUTS:
        filename - name of FITS file, or array of filenames, or file w/list of filenames
        flatname - name of FITS master flat file
    OPTIONAL KEYWORDS:
        flatminval - if not set to 0 below this value will set to NaNs
        flatmaxval - if not set to 0 above this value will set to NaNs
    EXAMPLE:
        flatpipeproc(filename, flatname, flatminval=0.3)
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
        files = filename

    flat = pf.getdata(flatname)

    med = np.nanmedian(flat)
    if (med < 0.5) or (med > 2.0):
        print('Warning: flat is not normalized to one')

    for fname in files:
        f = pf.open(fname)
        data = f[0].data
        head = f[0].header
        f.close()

        if np.shape(data) != np.shape(flat):
            print(fname + ' could not be dark subtracted because it is not the same' +
                  ' size as the master dark, remove file to avoid confusion')
            return

            # Set values too low/high to NaNs
        if flatminval > 0:
            flat[flat < flatminval] = float('NaN')
        goodsignal = np.where(flat - 1.0 < 0.1)

        if flatmaxval > 0:
            flat[flat > flatminval] = float('NaN')

        # Divides out flattened field and adds keywords to header to show change
        fdata = data / flat

        head['FLATFLD'] = flatname
        skycts = np.nanmedian(fdata[goodsignal])
        head['SKYCTS'] = (skycts, 'Sky counts')

        try:
            head['CTRATE'] = (skycts / head['EXPTIME'], 'Sky counts per second')
        except:
            print('No EXPTIME keyword')

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by flatproc ' + date)

        fileroot = os.path.basename(fname)
        filedir = os.path.dirname(fname)
        outnameim = filedir + '/f' + fileroot

        apd.write_fits(outnameim, fdata, head)
