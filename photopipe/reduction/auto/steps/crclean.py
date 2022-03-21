import glob
import os
import astropy.io.fits as pf
import numpy as np
import datetime
from photopipe.reduction.auto import cosmics
import astroscrappy


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
        except KeyError:
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
        cosmiczap(f, outfile, sigclip=6.0, maxiter=3, verbose=pipevar['verbose'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')


def cosmiczap_slow(filename, outname, sigclip=6.0, maxiter=3, verbose=True):
    """
    NAME:
        cosmiczap_slow
    PURPOSE:
        Removes cosmic rays using Laplacian cosmic ray identification written in cosmics.py
    INPUTS:
        filename - file or list of files to be cosmic ray zapped
        outfile  - name of output file
    OPTIONAL KEYWORDS:
        sigclip  - sigma to clip
        maxiter  - maximum number of times to iterate loop
        verbose  - quiet?
    EXAMPLE:
        cosmiczap(filename, outname)
    DEPENDENCIES:
        crclean.py (described in http://arxiv.org/pdf/1506.07791v3.pdf)
    FUTURE IMPROVEMENTS:
        Read readnoise from header?
    """

    data, head = cosmics.fromfits(filename, verbose=False)

    gain = head['GAIN']
    c = cosmics.cosmicsimage(
        data, gain=gain, readnoise=18, sigclip=sigclip, sigfrac=0.5, objlim=5.0, verbose=False)

    tot = c.run(maxiter=maxiter, verbose=False)

    head['NPZAP'] = (tot, "Num. of pixels zapped by cosmiczap")
    date = datetime.datetime.now().isoformat()
    head.add_history('Processed by cosmiczap ' + date)

    if verbose:
        print('  Zapped %d total affected pixels (%.3f%% of total)' % (tot, tot * 100.0 / np.size(data)))

    cosmics.tofits(outname, c.cleanarray, head, verbose=False)


def cosmiczap(filename, outname, sigclip=6.0, maxiter=4, verbose=True):
    """
    NAME:
        cosmiczap
    PURPOSE:
        Removes cosmic rays using Laplacian cosmic ray identification using the
        astroscrappy package https://github.com/astropy/astroscrappy
        cite https://zenodo.org/record/1482019
        cite http://adsabs.harvard.edu/abs/2001PASP..113.1420V
    INPUTS:
        filename - file or list of files to be cosmic ray zapped
        outfile  - name of output file
    OPTIONAL KEYWORDS:
        sigclip  - sigma to clip
        maxiter  - maximum number of times to iterate loop
        verbose  - higher verbosity
    EXAMPLE:
        cosmiczap(filename, outname)
    DEPENDENCIES:
        crclean.py (described in http://arxiv.org/pdf/1506.07791v3.pdf)
    FUTURE IMPROVEMENTS:
        Read readnoise from header? Not present.
        It would be useful to have the median FWHM in the header too
    """
    data, head = cosmics.fromfits(filename, verbose=False)

    # Gain, saturation level from the header
    gain = head['GAIN']
    satlevel = head['SATURATE']

    # Run astroscrappy
    mask, data_clean = astroscrappy.detect_cosmics(data, inmask=None, inbkg=None, invar=None, sigclip=sigclip, sigfrac=0.3, objlim=5.0, gain=gain, readnoise=18, satlevel=satlevel, niter=maxiter, sepmed=True, cleantype='meanmask', fsmode='median', psfmodel='gauss', psffwhm=2.5, psfsize=7, psfk=None, psfbeta=4.765, verbose=verbose)

    # Total number of zapped pixels
    tot = np.count_nonzero(mask)
    head['NPZAP'] = (tot, "Num. of pixels zapped by cosmiczap")

    # Time stamp
    date = datetime.datetime.now().isoformat()
    head.add_history('Processed by cosmiczap_astroscrappy ' + date)

    if verbose:
        print('  Zapped %d total affected pixels (%.4f%% of total)' % (tot, tot * 100.0 / np.size(data)))

    # Write to disk
    cosmics.tofits(outname, data_clean, head, verbose=False)
