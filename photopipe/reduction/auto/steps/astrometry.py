import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.steps.autoproc_depend as apd
from photopipe.reduction.astrom import vlt_autoastrometry as autoastro

inpipevar = {
    'autoastrocommand': 'autoastrometry', 'getsedcommand': 'get_SEDs', 'sexcommand': 'sex', 'swarpcommand': 'swarp',
    'rmifiles': 0, 'prefix': '', 'datadir': '', 'imworkingdir': '', 'overwrite': 0, 'verbose': 1, 'flatfail': '',
    'fullastrofail': '',	'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}


def autopipeastrometry(pipevar=None):
    """
    NAME:
        autopipepipeastrometry
    PURPOSE:
        Calculate astrometry of image files to fix WCS coordinates (shift and rotation)
        in header. Using fast astrometry solver (vlt_autoastrometry.py) that using
        pair-distance matching and asterism matching.  Returns file with corrected WCS
        coordinates saved as 'a'+fitsfile. Run Scamp for additional astrometry
        corrections, twice, once for basic individual LOOSE correction, second correct all
        together.  Uses distortion of 3 as default, but uses 7 if distortion parameters
        high (i.e. RATIR H2RG)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeastrometry(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.astrometry, vlt_autoastrometry.py, scamp, sextractor
    FUTURE IMPROVEMENTS:
        Better distinction between first and second scamp run
    """

    print('ASTROMETRY')
    if pipevar is None:
        pipevar = inpipevar
    files = glob.glob(pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')

    # If no files, look for those that were not cosmic ray zapped
    if len(files) == 0:
        files = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    # Calculate relative astrometric solution
    for f in files:

        fileroot = os.path.basename(f)
        outfile = pipevar['imworkingdir'] + 'a' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping astrometry for ' + f + '. File already exists')
            continue

        head = pf.getheader(f)

        targ = head['TARGNAME']
        sat = head['SATURATE']
        # ascen = head['OBSRA']
        # decl = head['OBSDEC']
        ascen = head['CRVAL1']
        decl = head['CRVAL2']

        if 'flat' in targ:
            continue
        # filename, pixelscale=-1, pa=-999, inv=0, uncpa=-1, userra=-999, userdec=-999, minfwhm=1.5, maxfwhm=20,
        # maxellip=0.5, boxsize=-1, maxrad=-1, tolerance=0.010, catalog='', nosolve=0, overwrite=False, outfile='',
        # saturation=-1, quiet=False
        # cmd = 'python ' + pipevar['autoastrocommand'] + ' ' + f + ' -l ' + str(sat) + ' -r ' + str(ascen) + ' -d ' + str(decl)
        # cmd = 'python ' + pipevar['autoastrocommand'] + ' ' + f + ' -l ' + str(sat)
        # Run direct astrometry
        autoastro.autoastrometry(f, saturation=sat, userdec=decl, userra=ascen, quiet=(not pipevar['verbose']))
        # if pipevar['verbose'] > 0:
        #     os.system(cmd)
        #     print(cmd)
        # else:
        #     os.system(cmd + ' -q')
        if not os.path.isfile(outfile):
            pipevar['fullastrofail'] += ' ' + f

    if not os.path.isfile('astrom.param'):
        os.system('cp ' + pipevar['defaultspath'] + '/astrom.param .')
    if not os.path.isfile('astrom.conv'):
        os.system('cp ' + pipevar['defaultspath'] + '/astrom.conv .')
    if not os.path.isfile('default.sex'):
        os.system('cp ' + pipevar['defaultspath'] + '/default.sex .')
    if not os.path.isfile('default.missfits'):
        os.system('cp ' + pipevar['defaultspath'] + '/default.missfits .')
    if not os.path.isfile('scamp.conf'):
        os.system('cp ' + pipevar['defaultspath'] + '/scamp.conf .')

        # Calculate astrometry again using Scamp. First identify objects using sextractor,
    # then Scamp will solve by comparing reference catalog (currently set by default to
    # SDSS) to sources found by sextractor. Adds WCS corrections and second astrometry
    # parameters to header
    afiles = glob.glob(pipevar['imworkingdir'] + 'azsfp' + pipevar['prefix'] + '*.fits')

    # If no files, look for those that were not cosmic ray zapped
    if len(afiles) == 0:
        afiles = glob.glob(pipevar['imworkingdir'] + 'asfp' + pipevar['prefix'] + '*.fits')

    if len(afiles) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    afiletarg = []
    afilefilt = []

    for afile in afiles:
        head = pf.getheader(afile)
        afiletarg += [head['TARGNAME']]
        afilefilt += [head['FILTER']]

    atargets = set(afiletarg)
    afilters = set(afilefilt)

    afiletarg = np.array(afiletarg)
    afilefilt = np.array(afilefilt)
    afiles = np.array(afiles)

    for atarg in atargets:
        for afilt in afilters:

            thisatarget = np.where(np.logical_and(atarg == afiletarg, afilt == afilefilt))
            atfimages = afiles[thisatarget]

            try:
                head = pf.getheader(atfimages[0])
            except:
                continue

            # If scamp has already been run, skip
            try:
                test = head['ASTIRMS1']
                print('Skipping scamp astrometry for: ', atarg, afilt, ' Files already exist')
                print('head["ASTIRMS1"]: {}'.format(test))
                continue
            except:
                # Run sextractor to find sources, then use those catalogs to run scamp
                # with loose fitting constraints
                astrometry(atfimages, scamprun=1, pipevar=pipevar)

                # Do same thing again but with more stringent scamp parameters
                astrometry(atfimages, scamprun=2, pipevar=pipevar)

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits, zsfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')


def astrometry(atfimages, scamprun=1, pipevar=None):
    """
    NAME:
        astrometry
    PURPOSE:
        Run sextractor and scamp to refine astrometric solution
    INPUTS:
        atfimages - list of files to run through scamp
        scamprun  - the first run does a LOOSE run with distortion degree 1, any
                    other run will look for high distortion parameters, if it
                    finds it will use distortion degree 7, otherwise 3 (will also cut out
                    FLXSCALE on runs after 1)
    EXAMPLE:
        astrometry(atfimages, scamprun=2, pipevar=pipevar)
    FUTURE IMPROVEMENTS:
        Better difference between scamp runs.
    """
    acatlist = ''
    scat = {'sdss': 'SDSS-R7', 'tmpsc': '2MASS', 'tmc': '2MASS', 'ub2': 'USNO-B1'}
    for cfile in atfimages:
        head = pf.getheader(cfile)
        pixscale = head['PIXSCALE']
        sourcecat = head['ASTR_CAT']

        trunfile = os.path.splitext(cfile)[0]

        if pipevar['verbose'] > 0:
            sexcom = pipevar['sexcommand'] + ' -CATALOG_NAME ' + trunfile + \
                     '.cat -CATALOG_TYPE FITS_LDAC -FILTER_NAME astrom.conv ' + \
                     '-PARAMETERS_NAME astrom.param -DETECT_THRESH 2.0 ' + \
                     '-ANALYSIS_THRESH 2.0 -PIXEL_SCALE ' + str(pixscale) + \
                     ' ' + cfile
            print(sexcom)
        else:
            sexcom = pipevar['sexcommand'] + ' -CATALOG_NAME ' + trunfile + \
                     '.cat -CATALOG_TYPE FITS_LDAC -FILTER_NAME astrom.conv ' + \
                     '-PARAMETERS_NAME astrom.param -DETECT_THRESH 2.0 ' + \
                     '-ANALYSIS_THRESH 2.0 -VERBOSE_TYPE QUIET -PIXEL_SCALE ' + \
                     str(pixscale) + ' ' + cfile

        os.system(sexcom)

        if head['ASTR_NUM'] > 0:
            acatlist += ' ' + trunfile + '.cat'

        if sourcecat in scat:
            cat_u = scat[sourcecat]
        else:
            print('No valid catalogs available for SCAMP, check that vlt_autoastrometry.py ran correctly')
            return

    if scamprun == 1:
        loose = ' -MOSAIC_TYPE LOOSE'
        distdeg = 1
    else:
        loose = ' '
        try:
            distort = head['PV1_37']
            print("head['PV1_37']={}".format(distort))
            distdeg = 7
        except:
            distdeg = 3

    if pipevar['verbose'] > 0:
        scampcmd = "scamp -POSITION_MAXERR 0.2 -DISTORT_DEGREES " + str(distdeg) + \
                   loose + " -ASTREF_CATALOG " + cat_u + \
                   " -SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 " + \
                   "-CHECKPLOT_DEV NULL -WRITE_XML N -VERBOSE_TYPE FULL " + \
                   acatlist
        print(scampcmd)
    else:
        scampcmd = "scamp -POSITION_MAXERR 0.2 -DISTORT_DEGREES " + str(distdeg) + \
                   loose + " -ASTREF_CATALOG " + cat_u + \
                   " -SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 " + \
                   "-CHECKPLOT_DEV NULL -WRITE_XML N -VERBOSE_TYPE QUIET " + \
                   acatlist

    os.system(scampcmd)
    os.system('rm ' + acatlist)

    # Adds header information to file and delete extra files
    for cfile in atfimages:
        trunfile = os.path.splitext(cfile)[0]

        if pipevar['verbose'] > 0:
            os.system('missfits -WRITE_XML N ' + cfile)
        else:
            os.system('missfits -WRITE_XML N -VERBOSE_TYPE QUIET' + cfile)

        os.system('rm ' + trunfile + '.head ' + cfile + '.back')

        if scamprun != 1:
            him = pf.getheader(cfile)
            data = pf.getdata(cfile)
            del him['FLXSCALE']
            apd.write_fits(cfile, data, him)

