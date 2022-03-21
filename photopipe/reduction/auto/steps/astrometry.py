import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.steps.autoproc_depend as apd
from photopipe.reduction.astrom import vlt_autoastrometry as autoastro
from photopipe.reduction.dependencies import ldac
from photopipe.reduction.astrom.astrometrysources import GaiaAstrometry
from photopipe.reduction.astrom.vlt_autoastrometry import writeregionfile

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
        #cmd = 'python ' + pipevar['autoastrocommand'] + ' ' + f + ' -l ' + str(sat)
        # Run direct astrometry

        if pipevar['nogaia'] is True:
            # FIXME Joe: this breaks, as it is, you may want to look into it
            # Note that we bypass this step when using Gaia
            autoastro.autoastrometry(f, saturation=sat, userdec=decl,
                                     userra=ascen,
                                     quiet=(not pipevar['verbose']))

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
    afiles = glob.glob(pipevar['imworkingdir'] + 'azsfp' + pipevar['prefix'] +\
                       '*.fits')
    # If no files, look for those that were not cosmic ray zapped
    if len(afiles) == 0:
        afiles = glob.glob(pipevar['imworkingdir'] + 'asfp' + \
                           pipevar['prefix'] + '*.fits')

    if len(afiles) == 0 and pipevar['nogaia'] is False:
        # Create a copy of the files if using Gaia
        afiles = glob.glob(pipevar['imworkingdir'] + 'zsfp' + \
                           pipevar['prefix'] + '*.fits')
        if len(afiles) > 0:
            for af in afiles:
                os.system("cp " + af + " " + af.replace("zsfp", "azsfp"))
        else:
            afiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + \
                               pipevar['prefix'] + '*.fits')
            if len(afiles) > 0:
                for af in afiles:
                    os.system("cp " + af + " " + af.replace("sfp", "asfp"))

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
            except KeyError:
                # Run sextractor to find sources, then use those catalogs to run scamp
                # with loose fitting constraints
                astrometry(atfimages, scamprun=1, pipevar=pipevar,
                           nogaia=pipevar['nogaia'])

                # Do same thing again but with more stringent scamp parameters
                astrometry(atfimages, scamprun=2, pipevar=pipevar,
                           nogaia=pipevar['nogaia'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # sky-*.fits, sfp(PREFIX)*.fits, zsfp(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')


def astrometry(atfimages, scamprun=1, pipevar=None, nogaia=False):
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
        nogaia    - if False, the Gaia catalog will be downloaded and used for
                    the astrometric calibration. If True, the default scamp internal
                    catalogs will be used instead
    EXAMPLE:
        astrometry(atfimages, scamprun=2, pipevar=pipevar)
    FUTURE IMPROVEMENTS:
        Better difference between scamp runs.
    """
    acatlist = []
    gaiacatlist = []
    scat = {'sdss': 'SDSS-R7', 'tmpsc': '2MASS', 'tmc': '2MASS', 'ub2': 'USNO-B1'}
    for cfile in atfimages:
        head = pf.getheader(cfile)
        pixscale = head['PIXSCALE']
        try:
            sourcecat = head['ASTR_CAT']
        except KeyError:
            sourcecat = ''

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

        # FIXME what is this ASTR_NUM? Are we sure we want this condition?
        ##if head['ASTR_NUM'] > 0:
        acatlist.append(trunfile + '.cat')

        if nogaia is True:
            # Catalog to use, if not Gaia
            if sourcecat in scat:
                cat_u = scat[sourcecat]
            else:
                print('No valid catalogs available for SCAMP, check that \
vlt_autoastrometry.py ran correctly')
                return
        else:
            # Prepare the Gaia catalog
            ra_center = head["CRVAL1"]
            dec_center = head["CRVAL2"]
            box = np.max([head["NAXIS1"], head["NAXIS2"]]) * head["PIXSCALE"]
            # Add 10% to the box size for the catalog search
            box += box * 0.1
            gaiacat = cfile.replace(".fits", "_gaia.ldac")
            gaiacatlist.append(gaiacat)
            if pipevar['verbose'] > 0:
                print("Preparing the Gaia catalog for " + cfile)
            prepare_gaia_catalog(ra_center, dec_center, box, gaiacat)
            if pipevar['verbose'] > 0:
                print("Done, Gaia catalog prepared.")

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

    # Build up the scamp command depending on a number of conditions
    scampcmd = "scamp -POSITION_MAXERR 0.2 -DISTORT_DEGREES " + \
                str(distdeg) + loose + \
                "-SOLVE_PHOTOM N -SN_THRESHOLDS 3.0,10.0 " + \
                "-CHECKPLOT_DEV NULL -WRITE_XML N"
    if pipevar['verbose'] > 0:
        scampcmd += " -VERBOSE_TYPE FULL "
    else:
        scampcmd += " -VERBOSE_TYPE QUIET "
    if nogaia is True:
        # Use scamp internal catalogs
        scampcmd += " -ASTREF_CATALOG " + cat_u + " "
        scampcmd += " " + " ".join(acatlist)
        if pipevar['verbose'] > 0:
            print(scampcmd)

        os.system(scampcmd)
        os.system('rm ' + " ".join(acatlist))
    else:
        # Use gaia downloaded catalog
        scampcmd_base = scampcmd
        # Better results if we iterate over the images
        for gaiacat, acat in zip(gaiacatlist, acatlist):
            scampcmd = scampcmd_base + " -ASTREF_CATALOG FILE " + \
                    "-ASTREFCAT_NAME " + gaiacat + \
                    " -ASTREFCENT_KEYS RA_ICRS,DE_ICRS " + \
                    " -ASTREFERR_KEYS e_RA_ICRS,e_DE_ICRS -ASTREFMAG_KEY Gmag "
            scampcmd += " " + acat
            if pipevar['verbose'] > 0:
                print("Scamp run #" + str(scamprun) + " for image " + cfile)
                print(scampcmd)
            os.system(scampcmd)
        # Remove temporary catalogs
        os.system('rm ' + " ".join(acatlist))

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
            try:
                del him['FLXSCALE']
            except KeyError:
                pass
            apd.write_fits(cfile, data, him)


def prepare_gaia_catalog(ra, dec, field_size, out_file, write_region=True):
    """Query Gaia DR2 and save a table to do the astrometric calibration

    Parameters
    ----------
    ra float
        Right Ascension of the center of the field in degrees
    dec float
        Declination of the center of the field in degrees
    field_size float
        side of the field, in arcsec
    out_file str
        file name for the FITS LDAC table to be saved

    The table resulting from the query is saved in FITS LDAC format
    """
    # Create a GaiaAstrometry object
    gaia = GaiaAstrometry((ra, dec), field_size)
    # Gaia sources table for astrometry
    t_gaia = gaia.query_gaia_astrom()
    # Save the table in FITS LDAC format
    ldac.save_table_as_ldac(t_gaia, out_file, overwrite=True)
    # Write a region file
    if write_region is True:
        filename = out_file.replace(".ldac", ".reg")
        objlist = []
        for l in t_gaia:
            objlist.append(SourceCoords(l['RA_ICRS'], l['DE_ICRS']))
        writeregionfile(filename, objlist)

class SourceCoords:
    # Simple class for the source coordinates
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec 
