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
from photopipe.photometry.dependencies import get_SEDs


inpipevar = {
    'autoastrocommand': 'autoastrometry', 'getsedcommand': 'get_SEDs', 'sexcommand': 'sex', 'swarpcommand': 'swarp',
    'rmifiles': 0, 'prefix': '', 'datadir': '', 'imworkingdir': '', 'overwrite': 0, 'verbose': 1, 'flatfail': '',
    'fullastrofail': '',	'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}
def autopipestack(pipevar=None, customcat=None, customcatfilt=None):
    print('STACK')
    if pipevar is None:
        pipevar = inpipevar
    if customcatfilt is None:
        customcatfilt = []

    # If swarp configuration file ('default.swarp') does not exist, move swarp
    # output default configuration file
    if not os.path.isfile('default.swarp'):
        os.system(pipevar['swarpcommand'] + ' -d > default.swarp')

    qtcmd = 'True'
    quiet = 1
    if pipevar['verbose'] > 0:
        quiet = 0
        qtcmd = 'False'

    # Find files that have had zpoint performed on them, stop program if don't exist
    files = glob.glob(pipevar['imworkingdir'] + 't**sfp' + pipevar['prefix'] + '*.fits')
    print(pipevar['imworkingdir'] + 't**sfp' + pipevar['prefix'] + '*.fits')
    print(files)
    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    filetargs = []
    fileexpos = []
    filefilts = []
    fileairmv = []
    # filesatvs = []; filearms1 = []; filearms2 = []; filetime  = []
    filesatvs = []
    filetime = []

    # Grab information in the headers of zeropoint + Fluxscale corrected file and save to array
    for i, f in enumerate(files):
        head = pf.getheader(f)
        obstime = Time(head['DATE-OBS'], format='isot', scale='utc')

        # Strip target name of whitespace
        filetargs += [re.sub(r'\s+', '', head['TARGNAME'])]
        fileexpos += [head['EXPTIME']]
        filefilts += [head['FILTER']]
        fileairmv += [head['AIRMASS']]
        filesatvs += [head['SATURATE']]
        # filearms1 += [head['ASTRRMS1']]; filearms2 += [head['ASTRRMS2']]
        filetime += [obstime.jd]

    files = np.array(files)
    filetargs = np.array(filetargs)
    fileexpos = np.array(fileexpos)
    filefilts = np.array(filefilts)
    # filesatvs = np.array(filesatvs)
    fileairmv = np.array(fileairmv)
    # filearms1 = np.array(filearms1); filearms2 = np.array(filearms2)
    filetime = np.array(filetime)
    targets = set(filetargs)

    # Dictionary of corresponding columns for catalog file
    catdict = {'u': 2, 'g': 3, 'r': 4, 'i': 5, 'z': 6, 'y': 7,
               'B': 8, 'V': 9, 'R': 10, 'I': 11, 'J': 12, 'H': 13, 'K': 14,
               'ue': 15, 'ge': 16, 're': 17, 'ie': 18, 'ze': 19, 'ye': 20,
               'Be': 21, 'Ve': 22, 'Re': 23, 'Ie': 24, 'Je': 25, 'He': 26, 'Ke': 27, 'mode': 28}

    # Finds files with same target and the filters associated with this target
    for targ in targets:

        thistarget = np.where(filetargs == targ)
        if len(thistarget) == 0:
            continue

        thistargetfilts = set(filefilts[thistarget])

        # Find files that have the same target and same filter and store information
        # on the exposure times and airmass. Only use good Scamp astrometric fit files
        for thistargetfilter in thistargetfilts:
            stacki = (filetargs == targ) & (filefilts == thistargetfilter)
            # stacki = (filetargs == targ) & (filefilts == filter) &\
            #          (filearms1 < 2.0e-4) & (filearms1 > 5.0e-6) &\
            #          (filearms2 < 2.0e-4) & (filearms2 > 5.0e-6)


            if sum(stacki) == 0:
                continue

            stacklist = files[stacki]
            stackexps = fileexpos[stacki]
            stackairm = fileairmv[stacki]
            stacktime = filetime[stacki]
            pipevar['stackexps'] = stackexps

            medexp = np.median(stackexps)
            medair = np.median(stackairm)
            minair = min(stackairm)
            maxair = max(stackairm)
            totexp = sum(stackexps)
            nstack = len(stacklist)
            firsttime = Time(stacktime[0], format='jd', scale='utc').isot
            lasttime = Time(stacktime[-1], format='jd', scale='utc').isot
            medtime = Time(np.median(stacktime), format='jd', scale='utc').isot

            newtextslist = ' '.join(stacklist)

            stackcmd = pipevar['swarpcommand']

            # Keywords to carry through
            stackcmd += ' -COPY_KEYWORDS OBJECT,TARGNAME,FILTER,' + \
                        'INSTRUME,PIXSCALE,WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE '

            # Create output variables that will be used by SWarp
            outfl = pipevar['imworkingdir'] + 'coadd' + targ + '_' + re.sub(r'[^\w]', '', medtime) + '_' + \
                    thistargetfilter + '.fits'
            outwt = pipevar['imworkingdir'] + 'coadd' + targ + '_' + re.sub(r'[^\w]', '', medtime) + '_' + \
                    thistargetfilter + '.weight.fits'

            if pipevar['verbose'] > 0:
                stackcmd += ' -VERBOSE_TYPE NORMAL '
            else:
                stackcmd = stackcmd + ' -VERBOSE_TYPE QUIET '

            # Coadd with flux scale
            stackcmd = stackcmd + ' -SUBTRACT_BACK N -WRITE_XML N -IMAGEOUT_NAME ' + \
                       outfl + ' -WEIGHTOUT_NAME ' + outwt + \
                       ' -FSCALE_KEYWORD NEWFLXSC ' + newtextslist

            if pipevar['verbose'] > 0:
                print(stackcmd)

            os.system(stackcmd)
            head = pf.getheader(outfl)
            print(outfl)
            pixscl = head['PIXSCALE']

            try:
                # apd.findsexobj(outfl, 10.0, pipevar, pix=pixscl, aperture=20.0,
                apd.findsexobj(outfl, 1.5, pipevar, pix=pixscl, aperture=20.0, wtimage=outwt, quiet=quiet)
            except:
                sys.exit('Problem opening coadd fits file, may need to coadd in smaller bin size')

            head = pf.getheader(outfl)
            # cpsfdi = 1.56 * float(head['SEEPIX'])
            cpsfdi = 1.34 * float(head['SEEPIX'])

            # Run sextractor again on new coadd file
            apd.findsexobj(outfl, 3.0, pipevar, pix=pixscl, aperture=cpsfdi,
            #apd.findsexobj(outfl, 1.5, pipevar, pix=pixscl, aperture=cpsfdi,
                           wtimage=outwt, quiet=quiet)

            head = pf.getheader(outfl)

            coaddvars = np.loadtxt(outfl + '.stars', unpack=True)
            xim = coaddvars[1, :]
            yim = coaddvars[2, :]
            mag = coaddvars[3, :]
            mage = coaddvars[4, :]
            flag = coaddvars[5, :]
            elon = coaddvars[8, :]

            # astropy does not like SWarp PV keywords or unicode, temporarily delete
            headcopy = head.copy()
            for key in headcopy.keys():
                for comp_key in ['MJD-OBS', 'DATE-OBS', 'PV1_', 'PV2_', 'COMMENT', 'HISTORY']:
                    if key.startswith(comp_key):
                        try:
                            del head[key]
                        except KeyError as error:
                            print(error)

            w = wcs.WCS(head)
            wrd = w.all_pix2world(np.transpose([xim, yim]), 0)
            imfile = outfl + '.im'
            catfile = outfl + '.cat'

            # Save stars from image
            np.savetxt(imfile, np.transpose([wrd[:, 0], wrd[:, 1], mag]))

            # Filter name correction:
            if thistargetfilter == 'Z' or thistargetfilter == 'Y':
                thistargetfilter = thistargetfilter.lower()
            if thistargetfilter == 'YISH':
                thistargetfilter = 'y'

            if 'SDSS' in thistargetfilter:
                thistargetfilter = thistargetfilter[-1].lower()

            nocustomcat = False
            # If custom catalog provided, match the same objects as the *.im file
            # and create refmag (reference magnitudes) that have the same matched
            # indices
            if customcat is not None and thistargetfilter in customcatfilt:

                print('USING CUSTOM CATALOG')
                in_data = np.loadtxt(imfile)
                input_coords = in_data[:, :2]
                cat_data = np.loadtxt(customcat, skiprows=1)
                cat_coords = cat_data[:, :2]

                cat_matches, tmp = apd.identify_matches(input_coords, cat_coords)

                refmag = np.zeros(len(mag)) + 99
                mode = np.zeros(len(mag)) + -1
                for i, i_ind in enumerate(cat_matches):
                    if i_ind > 0:
                        # print input_coords[i], cat_coords[i_ind]
                        refmag[i] = cat_data[i_ind, catdict[thistargetfilter]]
                        mode[i] = 4

                # If no matching indices, run with the standard catalogs
                if sum(refmag < 90.0) == 0:
                    nocustomcat = True
            else:
                nocustomcat = True

            # If custom catalog not provided, catalog doesn't include filter, or
            # no objects from catalog found in image then
            # use get_SEDs.py to make catalog using 2MASS + (SDSS or APASS or USNOB1)
            if nocustomcat:
                # Create catalog star file
                # (python get_SEDs.py imfile filter catfile USNOB_THRESH alloptstars)
                sedcmd = 'python ' + pipevar['getsedcommand'] + ' ' + imfile + ' ' + \
                         thistargetfilter + ' ' + catfile + " 15 True " + qtcmd

                if pipevar['verbose'] > 0:
                    print(sedcmd)
                os.system(sedcmd)

                if not os.path.isfile(catfile):
                    #zpts += [float('NaN')]
                    continue

                # Read in catalog file
                cvars = np.loadtxt(catfile, unpack=True)
                refmag = cvars[catdict[thistargetfilter], :]
                mode = cvars[catdict['mode'], :]

            # Find catalog filter values and only cutoff values of actual detections
            goodind = (mode != -1) & (refmag < 90.0) & (flag < 8) & (elon <= 1.3)

            refmag = refmag[goodind]
            obsmag = mag[goodind]
            obserr = mage[goodind]
            ra_im = wrd[:,0][goodind]
            dec_im = wrd[:,1][goodind]
            #file = pipevar['imworkingdir'] + 'SATtest2_' + targ + '_' + thistargetfilter



            #find Sat sources
            satfile = pipevar['imworkingdir'] + 'SATcoords_' + targ + '_' + thistargetfilter + '.txt'
            var = np.loadtxt(satfile)
            if np.shape(var) != (0,):
                ra_sat = var[:, 0]
                dec_sat = var[:, 1]

                sat_ind = []
                thres = 0.001
                for i in range(len(ra_im)):
                    key = False
                    for j in range(len(ra_sat)):
                        if all(np.isclose([ra_sat[j],dec_sat[j]],[ra_im[i],dec_im[i]],atol=thres)):
                        #if (ra_sat[j] - thres < ra_im[i] < ra_sat[j] + thres) and (dec_sat[j] - thres < dec_im[i] < dec_sat[j] + thres):
                            key = True
                            continue
                    if key: sat_ind += [False]
                    else: sat_ind += [True]

                refmag = refmag[sat_ind]
                obsmag = obsmag[sat_ind]
                obserr = obserr[sat_ind]

            obswts = np.zeros(len(obserr))
            obskpm = np.zeros(len(obsmag))

            # Store magnitudes and weights (with minimum magnitude error of 0.01)
            for i in np.arange(len(obsmag)):
                obskpm[i] = obsmag[i]
                obswts[i] = 1.0 / (max(obserr[i], 0.01) ** 2)
                # if obserr[i] < 0.1:
                #     obskpm[i] = obsmag[i]
                #     obswts[i] = 1.0 / (max(obserr[i], 0.01) ** 2)

            czpts, cscats, crmss = apd.calc_zpt(
                np.array([refmag]), np.array([obskpm]), np.array([obswts]), sigma=1.0,
                plotter=pipevar['imworkingdir'] + 'zpt_COADD_' + targ + '_' + thistargetfilter + '.png')

            chead = pf.getheader(outfl)
            chead['SPIX'] = (cpsfdi, 'Final aperture size')
            chead['ABSZPT'] = (czpts[0] + 25.0, 'Absolute zeropoint from calc_zpt')
            chead['ABSZPTSC'] = (cscats[0], 'Robust scatter of absolute zeropoint')
            chead['ABSZPRMS'] = (crmss[0], 'RMS of absolute zeropoint')

            # Add summary of stack information to header
            chead['DATE1'] = (firsttime, 'First frame time')
            chead['DATEN'] = (lasttime, 'Last frame time')
            chead['DATE'] = (medtime, 'Median frame time')
            chead['NSTACK'] = nstack
            chead['AIRMASS'] = (medair, 'Median exposure airmass')
            chead['AIRMIN'] = (minair, 'Minimum exposure airmass')
            chead['AIRMAX'] = (maxair, 'Maximum exposure airmass')
            chead['EXPTIME'] = (medexp, 'Effective rescaled exposure time')
            chead['TOTALEXP'] = (totexp, 'Total summed integration time')
            chead['MAXEXP'] = (max(stackexps), 'Length of longest exposure')
            chead['MINEXP'] = (min(stackexps), 'Length of shortest exposure')

            for i, f in enumerate(stacklist):
                chead['STACK' + str(i)] = f

            cdata = pf.getdata(outfl)
            pf.update(outfl, cdata, chead)


            # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
            # sky-*.fits, sfp(PREFIX)*.fits, zsfp(PREFIX)*.fits files
            if pipevar['rmifiles'] != 0:
                os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + 'zsfp' + pipevar['prefix'] + '*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + 'azsfp' + pipevar['prefix'] + '*.fits')
                os.system('rm -f ' + pipevar['imworkingdir'] + 't**fp' + pipevar['prefix'] + '*.im')
                os.system('rm -f ' + pipevar['imworkingdir'] + 't**fp' + pipevar['prefix'] + '*.stars')
                os.system('rm -f ' + pipevar['imworkingdir'] + 't**fp' + pipevar['prefix'] + '*.cat')
