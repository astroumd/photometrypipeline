import glob
import os
import astropy.io.fits as pf
import numpy as np
import photopipe.reduction.auto.autoproc_depend as apd
from astropy import wcs
import re
import matplotlib.pyplot as plt
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


def autopipezpoint(pipevar=None, customcat=None, customcatfilt=None):
    """
    NAME:
        autopipepipestack
    PURPOSE:
        Does zeropoint correction on each individual frame using sextractor and get_SEDs.
        Creates flux scale (newflxsc) from how close to median of zeropoint values.  Uses
        flux scale to stack images in Swarp (has moved bad zeropoint values and bad newflxsc
        values to marked folders - badzptfit/ and badflxsc/) and calculates absolute zeropoint
        correction of coadd. Saves zeropoint plot as zpt_(FILTER).ps
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipestack(pipevar=inpipevar)
    DEPENDENCIES:
        SWarp, get_SEDs, calc_zpt, findsexobj (sextractor)
    """

    print('ZPOINT')
    if pipevar is None:
        pipevar = inpipevar
    if customcatfilt is None:
        customcatfilt = []

    os.system('export CDSCLIENT=http')  # Fix for problem with timeout with CDSCLIENT

    qtcmd = 'True'
    quiet = 1
    if pipevar['verbose'] > 0:
        quiet = 0
        qtcmd = 'False'

    # If swarp configuration file ('default.swarp') does not exist, move swarp
    # output default configuration file
    if not os.path.isfile('default.swarp'):
        os.system(pipevar['swarpcommand'] + ' -d > default.swarp')

    # Find files that have had astrometry performed on them, stop program if don't exist
    files = glob.glob(pipevar['imworkingdir'] + 'a*sfp' + pipevar['prefix'] + '*.fits')
    print(pipevar['imworkingdir'] + 'a*sfp' + pipevar['prefix'] + '*.fits')
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

    # Grab information in the headers of astrometry corrected file and save to array
    for i, f in enumerate(files):
        head = pf.getheader(f)
        data = pf.getdata(f)
        obstime = Time(head['DATE-OBS'], format='isot', scale='utc')

        # Strip target name of whitespace
        filetargs += [re.sub(r'\s+', '', head['TARGNAME'])]
        fileexpos += [head['EXPTIME']]
        filefilts += [head['FILTER']]
        fileairmv += [head['AIRMASS']]
        filesatvs += [head['SATURATE']]
        # filearms1 += [head['ASTRRMS1']]; filearms2 += [head['ASTRRMS2']]
        filetime += [obstime.jd]
        fileroot = os.path.basename(f)
        print(fileroot)
        zfile = pipevar['imworkingdir'] + 't' + fileroot
        apd.write_fits(zfile, data, head)

    files = glob.glob(pipevar['imworkingdir'] + 't**sfp' + pipevar['prefix'] + '*.fits')
    print(pipevar['imworkingdir'] + 't**sfp' + pipevar['prefix'] + '*.fits')
    print(files)

    files = np.array(files)
    filetargs = np.array(filetargs)
    filefilts = np.array(filefilts)
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

            zpts = []
            it_num = 0 #for testing

            # Find stars for each individual frame and try to find matches with coadded
            # frame with each frame optimized with PSF size

            for sfile in stacklist:
                print("working")
                head = pf.getheader(sfile).copy()
                ipixscl = head['PIXSCALE']
                apd.findsexobj(sfile, 3.0, pipevar,pix=ipixscl,aperture=20.0, quiet=quiet) ### Change this
                #apd.findsexobj(sfile, 1.5, pipevar, pix=ipixscl, aperture=20.0, quiet=quiet)
                starfile = sfile + '.stars'

                svars = np.loadtxt(starfile, unpack=True)
                xim = svars[1, :]
                yim = svars[2, :]
                mag = svars[3, :]
                mage = svars[4, :]
                flag = svars[5, :]
                elon = svars[8, :]
                fwhm = svars[10, :]

                # astropy does not like SWarp PV keywords or unicode, temporarily delete
                headcopy = head.copy()
                for key in headcopy.keys():
                    for comp_key in ['PV1_', 'PV2_', 'COMMENT', 'HISTORY']:
                        if key.startswith(comp_key):
                            try:
                                del head[key]
                            except KeyError as error:
                                print(error)

                    # if any(x in key for x in ['PV1_', 'PV2_', 'COMMENT', 'HISTORY']):
                    #     try:
                    #         del head[key]
                    #     except KeyError as error:
                    #         print(error)

                w = wcs.WCS(head)
                wrd = w.all_pix2world(np.transpose([xim, yim]), 0)
                imfile = sfile + '.im'
                catfile = sfile + '.cat'

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
                        zpts += [float('NaN')]
                        continue

                    # Read in catalog file
                    cvars = np.loadtxt(catfile, unpack=True)
                    refmag = cvars[catdict[thistargetfilter], :]
                    mode = cvars[catdict['mode'], :]

                # Find catalog filter values and only cutoff values of actual detections
                goodind = (mode != -1) & (refmag < 90.0) & (flag < 8) & (elon <= 1.5)

                refmag = refmag[goodind]
                obsmag = mag[goodind]
                obserr = mage[goodind]
                obswts = np.zeros(len(obserr))
                obskpm = np.zeros(len(obsmag))

                # Store magnitudes and weights (with minimum magnitude error of 0.01)
                for i in np.arange(len(obsmag)):
                    if obserr[i] < 0.1:
                        obskpm[i] = obsmag[i]
                        obswts[i] = 1.0 / (max(obserr[i], 0.01) ** 2)

                if len(refmag) > 0 and len(obskpm) > 0 and len(obswts) > 0:
                    #zpt, scats, rmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]),np.array([obswts]), sigma=3.0)
                    zpt, scats, rmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]),
                                                    np.array([obswts]), sigma=3.0,
                                                    plotter=pipevar['imworkingdir'] + 'zpt_' + targ + '_' + thistargetfilter + '_{}.png'.format(it_num))
                    it_num += 1 #testing
                    # Reload because we had to remove distortion parameters before
                    head = pf.getheader(sfile)
                    data = pf.getdata(sfile)
                    head['ABSZPT'] = (zpt[0] + 25.0, 'Relative zeropoint from calc_zpt')
                    head['ABSZPTSC'] = (scats[0], 'Robust scatter of relative zeropoint')
                    head['ABSZPRMS'] = (rmss[0], 'RMS of relative zeropoint')

                    pf.update(sfile, data, head)
                    zpts += zpt
                else:
                    zpts = [np.inf]
                    it_num += 1 #tesing

            # Move files with bad zeropoint calculations to folder 'badzptfit'
            # and do not use those frames
            zpts = np.array(zpts)
            goodframes = np.isfinite(zpts)
            badframes = ~np.isfinite(zpts)

            if len(zpts[badframes]) != 0:
                if not os.path.exists(pipevar['imworkingdir'] + '/badzptfit'):
                    os.makedirs(pipevar['imworkingdir'] + '/badzptfit')
                for f in stacklist[badframes]:
                    os.system('mv ' + f + ' ' + pipevar['imworkingdir'] + '/badzptfit/')
                zpts = zpts[goodframes]
                newstack = stacklist[goodframes]
            else:
                newstack = stacklist

            badnewflxsc = []
            # Add relative zeropoint values to headers and calculate flux scale.
            # Remove unphysical fluxscale files
            medzp = np.median(zpts)
            for i, f in enumerate(newstack):
                head = pf.getheader(f)
                head['NEWFLXSC'] = (1.0 / (10.0 ** ((zpts[i] - medzp) / 2.5)),
                                    'Flux scaling based on median zp')

                if 1.0 / (10.0 ** ((zpts[i] - medzp) / 2.5)) < 0.1:
                    badnewflxsc += [f]

                data = pf.getdata(f)
                pf.update(f, data, head)

            removedframes = []
            # Removes files that have bad newflxsc values and removes from stack list
            if len(badnewflxsc) > 0:
                if not os.path.exists(pipevar['imworkingdir'] + '/badflxsc'):
                    os.makedirs(pipevar['imworkingdir'] + '/badflxsc')

                for ibad in badnewflxsc:
                    os.system('mv ' + ibad + ' ' + pipevar['imworkingdir'] + 'badflxsc/')

                removedframes += badnewflxsc

            if len(removedframes) > 0:
                print('Removed frames with bad zeropoint fits: ')
                print(removedframes)

            #### Testing ####
            print(np.shape(zpts))
            plt.hist(zpts.flatten())
            # for row in zpts:
            #     xx = np.zeros(len(row))+it_num
            #     yy = np.array(row)+25.0
            #     plt.plot(xx,yy, '*')
            #     it_num += 1
            #     #plt.errorbar(xx, yy, yerr=1.0 / np.sqrt(wts_init), fmt='.')
            plt.title('Abs. Zeropoint Histogram')
            plt.xlabel('Abs. Zeropoint')
            plt.ylabel('Counts')
            plt.savefig(pipevar['imworkingdir'] + 'zpoint_values.png')
            plt.clf()