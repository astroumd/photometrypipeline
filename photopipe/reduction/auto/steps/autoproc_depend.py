import os
import shutil
import glob

import astropy.io.fits as pf
import numpy as np
import matplotlib.pyplot as plt

# Disable interactive mode
plt.ioff()


def findcals(pipevar, file_format_str):
    calfiles = glob.glob(os.path.join(pipevar['caldir'], file_format_str))
    if len(calfiles) == 0:
        calfiles = glob.glob(os.path.join(pipevar['datadir'], file_format_str))
    if len(calfiles) == 0:
        calfiles = glob.glob(os.path.join(pipevar['imworkingdir'], file_format_str))
    return calfiles


def write_fits(filename, data, header):
    try:
        pf.writeto(filename, data, header, overwrite=True)
        # print "write_fits.pf.writeto(filename, data, header, overwrite=True) worked!"
    except:
        temp_filename = filename + '.tmp'
        print("saving to temporary file: {}".format(temp_filename))
        pf.writeto(temp_filename, data, header, overwrite=True)
        try:
            os.remove(filename)
            print("deleted {}".format(filename))
        except:
            print("couldn't delete {}".format(filename))
            os.listdir(os.path.dirname(filename))
            pass
        print('attempting to overwrite {} --> {}'.format(temp_filename, filename))
        shutil.move(temp_filename, filename)
        print("write_fits.shutil.move(temp_filename, filename) worked!")


def findsexobj(filename, sigma, pipevar, masksfx=None, zeropt=25.0, maptype='MAP_WEIGHT',
               wtimage=None, fwhm=1.5, pix=0.3787, aperture=5.0, elong_cut=1.5, 
               quiet=0):
    """
    NAME:
        findsexobj
    PURPOSE:
        Finds sextractor objects with optional inputs. Estimates seeing from stars found. 
    INPUTS:
        file    - fits file to run sextractor on
        sigma   - detection threshold and analysis threshold for sextractor
        pipevar - pipeline parameters (typically set in autopipedefaults or autoproc)

    OPTIONAL KEYWORDS:
        masksfx   - text string identifier for sextractor CHECKIMAGE_NAME
        zeropt    - input value for sextractor MAG_ZEROPOINT
        wtimage   - file for input for sextractor WEIGHT_IMAGE
        fwhm      - input value for sextractor SEEING_FWHM
        pix       - input value for sextractor PIXEL_SCALE
        aperture  - input value for sextractor PHOT_APERTURES
        elong_cut - cutoff limit for FWHM calculation of elongation to eliminate non-stars
        quiet     - no output from sextractor if set
    EXAMPLE:
        findsexobj(file, 3.0, pipevar, aperture=20.0)
    DEPENDENCIES:
        sextractor
    FUTURE IMPROVEMENTS:
        More keywords to sextractor?
    """
    
    # Move necessary sextractor configuration files if they are not in current directory
    if not os.path.isfile('coadd.param'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.param .')
    if not os.path.isfile('coadd.conv'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.conv .') 
    if not os.path.isfile('coadd.config'): 
        os.system('cp ' + pipevar['defaultspath'] + '/coadd.config .') 
    if not os.path.isfile('default.nnw'): 
        os.system('cp ' + pipevar['defaultspath'] + '/default.nnw .')  

    if quiet > 0: 
        verbosetype = 'QUIET'
    else:
        verbosetype = 'NORMAL'
        
    # Run sextractor with given input parameters. Saves temp.cat as 
    # starfile, saves starmask, and calculates seeing from starlike objects. Saves 
    # necessary parameters to header
    if filename == '':
        return
    
    if not os.path.isfile(filename):
        return
    starfile = filename + '.stars'
        
    trunfile = os.path.splitext(filename)[0]
        
    sexcmd = pipevar['sexcommand'] + ' -c coadd.config -DETECT_THRESH ' +\
        str(sigma) + ' -ANALYSIS_THRESH ' + str(sigma) + ' -PHOT_APERTURES ' +\
        str(aperture) + ' -MAG_ZEROPOINT ' + str(zeropt) + ' -PIXEL_SCALE ' +\
        str(pix) + ' -SEEING_FWHM ' + str(fwhm) + ' -VERBOSE_TYPE ' + verbosetype
    
    if masksfx is not None:
        mskimg = trunfile + '_' + masksfx
        sexcmd += ' -CHECKIMAGE_TYPE OBJECTS' + ' -CHECKIMAGE_NAME ' + mskimg
    else:
        mskimg = "NOT USED"
        
    if wtimage is not None:
        sexcmd += ' -WEIGHT_TYPE '+maptype+' -WEIGHT_IMAGE ' + wtimage + ' '
        
    sexcmd += ' ' + filename
    print(sexcmd)
    if quiet == 0:
        print(sexcmd)
    os.system(sexcmd)
        
    if quiet == 0:
        print('mv -f test.cat ' + starfile)
    os.system('mv -f test.cat ' + starfile)
    
    num = 0    
    # Calculates seeing with starlike objects
    print(str(starfile))
    if os.path.isfile(starfile):
        variables = np.loadtxt(starfile, unpack=True)
        num = variables[0, :]
        flag = variables[5, :]
        elon = variables[8, :]
        fwhmim = variables[9, :]
        keep = (flag == 0) & (elon < elong_cut) & (fwhmim > fwhm) & (fwhmim < 20.0)

        if sum(keep) <= 1: 
            seepix = None
        else:
            seepix = np.nanmedian(fwhmim[keep])
    else:
        print('Failed to find Sextractor output file!')
        seepix = None
    head = pf.getheader(filename)
    
    if masksfx is not None:
        head['MASKNAME'] = (mskimg, "Object mask image from Sextractor")
    
    head['STARFILE'] = (starfile, "Objects file from Sextractor")
    head['ZEROPT'] = (zeropt, "Photometric zero-point used for Sextractor")
    if seepix is not None:
        head['SEEPIX'] = (seepix, "Estimated seeing from Sextractor objects (pix)")
    head['NSTARS'] = (len(num), "Estimated number of objects from Sextractor")
    
    data = pf.getdata(filename)
    write_fits(filename, data, head)
    
    # Removes config files after done
    os.system('rm -f coadd.param')
    os.system('rm -f coadd.conv')
    os.system('rm -f coadd.config')
    os.system('rm -f default.nnw')


def robust_scat(diff, wts, nobs, nstars, sigma):
    """
    NAME:
        robust_scat
    PURPOSE:
        Calculate robust scatter and set the weight of those above this limit to 0
    INPUTS:
        diff   - values to calculate robust scatter over
        wts    - weighting (0 is bad)
        nobs   - number of observations to iterate over
        nstars - number of stars to iterate over
        sigma  - sigma*robust scatter that is acceptable
    OUTPUTS:
        scats  - robust scatter of each observation
        rmss   - standard deviation (without bad weight points) of each observation
    EXAMPLE:
        robust_scat(diff, wts, 1, 12, 3)
    """

    scats = np.zeros(nobs)
    rmss = np.zeros(nobs)
    for i in np.arange(nobs):
        goodwts = np.where(wts[i, :] > 0)
        if len(goodwts[0]) == 0:
            continue
        gooddiff = diff[i, goodwts]

        # Median absolute deviation
        scat = 1.48 * np.nanmedian(abs(gooddiff - np.nanmedian(gooddiff)))
        for j in np.arange(nstars):
            if abs(diff[i, j] - np.nanmedian(gooddiff)) > (sigma * scat):
                wts[i, j] = 0
        scats[i] = scat
        rmss[i] = np.nanstd(gooddiff)
    return scats, rmss


def calc_zpt(catmag, obsmag, wts, sigma=3.0, plotter=None):
    """
    NAME:
        calc_zpt
    PURPOSE:
        Find zeropoint using robust scatter
    INPUTS:
        catmag  - 2d array with catalog magnitudes catmag[nobs,nstar]
        obsmag  - 2d array with observed magnitudes obsmag[nobs,nstar]
        wts     - 2d array with weights wts[nobs,nstar]
    OPTIONAL KEYWORDS:
        sigma   - sigma value for how far values can be from robust scatter value
        plotter - filename to save zeropoint plot
    OUTPUTS:
        z2     - zeropoint correction
        scats  - robust scatter of each observation
        rmss   - standard deviation (without bad weight points) of each observation
    EXAMPLE:
        zpt,scats,rmss = calc_zpt(catmag,obsmag,wts, sigma=3.0)
    """

    # Find difference between catalog and observed magnitudes
    diff = catmag - obsmag
    # print(np.shape(obsmag))

    keep = np.where(wts != 0)
    diff = diff[keep]
    obsmag = obsmag[keep]
    catmag = catmag[keep]
    wts = wts[keep]
    catmag_0 = np.copy(catmag)
    diff_0 = np.copy(diff)
    wts_0 = np.copy(wts)

    # # Sigma clip
    # med, sig = medclip(diff, 3.0, 2)
    # #print("Med: {}, Sigma: {}".format(med, sigma))
    # keep = np.where(abs(diff - med) < 3 * sig)
    # diff = diff[keep]
    # obsmag = obsmag[keep]
    # catmag = catmag[keep]
    # wts = wts[keep]
    # catmag_1 = np.copy(catmag)
    # diff_1 = np.copy(diff)
    # wts_1 = np.copy(wts)

    # Remove dim sources above provided threshold
    err_thresh = 0.1
    wt_thresh = 1 / (err_thresh ** 2)
    keep = np.where(wts > wt_thresh)
    diff = np.array([diff[keep]])
    obsmag = np.array([obsmag[keep]])
    catmag = np.array([catmag[keep]])
    wts = np.array([wts[keep]])
    catmag_2 = np.copy(catmag[0])
    diff_2 = np.copy(diff[0])
    wts_2 = np.copy(wts[0])

    # For each observation (i.e. frame) find the weighted difference and store zeropoint
    # and new magnitude with zeropoint correction
    nobs, nstars = np.shape(diff)
    z = []
    modmag = np.copy(obsmag)
    for i in np.arange(nobs):
        indz = sum(diff[i, :] * wts[i, :]) / sum(wts[i, :])
        z += [indz]
        modmag[i, :] = obsmag[i, :] + indz

    # Find difference of catalog and zeropoint corrected values. Remove any values with
    # weights set to 0 or lower.  Calculate robust scatter on these values.  If difference
    # with these weights is not within sigma*robust scatter then set weight to 0
    adiff1 = catmag - modmag
    scats, rmss = robust_scat(adiff1, wts, nobs, nstars, sigma)

    z2 = []
    # Recalculate zeropoint using corrected weights (difference still same)
    modmag2 = np.copy(obsmag)
    for i in np.arange(nobs):
        indz = sum(diff[i, :] * wts[i, :]) / sum(wts[i, :])
        z2 += [indz]
        modmag2[i, :] = obsmag[i, :] + indz

    adiff2 = catmag - modmag2
    # Recalculate robust scatter and rms scatter value on twice zeropoint corrected mags
    scats, rmss = robust_scat(adiff2, wts, nobs, nstars, sigma)

    if plotter is not None:
        plt.plot(catmag_0, diff_0, '*')
        plt.errorbar(catmag_0, diff_0, yerr=1.0 / np.sqrt(wts_0), fmt='.')
        plt.title('Before Robust Scatter')
        plt.ylabel('Difference between Catalog and Observed')
        plt.xlabel('Catalog magnitude')
        filename = plotter.replace('zpt', 'zptall')
        plt.savefig(filename)
        plt.clf()

        plt.hist(1.0 / np.sqrt(wts_0), bins=30)
        plt.title('Error Hist for All Sources')
        plt.ylabel('Counts')
        plt.xlabel('Error in Diff (Obs-Cat)')
        filename = plotter.replace('zpt', 'zptall_hist')
        plt.savefig(filename)
        plt.clf()

        # plt.plot(catmag_1, diff_1, '*')
        # plt.errorbar(catmag_1, diff_1, yerr=1.0 / np.sqrt(wts_1), fmt='.')
        # plt.title('Before Robust Scatter with Sigma Clipping')
        # plt.ylabel('Difference between Catalog and Observed')
        # plt.xlabel('Catalog magnitude')
        # filename = plotter.replace('zpt', 'zptall_sigmaclip')
        # plt.savefig(filename)
        # plt.clf()

        # plt.plot(catmag_2, diff_2, '*')
        # plt.errorbar(catmag_2, diff_2, yerr=1.0 / np.sqrt(wts_2), fmt='.')
        # plt.title('Before Robust Scatter with sigma clip and Dim Removal')
        # plt.ylabel('Difference between Catalog and Observed')
        # plt.xlabel('Catalog magnitude')
        # filename = plotter.replace('zpt', 'zptall_sigmaclip_dimrmv')
        # plt.savefig(filename)
        # plt.clf()

        keep = np.where(wts != 0)
        print(np.shape(catmag[keep]))
        plt.plot(catmag[keep], adiff2[keep], '*')
        plt.errorbar(catmag[keep], adiff2[keep], yerr=1.0 / np.sqrt(wts[keep]), fmt='.')
        plt.title('Post Robust Scatter')
        plt.ylabel('Difference between Catalog and Observed')
        plt.xlabel('Catalog magnitude')
        plt.savefig(plotter)
        plt.clf()

    return z2, scats, rmss


def parse_header_file(header_file, remove_comments=True, remove_history=True):
    fix_headers = ['CTYPE1', 'CTYPE2']
    with open(header_file, "r") as fin:
        lines = [line.strip() for line in fin]
    # Don't include END (or later lines)
    end = lines.index('END') if 'END' in lines else len(lines)
    lines = lines[:end]
    # Later pyfits versions changed this to a class method, so you can write
    # pyfits.Card.fromstring(text).  But in older pyfits versions, it was
    # a regular method.  This syntax should work in both cases.
    cards = [pf.Card().fromstring(line) for line in lines]
    header = pf.Header(cards)
    for card in fix_headers:
        header[card] = header[card].replace('TAN', 'TPV')
    if remove_comments:
        del header['COMMENT']
    if remove_history:
        del header['HISTORY']
    return header


def combine_header_and_fits(header_file, fits_file, remove_header_file=False):
    hdr = parse_header_file(header_file)
    with pf.open(fits_file) as fin:
        fin[0].header.update(hdr)
    if remove_header_file:
        os.remove(header_file)


def combine_header_and_fits_list(fits_file_paths, remove_header_file=False, header_extension='.head'):
    header_file_paths = [os.path.splitext(f)[0]+header_extension for f in fits_file_paths]
    for header_file_path, fits_file_path in zip(header_file_paths, fits_file_paths):
        combine_header_and_fits(header_file_path, fits_file_path, remove_header_file)
