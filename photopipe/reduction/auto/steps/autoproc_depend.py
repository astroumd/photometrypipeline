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
            seepix = np.median(fwhmim[keep])        
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
