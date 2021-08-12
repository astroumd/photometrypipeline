"""
Purpose:    Creates master frames for each mtype (bias, dark, flat)
"""
import os
import pickle
import sys

# installed modules
import astropy.io.fits as pf
import numpy as np

# custom modules/functions
from photopipe.reduction.dependencies import astro_functs as af
from photopipe.instruments.specific_instruments import instrument_dict

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n))
# function to make FITS keywords to store file names of combined frames


def mkmaster(instrument, fn_dict, mtype, fmin=5, master_dir='./'):
    """
    PURPOSE:
        Make master calibration frames (bias, dark, flat)
        * currently no outlier rejection other than median combine
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        fn_dict     - dictionary output by choose_calib() containing organized
                      fits file names.  can also provide file name of pickled dictionary.
        mtype       - type of master frame. should be either 'flat', 'dark' or 'bias'
        fmin        - minimum number of files needed to make a master frame
    EXAMPLE:
        mkmaster('ratir', fn_dict=output from choose_calib(), mtype = bias, dark or flat name)
    FUTURE IMPROVEMENTS:
        - Better outlier rejection
    """
    instrum = instrument_dict[instrument]

    # check if input is a file name
    if type(fn_dict) is str:
        if fn_dict.split('.')[-1] == 'p':
            af.print_bold("Loading pickled dictionary from file.")
            fn_dict = pickle.load(open(fn_dict, 'rb'))
        else:
            af.print_err("Invalid pickle file extension detected. Exiting...")
            return

    # check for valid mtype
    if mtype not in [instrum.flatname, instrum.biasname, instrum.darkname]:
        af.print_err(
            "Error: valid arguments for mtype are {}, {} and {}. Exiting...".format(instrum.flatname, instrum.biasname,
                                                                                    instrum.darkname))
        return

    bands = fn_dict.keys()

    sorttype = 'BAND'

    if mtype in [instrum.biasname, instrum.darkname]:
        sorttype = 'CAMERA'

    d = os.getcwd().split('/')[-1]  # name of current directory
    af.print_head("\nMaking master {} frame in {}:".format(mtype, d))

    # work on FITs files for specified photometric bands
    for band in bands:

        # print current band
        af.print_under("\n{:^50}".format('{} {}'.format(band, sorttype)))

        first_file = fn_dict[band][0]
        ind_C = first_file.index('.fits') - 3
        cam = first_file[ind_C:ind_C + 2]
        cam_i = int(cam[1])

        # check if required files are present
        if instrum.has_cam_bias(cam_i):
            mbias_fn = '{}_{}.fits'.format(instrum.biasname, cam)
            mbias_fn = os.path.join(master_dir, mbias_fn)
        else:
            mbias_fn = None

        if instrum.has_cam_dark(cam_i):
            mdark_fn = '{}_{}.fits'.format(instrum.darkname, cam)
            mdark_fn = os.path.join(master_dir, mdark_fn)
        else:
            mdark_fn = None

        if mtype is not instrum.biasname:
            print(mbias_fn)
            print(mdark_fn)
            if mbias_fn is not None:
                if not os.path.exists(mbias_fn):
                    af.print_err(
                        'Error: {} not found.  Move master bias file to working directory to proceed.'.format(mbias_fn))
                    continue
            if mdark_fn is not None:
                if (mtype is instrum.flatname) and (not os.path.exists(mdark_fn)):
                    af.print_err(
                        'Error: {} not found.  Move master dark file to working directory to proceed.'.format(mdark_fn))
                    continue

        # check dictionary entries
        fns = fn_dict[band]
        if len(fns) < fmin:
            if len(fns) == 0:
                af.print_err(
                    'Error: no frames available to make master {} for {} {}.'.format(mtype, band, sorttype.lower()))
                continue
            else:
                temp = input(
                    af.bcolors.WARNING +
                    "Only {} frames available to make master {} for {} {}.  Continue? (y/n): ".format(
                        len(fns), mtype, band, sorttype.lower()) + af.bcolors.ENDC)
                if temp.lower() != 'y' and temp.lower() != 'yes':
                    af.print_warn("Skipping {}...".format(band))
                    continue

        # load calibration data
        hdu = pf.PrimaryHDU()
        filter_arr = []  # to check that all frames used the same filter
        exptime_arr = []  # to check that all frames have the same exposure time (where relevant)
        data_arr = []
        i = 0
        for fn in fns:
            print(fn)
            hdu.header[FITS_IN_KEY(i)] = fn  # add flat fn to master flat header
            hdulist = pf.open(fn)
            data_arr.append(hdulist[0].data)
            filter_arr.append(hdulist[0].header['FILTER'])
            exptime_arr.append(hdulist[0].header['EXPTIME'])
            i += 1
        data_arr = np.array(data_arr, dtype=np.float)

        # check that frames match
        for i in range(len(fns) - 1):
            if (filter_arr[i + 1] != filter_arr[0]) and (mtype is instrum.flatname):
                af.print_err("Error: cannot combine flat frames with different filters. Skipping {} {}...".format(
                    band, sorttype.lower()))
                continue
            if (exptime_arr[i + 1] != exptime_arr[0]) and (mtype is instrum.darkname):
                af.print_err(
                    "Error: cannot combine dark frames with different exposure times. Skipping {} {}...".format(
                        band, sorttype.lower()))
                continue
        if instrum.flatname:
            hdu.header['FILTER'] = filter_arr[0]  # add filter keyword to master frame
        if instrum.darkname:
            hdu.header['EXPTIME'] = exptime_arr[0]  # add exposure time keyword to master frame

        # add CAMERA header keyword
        hdu.header['CAMERA'] = cam_i  # add camera keyword to master frame

        # crop bias frames
        if mtype is instrum.biasname:
            data_arr = data_arr[(np.s_[:], instrum.slice(cam)[0], instrum.slice(cam)[1])]
        # crop dark frames and perform calculations
        elif mtype is instrum.darkname:
            data_arr = data_arr[(np.s_[:], instrum.slice(cam)[0], instrum.slice(cam)[1])]
            data_arr = (data_arr - pf.getdata(mbias_fn)) / hdu.header['EXPTIME']  # calculate dark current
        # crop flat frames and perform calculations
        elif mtype is instrum.flatname:
            if mbias_fn is not None:
                mbd = pf.getdata(mbias_fn)
                mbd = mbd
            if mdark_fn is not None:
                mdd = pf.getdata(mdark_fn)
                mdd = mdd

            if instrum.is_cam_split(cam_i):
                pass  # split data is already cropped
            else:
                data_arr = data_arr[(np.s_[:], instrum.slice(cam)[0], instrum.slice(cam)[1])]

            for i in range(len(exptime_arr)):
                if mbias_fn is not None:
                    data_arr[i] -= mbd
                if mdark_fn is not None:
                    data_arr[i] -= mdd * exptime_arr[i]
                data_arr[i] /= np.median(data_arr[i])

        # make master frame
        master = af.imcombine(data_arr, type='median').astype(np.float)

        # add master to hdu
        if mtype is instrum.flatname:
            hdu.data = master / np.median(master)  # normalize master flat
        else:
            hdu.data = master

        # save master to fits
        hdulist = pf.HDUList([hdu])
        try:
            hdulist.writeto('{}{}_{}.fits'.format(master_dir, mtype, band), overwrite=True)
        except IOError:
            os.mkdir(master_dir)
            hdulist.writeto('{}{}_{}.fits'.format(master_dir, mtype, band), overwrite=True)
