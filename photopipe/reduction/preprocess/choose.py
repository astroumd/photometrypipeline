"""
Purpose:    Automatically selects individual types of frames from Image Directory
"""
import os
from matplotlib.patches import Rectangle
import shutil
from glob import glob
import datetime
from six.moves.collections_abc import Iterable

# installed modules
import astropy.io.fits as pf
import matplotlib.pylab as pl
import numpy as np

# custom modules/functions
# from photopipe.reduction.dependencies.zscale import zscale
from photopipe.reduction.dependencies import astro_functs as af
from photopipe.instruments.specific_instruments import instrument_dict
from photopipe.reduction.preprocess import json_helper

# Preprocessing constants
FITS_IN_KEY = lambda n: 'IMCMB{:03}'.format(int(n))
# function to make FITS keywords to store file names of combined frames


def choose_calib(instrument, ftype, workdir='.', cams=(0, 1, 2, 3), auto=False, reject_sat=True, amin=0.2, amax=0.8,
                 save_select=True, figsize=(8, 5), noplot=False, yes=False):
    """
    NAME:
        choose_calib
    PURPOSE:
        Either auto-select or display calibration images for user verification
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        ftype       - type of calibration frames (ex. 'flat', 'bias', 'dark')
        workdir     - directory where function executed
        cams        - camera numbers (default is all)
        auto        - automated frame selection. If 'bias', will select all, if 'flat'
                      will select non-saturated frames with sufficient counts
        reject_sat  - reject frames with saturated pixels
        amin        - minimum fraction of saturation value for median (automated)
        amax        - maximum fraction of saturation value for median (automated)
        save_select - save dictionary of selected frames to python json file
    EXAMPLE:
        file_dict = choose_calib('ratir', ftype = bias, dark or flat name,
            workdir = 'path/to/data/', cams = [#,#,...])
        *** call mkmaster using dictionary or json file ***
    """
    instrum = instrument_dict[instrument]

    if auto and (ftype is instrum.flatname) and not yes:
        temp = input(
            af.bcolors.WARNING + "Warning: automated selection of flats is not recommended! Continue? (y/n): " +
            af.bcolors.ENDC
        )
        if (temp.lower() != 'y') and (temp.lower() != 'yes'):
            af.print_bold("Exiting...")
            return

    # check for non-list camera argument
    if not isinstance(cams, Iterable):
        cams = [cams]  # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return

    if not noplot:
        pl.ion()  # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)
    # wrkdir = os.getcwd()  # TODO

    d = os.getcwd().split('/')[-1]  # name of current directory
    if not auto:
        af.print_head("\nDisplaying {} frames in {} for selection:".format(ftype, d))
    else:
        af.print_head("\nAutomatically selecting {} frames in {}:".format(ftype, d))

    # dictionary to store selected fits files by camera or band_filter
    fits_list_dict = {}

    # if len(fits_check) == 0:
    if len(glob(instrum.original_file_format())) != 0:
        print("Changing file names:")
        instrum.change_file_names(glob(instrum.original_file_format()))

        fits_check_2 = glob('????????T??????C??.fits')

        if len(fits_check_2) == 0:
            af.print_err("Error: no files with correct format after file name changes")
            return

    # open figure for images if not auto
    # if not auto and not noplot:
    fig = pl.figure(figsize=figsize)

    # work on FITs files for specified cameras
    for cam_i in cams:

        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))

        # check for valid calibration request
        if not instrum.has_cam_bias(cam_i) and ftype is instrum.biasname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  " +
                          "Skipping...".format(cam_i, instrum.biasname))
            continue
        if not instrum.has_cam_dark(cam_i) and ftype is instrum.darkname:
            af.print_warn("Warning: Camera C{} does not have {} frames.  " +
                          "Skipping...".format(cam_i, instrum.darkname))
            continue

        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i,
                                                            instrum.ftype_post[ftype]))

        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:

            fits_id = fits_fn.split('.')[0]  # fits file name with extention removed
            print('{}'.format(fits_fn))

            # open data
            hdulist = pf.open(fits_fn, mode='update')
            im = hdulist[0].data
            h = hdulist[0].header

            sat_pt = instrum.get_cam_sat(h, cam_i)

            if reject_sat:
                if np.any(im == sat_pt):
                    af.print_warn("Warning: saturated pixels in frame.  Skipping frame {}.".format(fits_fn))
                    continue

            # check if instrument camera is split, if it is make sure correct specs being used
            if instrum.is_cam_split(cam_i):
                print('\t* Split band_filter used')
            else:
                if (instrum.get_filter(h, 'C{}'.format(cam_i)) not in instrum.possible_filters()) and (
                        ftype is instrum.flatname):
                    af.print_warn("Warning: invalid band_filter detected.  Skipping {} band.".format(
                        instrum.get_filter(h, 'C{}'.format(cam_i))))
                    continue

                if ftype is instrum.flatname:
                    print('\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}'.format(cam_i))))

            h = instrum.change_header_keywords(h, 'C{}'.format(cam_i))

            # return quick image summary
            [im1, m, s, sfrac] = image_summary(im, sat_pt, cam_i, instrum, split=instrum.is_cam_split(cam_i))

            if instrum.is_cam_split(cam_i):
                [m1, m2] = m
                [s1, s2] = s
                [im1, im2] = im1
                [sfrac1, sfrac2] = sfrac
            else:
                m1 = m2 = m
                s1 = s2 = s
                im1 = im2 = im1
                sfrac1 = sfrac2 = sfrac

                # if auto select then find values with correct ranges
            if auto:

                # all bias and dark frames are selected
                if ftype in [instrum.biasname, instrum.darkname]:
                    addtodict(dictionary=fits_list_dict, key='C{}'.format(cam_i), value='{}{}'.format(workdir, fits_fn))

                # flats are selected based on median value
                elif ftype is instrum.flatname:

                    vmin = amin * sat_pt
                    vmax = amax * sat_pt

                    if instrum.is_cam_split(cam_i):

                        # check whether median values are in specified range
                        # bottom side
                        if vmin < m1 < vmax:
                            print('\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}a'.format(cam_i))))
                            af.print_blue("\t* Bottom side selected.")
                            imfits_1 = savefile(fits_id, im1, instrum.get_filter(h, 'C{}a'.format(cam_i)), h)
                            addtodict(dictionary=fits_list_dict,
                                      key=instrum.get_filter(h, 'C{}a'.format(cam_i)),
                                      value='{}{}'.format(workdir, imfits_1))

                        else:
                            if m1 < vmin:
                                af.print_warn("\t* Bottom side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Bottom side rejected:\tSATURATED.")

                        # top side
                        if vmin < m2 < vmax:
                            print('\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}b'.format(cam_i))))
                            af.print_blue("\t* Top side selected.")
                            imfits_2 = savefile(fits_id, im2, instrum.get_filter(h, 'C{}b'.format(cam_i)), h)
                            addtodict(dictionary=fits_list_dict,
                                      key=instrum.get_filter(h, 'C{}b'.format(cam_i)),
                                      value='{}{}'.format(workdir, imfits_2))
                        else:
                            if m2 < vmin:
                                af.print_warn("\t* Top side rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Top side rejected:\tSATURATED.")

                    # Not split frame
                    else:

                        # check whether median value is in specified range
                        if vmin < m < vmax:
                            af.print_blue("\t* Frame selected.")
                            addtodict(dictionary=fits_list_dict,
                                      key=instrum.get_filter(h, 'C{}'.format(cam_i)),
                                      value='{}{}'.format(workdir, fits_fn))

                        else:
                            if m < vmin:
                                af.print_warn("\t* Frame rejected:\tUNDEREXPOSED.")
                            else:
                                af.print_warn("\t* Frame rejected:\tSATURATED.")

            # display image and prompt user
            else:

                if instrum.is_cam_split(cam_i):

                    if (sfrac1 < amin) or (sfrac1 > amax) or (sfrac2 < amin) or (sfrac2 > amax):
                        af.print_warn(
                            "Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in" +
                            " frame.  Skipping frame {}.".format(amin, amax, fits_fn))
                        continue

                    if not noplot:
                        # show top frame
                        ax1 = fig.add_subplot(221)
                        plot_params_calib(ax1, im1, m1, s1, sat_pt, hist=False)

                        # show pixel distribution
                        axhist = fig.add_subplot(222)
                        plot_params_calib(axhist, im1, m1, s1, sat_pt, hist=True)

                        # show bottom frame
                        ax2 = fig.add_subplot(223)
                        plot_params_calib(ax2, im2, m2, s2, sat_pt, hist=False)

                        # show pixel distribution
                        axhist = fig.add_subplot(224)
                        plot_params_calib(axhist, im2, m2, s2, sat_pt, hist=True)
                        fig.subplots_adjust(wspace=0.1, hspace=0.45)

                else:
                    if (sfrac < amin) or (sfrac > amax):
                        af.print_warn(
                            "Warning: median value outside specified range of {:.0%} - {:.0%} of saturation value in" +
                            " frame.  Skipping frame {}.".format(amin, amax, fits_fn)
                        )
                        continue

                    if not noplot:
                        # show frame
                        ax = fig.add_subplot(121)
                        plot_params_calib(ax, im1, m, s, sat_pt, hist=False)

                        # show pixel distribution
                        axhist = fig.add_subplot(122)
                        plot_params_calib(axhist, im1, m, s, sat_pt, hist=True)

                if not noplot:
                    fig.canvas.draw()

                # query user until valid response is provided
                valid_entry = False
                while not valid_entry:
                    if yes:
                        user = 'y'  # TODO: determine if we want this automated with other inputs

                    else:
                        user = input("\nType Y for YES, N for NO, Q for QUIT: ")

                    if user.lower() == 'y':

                        if instrum.is_cam_split(cam_i):
                            imfits_1 = savefile(fits_id, im1,
                                                instrum.get_filter(h, 'C{}a'.format(cam_i)), h)
                            addtodict(dictionary=fits_list_dict,
                                      key=instrum.get_filter(h, 'C{}a'.format(cam_i)),
                                      value='{}{}'.format(workdir, imfits_1))

                            imfits_2 = savefile(fits_id, im2,
                                                instrum.get_filter(h, 'C{}b'.format(cam_i)), h)
                            addtodict(dictionary=fits_list_dict,
                                      key=instrum.get_filter(h, 'C{}b'.format(cam_i)),
                                      value='{}{}'.format(workdir, imfits_2))

                        else:
                            if ftype is instrum.flatname:
                                fl_key = instrum.get_filter(h, 'C{}'.format(cam_i))
                            else:
                                fl_key = 'C{}'.format(cam_i)
                            addtodict(dictionary=fits_list_dict, key=fl_key, value='{}{}'.format(workdir, fits_fn))

                        valid_entry = True

                    elif user.lower() == 'q':  # exit function
                        af.print_bold("Exiting...")
                        os.chdir(start_dir)  # move back to starting directory
                        pl.close('all')  # close image to free memory
                        return

                    elif user.lower() == 'n':  # 'N' selected, skip
                        valid_entry = True

                    else:  # invalid case
                        af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto:
                fig.clear()  # clear image
            hdulist.close()  # close FITs file

    if auto:
        if not noplot:
            af.print_head("\nDisplaying automatically selected {} frames:".format(ftype))
            af.show_list(fits_list_dict)
    else:
        if not noplot:
            pl.close('all')  # close image to free memory

    if save_select:
        dt = datetime.datetime.now()
        fnout = '{}_'.format(ftype) + dt.isoformat().split('.')[0].replace('-', '').replace(':', '') + '.json'
        # python json extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        json_helper.save_dict_to_json(fits_list_dict, fnout)  # save dictionary to json

    os.chdir(start_dir)  # move back to starting directory

    return fits_list_dict


def image_summary(im, sat_pt, cam_i, instrum, split=False):
    """
    NAME:
        image_summary
    PURPOSE:
        Calculate median, robust scatter, and fraction of saturation point for slice.  If
        split array then will output information for both sides of array
    INPUTS:
        im      - data
        sat_pt  - saturation point
        cam_i   - camera that is being used
        instrum - module that contains instrument specific information about cameras and split
        split   - boolean that tells if camera has split filters
    """

    if split:
        im1 = im[instrum.slice('C{}a'.format(cam_i))]
        m1 = np.median(im1)
        s1 = af.robust_sigma(im1)
        sfrac1 = float(m1) / sat_pt
        im2 = im[instrum.slice('C{}b'.format(cam_i))]
        m2 = np.median(im2)
        s2 = af.robust_sigma(im2)
        sfrac2 = float(m2) / sat_pt
        print('\t* Median of left side is {} counts ({:.0%} of saturation level).'.format(m1, sfrac1))
        print('\t* Median of right side is {} counts ({:.0%} of saturation level).'.format(m2, sfrac2))

        return [[im1, im2], [m1, m2], [s1, s2], [sfrac1, sfrac2]]

    else:
        im1 = im[instrum.slice('C{}'.format(cam_i))]
        m = np.median(im1)
        print("Median: " + str(m))
        s = af.robust_sigma(im1)
        sfrac = float(m) / sat_pt
        print('\t* Median is {} counts ({:.0%} of saturation level).'.format(m, sfrac))

        return [im1, m, s, sfrac]


def addtodict(dictionary=None, key=None, value=None):
    """
    NAME:
        addtodict
    PURPOSE:
        Adds (key, value) to dictionary by initializing or appending
    INPUTS:
        dictionary  - dictionary to add to
        key   - key to add to dictionary
        value - value of key to add to dictionary
    EXAMPLE:
        addtodict(dictionary=dictionary, key='test', value='testing')
    """

    try:
        dictionary[key].append(value)
    except:
        dictionary[key] = [value]


def savefile(file, im, band_filter, h):
    """
    NAME:
        savefile
    PURPOSE:
        Adds band_filter keyword and saves file (usually for split band_filter cameras)
    INPUT:
        file   - filename without extension
        im     - data
        band_filter - band_filter name to add
        h      - header
    EXAMPLE:
        savefile('test', data, 'Z', header)
    """
    newfile = '{}_{}.fits'.format(file, band_filter)
    h['FILTER'] = band_filter
    if os.path.exists(newfile):
        os.remove(newfile)  # delete old copy
    pf.writeto(newfile, im, header=h, overwrite=True)  # save object frame
    return newfile


def plot_params_calib(ax, im, m, s, sat_pt, hist=False):
    """
    NAME:
        plot_params_calib
    PURPOSE:
        Plots calibration files (image and histogram) for user selection
    INPUTS:
        ax     - plot reference that we will be using to plot images
        im     - data to plot
        m      - median
        s      - robust scatter
        sat_pt - saturation point
    EXAMPLE:
        plot_params_calib(ax, im, m, s, sat_pt, hist=False)
    NOTE:
        Will not plot unless you have a show() or something to display
    """

    z1, z2 = af.zscale(im)
    if z2 <= z1:
        z1 = m - s
        z2 = m + s

    if hist:
        ax.hist(im.flatten(), bins=50, normed=True, log=True, range=(0, sat_pt))
        ax.set_xlim((0, sat_pt))
        ax.set_xticks([0, 0.5 * sat_pt, sat_pt])
        ax.set_xticklabels(['0%', '50%', '100%'])
        ax.set_yticks([])
        ax.grid()
        ax.set_title("Pixel distribution")
        ax.set_ylabel("Log Scale")
    else:
        ax.imshow(im, vmin=z1, vmax=z2, origin='lower', cmap=pl.cm.gray, interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(r"Median = {}, $\sigma$ = {:.1f}".format(int(m), s))
        ax.set_xlabel(r"Median is {:.0%} of saturation level.".format(float(m) / sat_pt))


def plot_params_science(ax, disp_im, band_filter, h, central=False, window_zoom=4):
    """
    NAME:
        plot_params_science
    PURPOSE:
        Plots science files (image and histogram) for user selection
    INPUTS:
        ax          - plot reference that we will be using to plot images
        disp_im     - data to plot
        band_filter      - band_filter for labeling plot
        h           - header
        central     - boolean to center image to zoomed window
        window_zoom - zoom level for closer look
    NOTE:
        Will not plot unless you have a show() or something to display
    """

    z1, z2 = af.zscale(disp_im)
    ax.set_xticks([])
    ax.set_yticks([])
    xm, ym = np.array(disp_im.shape, dtype=float) / 2
    xr = xm / float(window_zoom)
    yr = ym / float(window_zoom)

    if central:
        ax.imshow(disp_im[xm - xr:xm + xr, ym - yr:ym + yr], vmin=z1, vmax=z2, origin='lower',
                  cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im[xm - xr:xm + xr, ym - yr:ym + yr], levels=[z2],
                   origin='lower', colors='r')
        ax.set_title("Zoomed region")

    else:
        ax.imshow(disp_im, vmin=z1, vmax=z2, origin='lower',
                  cmap=pl.cm.gray, interpolation='none')
        ax.contour(disp_im, levels=[z2], origin='lower', colors='r')
        ax.add_patch(Rectangle((ym - yr, xm - xr), 2 * yr, 2 * xr, ec='b', fc='none', lw=2))
        ax.set_title(r"{} band".format(band_filter))


def choose_science(instrument, workdir='.', targetdir='.', cams=(0, 1, 2, 3), auto=False, save_select=True,
                   figsize=(10, 10), window_zoom=4, calibrate=False, noplot=False, yes=False):
    """
    PURPOSE:
        Display science images for verification by user
    INPUT:
        instrument  - instrument name defined in instrument_dict (ex. 'ratir')
        workdir     - directory where function is to be executed
        targetdir   - directory where selected frames and lists are output
        cams        - camera numbers to process data, all by default
        auto        - select all science frames
        save_select - save dictionary of selected frames to python json file
        figsize     - dimensions of figure used to display frames for selection
        window_zoom - zoom level for closer look
    EXAMPLE:
        file_dict = choose_science('ratir', workdir = 'path/to/data/',
            targetdir = 'path/to/processeddata/',cams = [#,#,...], calibrate=True)
    """

    instrum = instrument_dict[instrument]

    # check for non-list camera argument
    if not isinstance(cams, Iterable):
        cams = [cams]  # convert to list

    # check for non-integer camera designators
    if type(cams[0]) is not int:
        af.print_err("Error: cameras must be specified by an integer. Exiting...")
        return

    pl.ion()  # pylab in interactive mode

    # move to working directory
    start_dir = os.getcwd()
    os.chdir(workdir)

    d = os.getcwd().split('/')[-1]  # name of current directory
    if not auto:
        af.print_head("\nDisplaying science frames in {} for selection:".format(d))
    else:
        af.print_head("\nSelecting all science frames in {}:".format(d))

    # dictionary to store selected fits files by camera or band_filter
    fits_list_dict = {}

    # remove tailing / from target directory name if present
    if targetdir[-1] == '/':
        targetdir = targetdir[:-1]

    # make target directory if it does not exist
    if not os.path.exists(targetdir):
        af.print_blue("Creating target directory: {}".format(targetdir))
        os.makedirs(targetdir)
    # warn user if previous files may be overwritten
    else:
        if not yes:
            af.print_warn("Warning: Target directory exists. Existing files will be overwritten.")
            resp = input("Proceed? (y/n): ")
            if resp.lower() != 'y':
                af.print_bold("Exiting...")
                os.chdir(start_dir)  # move back to starting directory
                return
        shutil.rmtree(targetdir)
        os.makedirs(targetdir)

    fits_check = glob('????????T??????C??.fits')

    if len(fits_check) == 0:
        instrum.change_file_names(glob(instrum.original_file_format()))

        fits_check_2 = glob('????????T??????C??.fits')

        if len(fits_check_2) == 0:
            af.print_err("Error: no files with correct format after file name changes")
            return

    # open figure for images if not auto
    if not auto and not noplot:
        fig = pl.figure(figsize=figsize)

    # work on FITs files for specified cameras
    for cam_i in cams:

        # print current camera number
        af.print_under("\n{:^50}".format('CAMERA {}'.format(cam_i)))

        if calibrate:
            # get master dark and bias frames for current camera if required
            if instrum.has_cam_bias(cam_i):
                # if instrum.has_cam_bias(cam_i) == True:
                mbias_fn = '{}_C{}.fits'.format(instrum.biasname, cam_i)

                if not os.path.exists(mbias_fn):
                    af.print_err(
                        'Error: {} not found.  Move master bias file to working directory to proceed.'.format(mbias_fn))
                    continue
                # else:
                mbias_data = pf.getdata(mbias_fn)

            if instrum.has_cam_dark(cam_i):
                mdark_fn = '{}_C{}.fits'.format(instrum.darkname, cam_i)

                if not os.path.exists(mdark_fn):
                    af.print_err(
                        'Error: {} not found.  Move master dark file to working directory to proceed.'.format(mdark_fn))
                    continue
                # else:
                mdark_data = pf.getdata(mdark_fn)

        # find raw files of selected type for this camera
        fits_list = glob('????????T??????C{}{}.fits'.format(cam_i, instrum.ftype_post[instrum.objname]))
        if len(fits_list) == 0:
            af.print_warn("Warning: no fits files found.  Skipping camera {}.".format(cam_i))
            continue

        # look at FITs data sequentially
        for fits_fn in fits_list:

            fits_id = fits_fn.split('.')[0]  # fits file name with extension removed
            print('{}'.format(fits_fn))

            # open data
            hdulist = pf.open(fits_fn)
            im = hdulist[0].data
            h = hdulist[0].header

            if calibrate:
                # get master flat frame for current band_filter
                if instrum.is_cam_split(cam_i):
                    mflat_fn1 = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h, 'C{}a'.format(cam_i)))
                    if not os.path.exists(mflat_fn1):
                        af.print_err(
                            'Error: {} not found.  Move master flat file to working directory to proceed.'.format(
                                mflat_fn1))
                        continue
                    else:
                        mflat_data1 = pf.getdata(mflat_fn1)
                    mflat_fn2 = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h, 'C{}b'.format(cam_i)))
                    if not os.path.exists(mflat_fn2):
                        af.print_err(
                            'Error: {} not found.  Move master flat file to working directory to proceed.'.format(
                                mflat_fn2))
                        continue
                    else:
                        mflat_data2 = pf.getdata(mflat_fn2)

                else:
                    mflat_fn = '{}_{}.fits'.format(instrum.flatname, instrum.get_filter(h, 'C{}'.format(cam_i)))
                    if not os.path.exists(mflat_fn):
                        af.print_err(
                            'Error: {} not found.  Move master flat file to working directory to proceed.'.format(
                                mflat_fn))
                        continue
                    else:
                        mflat_data = pf.getdata(mflat_fn)

            # get image statistics
            if instrum.is_cam_split(cam_i):
                im1 = im[instrum.slice('C{}a'.format(cam_i))]
                im2 = im[instrum.slice('C{}b'.format(cam_i))]
            else:
                im1 = im2 = im[instrum.slice('C{}'.format(cam_i))]

            # display image and prompt user
            if not auto:

                if instrum.is_cam_split(cam_i):

                    disp_im1 = np.copy(im1)
                    disp_im2 = np.copy(im2)

                    if calibrate:
                        if instrum.has_cam_bias(cam_i):
                            disp_im1 -= mbias_data
                            disp_im2 -= mbias_data

                        if instrum.has_cam_dark(cam_i):
                            disp_im1 -= mbias_data * instrum.get_exptime(h)
                            disp_im2 -= mbias_data * instrum.get_exptime(h)

                        disp_im1 = np.divide(disp_im1, mflat_data1)
                        disp_im2 = np.divide(disp_im2, mflat_data2)

                    if not noplot:
                        # display top
                        ax1 = fig.add_subplot(221)
                        plot_params_science(ax1, disp_im1, instrum.get_filter(h, 'C{}a'.format(cam_i)),
                                            h, central=False)

                        # and central subregion
                        ax1s = fig.add_subplot(222)
                        plot_params_science(ax1s, disp_im1, instrum.get_filter(h, 'C{}a'.format(cam_i)),
                                            h, central=True)

                        # display bottom
                        ax2 = fig.add_subplot(223)
                        plot_params_science(ax2, disp_im2, instrum.get_filter(h, 'C{}b'.format(cam_i)),
                                            h, central=False)

                        # and central subregion
                        ax2s = fig.add_subplot(224)
                        plot_params_science(ax2s, disp_im2, instrum.get_filter(h, 'C{}b'.format(cam_i)),
                                            h, central=True)

                else:
                    disp_im1 = np.copy(im1)

                    if calibrate:
                        if instrum.has_cam_bias(cam_i):
                            disp_im1 -= mbias_data

                        if instrum.has_cam_dark(cam_i):
                            disp_im1 -= mdark_data * instrum.get_exptime(h)

                        disp_im1 = np.divide(disp_im1, mflat_data)

                    if not noplot:
                        ax = fig.add_subplot(121)
                        plot_params_science(ax, disp_im1, instrum.get_filter(h, 'C{}'.format(cam_i)),
                                            h, central=False)

                        # and central subregion
                        axs = fig.add_subplot(122)
                        plot_params_science(axs, disp_im1, instrum.get_filter(h, 'C{}'.format(cam_i)),
                                            h, central=True)

                if not noplot:
                    fig.set_tight_layout(True)
                    fig.canvas.draw()

            if instrum.is_cam_split(cam_i):
                if instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}a'.format(cam_i))) != 0:
                    print("\t* The target is focused on the {} band_filter.".format(
                        instrum.get_filter(h, 'C{}a'.format(cam_i))))
                elif instrum.get_centered_filter(h, cam_i).count(instrum.get_filter(h, 'C{}b'.format(cam_i))) != 0:
                    print("\t* The target is focused on the {} band_filter.".format(
                        instrum.get_filter(h, 'C{}b'.format(cam_i))))
                else:
                    af.print_warn(
                        ("\t* Warning: The target is NOT focused on an split band_filter. The target is focused on " +
                         "the {} band_filter.").format(instrum.get_centered_filter(h, cam_i)))
            else:
                # print band_filter name
                print('\t* Filter used: {}'.format(instrum.get_filter(h, 'C{}'.format(cam_i))))

            # query user until valid response is provided
            valid_entry = False
            while not valid_entry:

                # either select all if auto, or have user select
                if auto or yes:
                    user = 'y'
                else:
                    user = input("\nType Y for YES, N for NO, Q for QUIT: ")

                if user.lower() == 'y' and instrum.is_cam_split(cam_i):
                    # filterA = instrum.get_filter(h, 'C{}a'.format(cam_i)).lower()
                    filterA = instrum.get_filter(h, 'C{}a'.format(cam_i))
                    # filterB = instrum.get_filter(h, 'C{}b'.format(cam_i)).lower()
                    filterB = instrum.get_filter(h, 'C{}b'.format(cam_i))

                    if instrum.get_centered_filter(h, cam_i).count(filterB) != 0:
                        direction = 't'
                    elif instrum.get_centered_filter(h, cam_i).count(filterA) != 0:
                        direction = 'b'
                    elif instrum.get_centered_filter(h, cam_i).lower() == 'r' and (filterB == 'h' or filterA == 'z'):
                        # keeping frames 'r' centered with split band_filter 'Z/Y' or 'J/H'
                        direction = 'b'
                    else:
                        af.print_warn("\t* Warning: Skipping frame not centered on split band_filter.")
                        user = 'n'
                        direction = ''

                if user.lower() == 'y':

                    h_c = h.copy()

                    if instrum.is_cam_split(cam_i):

                        if direction.lower() == 'b':
                            f_img_post = 'a'
                            f_sky_post = 'b'
                        else:
                            f_img_post = 'b'
                            f_sky_post = 'a'

                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, instrum.objname,
                                                           instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_img_post)))
                        im_img = im[instrum.slice('C{}{}'.format(cam_i, f_img_post))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}{}'.format(cam_i, f_img_post))
                        pf.writeto(imfits, im_img, header=hnew, overwrite=True)  # save object frame
                        # 'has_key()' depreciated; changed  to 'in'
                        if instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_img_post)) in fits_list_dict:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_img_post))].append(imfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_img_post))] = [imfits]

                        # band_filter side with sky, now saved as object, but different list to keep track
                        skyfits = '{}/{}_{}_{}.fits'.format(
                            targetdir, fits_id, instrum.objname,
                            instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_sky_post))
                        )
                        im_sky = im[instrum.slice('C{}{}'.format(cam_i, f_sky_post))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}{}'.format(cam_i, f_sky_post))
                        pf.writeto(skyfits, im_sky, header=hnew, overwrite=True)  # save sky frame
                        # 'has_key()' depreciated; changed  to 'in'
                        if instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_sky_post)) in fits_list_dict:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_sky_post))].append(skyfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}{}'.format(cam_i, f_sky_post))] = [skyfits]

                        valid_entry = True

                    else:

                        imfits = '{}/{}_{}_{}.fits'.format(targetdir, fits_id, instrum.objname, cam_i)
                        im_img = im[instrum.slice('C{}'.format(cam_i))]
                        hnew = instrum.change_header_keywords(h_c, 'C{}'.format(cam_i))
                        pf.writeto(imfits, im_img, header=hnew, overwrite=True)
                        # 'has_key()' depreciated; changed  to 'in'
                        if instrum.get_filter(h_c, 'C{}'.format(cam_i)) in fits_list_dict:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}'.format(cam_i))].append(imfits)
                        else:
                            fits_list_dict[instrum.get_filter(h_c, 'C{}'.format(cam_i))] = [imfits]

                        valid_entry = True

                elif user.lower() == 'q':  # exit function
                    af.print_bold("Exiting...")
                    os.chdir(start_dir)  # move back to starting directory
                    pl.close('all')  # close image to free memory
                    return

                elif user.lower() == 'n':  # 'N' selected, skip
                    valid_entry = True

                else:  # invalid case
                    af.print_warn("'{}' is not a valid entry.".format(user))

            if not auto and not noplot:
                fig.clear()  # clear image
            hdulist.close()  # close FITs file

    if auto:
        if not noplot:
            af.print_head("\nDisplaying automatically selected science frames:")
            af.show_list(fits_list_dict)
    else:
        if not noplot:
            pl.close('all')  # close image to free memory

    if save_select:
        dt = datetime.datetime.now()
        fnout = 'object_' + dt.isoformat().split('.')[0].replace('-', '').replace(':',
                                                                                  '') + '.json'  # python json extension
        af.print_head("\nSaving selection dictionary to {}".format(fnout))
        json_helper.save_dict_to_json(fits_list_dict, fnout)  # save dictionary to json

    os.chdir(start_dir)  # move back to starting directory

    return fits_list_dict
