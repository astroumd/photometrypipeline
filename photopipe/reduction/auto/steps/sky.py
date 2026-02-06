import glob
import os
import re
import datetime
import time

import astropy.io.fits as pf
import numpy as np
import scipy
from scipy import ndimage

import photopipe.reduction.auto.steps.autoproc_depend as apd


inpipevar = {
    'autoastrocommand': 'autoastrometry', 'getsedcommand': 'get_SEDs', 'sexcommand': 'sex', 'swarpcommand': 'swarp',
    'rmifiles': 0, 'prefix': '', 'datadir': '', 'imworkingdir': '', 'overwrite': 0, 'verbose': 1, 'flatfail': '',
    'fullastrofail': '',	'pipeautopath': '', 'refdatapath': '', 'defaultspath': ''
}


def median_filter_masking(image, size=50):
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('start_nan_percentage: {}'.format(nan_percent))
    median_filter_image = ndimage.median_filter(image, size=size)
    image[np.isnan(image)] = median_filter_image[np.isnan(image)]
    nan_percent = 100 * np.count_nonzero(np.isnan(image)) / (image.shape[0] * image.shape[1])
    print('end_nan_percentage: {}'.format(nan_percent))
    return image


def autopipeskysubmed(pipevar=None):
    """
    NAME:
        autopipeskysubmed
    PURPOSE:
        Subtract median, does NOT use master sky
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.skypipeproc
    """

    print('SKY-SUBTRACT MEDIAN ONLY')
    if pipevar is None:
        pipevar = inpipevar
    # Find data that needs to be sky subtracted
    files = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    # For each file if output files don't exist or override set check if we have master
    # skyflat for band_filter, sky subtract if it exists using skypipeproc
    for f in files:
        print(f)
        fileroot = os.path.basename(f)
        outfile = pipevar['imworkingdir'] + 's' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping sky subtraction for ' + f + '. File already exists')
            continue

        head = pf.getheader(f)
        data = pf.getdata(f)
        good = np.isfinite(data)
        bad = ~good  # Opposite of boolean array good

        data -= np.median(data[~np.isnan(data)])
        data = np.nan_to_num(data)
        # x = np.arange(0, data.shape[1])
        # y = np.arange(0, data.shape[0])
        # data = np.ma.masked_invalid(data)
        # xx, yy = np.meshgrid(x, y)
        # x1 = xx[~data.mask]
        # y1 = yy[~data.mask]
        # newarr = data[~data.mask]
        # data = interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy), method='nearest')

        # head_filter = head['FILTER'] TODO: determine if this not being used is a bug

        if pipevar['verbose'] > 0:
            print('Sky subtracting median only')
            print("# of bad pix: {}".format(np.sum(bad)))

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeprocmed ' + date)

        pf.writeto(outfile, data, head, overwrite=pipevar['overwrite'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')


def autopipemakesky(pipevar=inpipevar, by_targ=False):
    """
    NAME:
        autopipemakesky
    PURPOSE:
        Combine sky flats based on band_filter type (sigma clipping for sources)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipemakesky(pipevar=inpipevar)
    DEPENDENCIES:
        astroproc_depend.skypipecombine_new, astroproc_depend.medclip, Sextractor
    FUTURE IMPROVEMENTS:
        skypipecombine slow, find better algorithm
    """

    print('MAKE SKY')

    # Copies necessary parameter file for sextractor if not in current working directory
    if not os.path.isfile('source.param'):
        os.system('cp ' + pipevar['defaultspath'] + '/source.param .')
    if not os.path.isfile('sex_source.config'):
        os.system('cp ' + pipevar['defaultspath'] + '/sex_source.config .')
    if not os.path.isfile('sex.conv'):
        os.system('cp ' + pipevar['defaultspath'] + '/sex.conv .')
    if not os.path.isfile('defaulf.nnw'):
        os.system('cp ' + pipevar['defaultspath'] + '/default.nnw .')

    # Finds files with given prefix
    files = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    filters = []
    targets = []
    for file in files:
        head = pf.getheader(file)
        band_filter = head['FILTER']
        target = re.sub(r'\s+', '', head['TARGNAME'])
        filters += [band_filter]
        targets += [target]
    filters = np.array(filters)
    targets = np.array(targets)

    # Unique list of filters
    filterlist = set(filters)
    files = np.array(files)

    # For each unique target and band_filter combination, combine sky files using skycombine if more than 2 files
    # Otherwise return list of unprocessed files
    for filt in filterlist:

        thisfilter = np.where(filters == filt)
        if len(thisfilter) == 0:
            continue

        thisfiltertargs = set(targets[thisfilter])

        if not by_targ:
            skyflats = np.where(filters == filt)
            outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '.fits'
            call_skycombine(files, skyflats, outflatname, filt, pipevar)

        else:
            for targ in thisfiltertargs:
                skyflats = np.where(np.logical_and(filters == filt, targets == targ))
                outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '_' + targ + '.fits'
                call_skycombine(files, skyflats, outflatname, filt, pipevar)

    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')


def call_skycombine(files, skyflats, outflatname, filt, pipevar):
    if len(skyflats[0]) >= 2:

        if os.path.isfile(outflatname) and pipevar['overwrite'] == 0:
            print('Skipping makesky for ' + filt + '. File already exists')
        else:
            if pipevar['verbose']:
                print(filt, '-band sky flats.')
                print(files[skyflats])

            skypipecombine(files[skyflats], outflatname, filt, pipevar, removeobjects=True, type='sky')
    else:
        print('Unable to produce a flat field for this setting: ' + filt)
        print(
            'Will not be able to further process ' + str(len(skyflats)) +
            ' image(s) without a flat from another source:'
        )


def autopipeskysub(pipevar=inpipevar, by_targ=False):
    """
    NAME:
        autopipeskysub_targets
    PURPOSE:
        Subtracts master sky flat from data and subtracts median.
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        skypipeproc
    """

    print('SKY-SUBTRACT')

    # Find data that needs to be sky subtracted
    files = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    skys = glob.glob(pipevar['imworkingdir'] + '*sky-*.fits')

    if len(skys) == 0:
        print('No master sky files found, cannot sky subtract')
        return

    # Find the associated band_filter of each master skyflat
    skyset = []
    skyfilts = []
    for sky in skys:
        head = pf.getheader(sky)
        band_filter = head['FILTER']
        target = [re.sub(r'\s+', '', head['TARGNAME'])]
        skyfilts += [band_filter]
        skyset += [(band_filter, target)]

    # For each file if output files don't exist or override set check if we have master
    # skyflat for band_filter, sky subtract if it exists using skypipeproc
    for file in files:

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping sky subtraction for ' + file + '. File already exists')
            continue

        head = pf.getheader(file)
        band_filter = head['FILTER']
        target = [re.sub(r'\s+', '', head['TARGNAME'])]

        # Find corresponding master skyflat
        try:
            if by_targ:
                skyloc = skyset.index((band_filter, target))
            else:
                skyloc = skyfilts.index(band_filter)
            skyfile = skys[skyloc]
        except:
            print('Sky field not found for ', file)
            pipevar['flatfail'] += ' ' + file
            continue
        if pipevar['verbose'] > 0:
            print('Sky subtracting', file, 'using', skyfile)

        skypipeproc(file, skyfile, outfile)

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')


def skypipecombine(
        filelist, outfile, filt, pipevar, removeobjects=None,
        objthresh=2.0, objthresh2=2.0, algorithm='median', trimlo=None, trimhi=None,
        mincounts=1, maxcounts=55000, satlevel=30000, type=None
):
    """
    NAME:
        skypipeindividual
    PURPOSE:
        Create sigma clipped sky flat.  Scales each file based on the overall
        sigma clipped median, then removes objects selected with sextractor (uses flux
        fraction radius) in each file. Removes saturated pixels.  Calculates the sky
        values of each pixel and saves anything with non-finite values (saturated or
        source) to the median of the entire frame.  Then subtract the sky from each
        individual frame


    INPUTS:
        filelist - files to be processed
        filt	 - band_filter of files
        pipevar  - pipeline parameters in dictionary
    OPTIONAL KEYWORDS:
        removeobjects 	- specifies if you want objects removed
        objthresh	- sets sigma in removeobjects (default is 2)
        objthresh2      - set sigma for the refined selection of the object in the sky imae (default is 0.5)
        algorithm       - algorithm to solve (mean or median, default is median)
        trimlo			- trim off bottom of data in mean algorithm mode (default is 25%)
        trimhi			- rim off top of data in mean algorithm mode (default is 25%)
        mincounts		- sets minimum counts allowed (default is 1)
        maxcounts		- sets maximum counts allowed (default is 55000)
        satlevel		- sets saturation level (default is 30000)
        type			- sets 'SKYTYPE' keyword in header of outfile to this string
    EXAMPLE:
        skypipecombine(filelist, 'sky-filt.fits', removeobjects=True, type='sky')
    DEPENDENCIES:
        medclip, Sextractor
    FUTURE IMPROVEMENTS:
        medclip slow find faster solution
        Need to take saturation level from header?
        Saved header is from middle file, maybe use blank?
    """

    # source position in pixel
    posx = 100
    posy = 100
    half_src_pixels = 20

    # Sets defaults for trimming (25% of list)
    if trimlo is not None:
        (len(filelist) + 1) / 4
    if trimhi is not None:
        trimlo

    # If given list, then grab all filenames, saved to files
    # if len(filelist) == 1:
    #    f = open(filelist,'r')
    #    files = f.read().splitlines()
    #    f.close()
    # else:
    #    files = filelist

    files = filelist
    nfiles = len(files)
    nmid = int(len(files) / 2)
    print(files[int(nmid)])

    # Read in middle file and initialize arrays
    f = pf.open(files[nmid])
    data_m = f[0].data
    head_m = f[0].header
    f.close()

    nx = head_m['NAXIS1']
    ny = head_m['NAXIS2']
    pixsc = head_m['PIXSCALE']

    data = np.zeros((nfiles, ny, nx)) + float('NaN')
    skymeds = []
    usefiles = []

    z = 0
    # For each file and make sure size matches middle file, calculate sigma clipped
    # median (3sig, 6 iter), then if within counts limit save data into 3d data cube
    # and save clipped median into skymed, and mark file as usable
    # Increment z by one when this is true
    for file in files:
        f = pf.open(file)
        data_i = f[0].data
        head_i = f[0].header
        f.close()

        good = np.isfinite(data_i)
        bad = ~good  # Opposite of boolean array good
        # print("Bad pix: {}".format(np.sum(bad)))

        inx = head_i['NAXIS1']
        iny = head_i['NAXIS2']

        if (inx != nx) or (iny != ny):
            print('File ' + file + ' has wrong dimensions (' + str(inx) +
                  ' x ' + str(iny) + '; should have ' + str(nx) + ' x ' + str(ny) + ')')

        # Perform 3 sigma clipped median and save to inmeds
        inmed, instd = medclip(data_i, clipsig=5, maxiter=3)

        # If median is within limits save data, otherwise exclude files
        if mincounts <= inmed <= maxcounts:
            if pipevar['verbose'] > 0:
                print(file + ' (' + str(inmed) + ' counts/pix)')

            skymeds += [inmed]
            usefiles += [file]
            data[z, :, :] = data_i
            z += 1
        else:
            if inmed < mincounts:
                print(file + ' (' + str(inmed) + ' counts/pix) - too few counts; excluding')
            if inmed > maxcounts:
                print(file + ' (' + str(inmed) + ' counts/pix) - too many counts; excluding')

                # if z < 2:
    #    print 'ERROR - Not enough counts to make a flat with these data!'
    #    return

    # Median of sigma clipped medians
    medsky = np.nanmedian(skymeds)

    # Scale each file by median of sigma clipped medians divided by median of data
    # Corrects for each flat's changing sky background
    for f in np.arange(z):
        # for f in np.arange(z-1):
        factor = medsky / skymeds[f]
        data[f, :, :] = data[f, :, :] * factor

    # Removes extraneous indexes in data for skipped files
    if z != nfiles:
        data = data[0:z, :, :]

    # Removes objects from field by calculating iterative median sigma clipping
    # (5 sigma, 5 iter) and using the calculated stddev to remove 2sigma (or non-default
    # object threshold) data from the median along with values above the saturation limit.

    if removeobjects is not None:
        if pipevar['verbose'] > 0:
            print('  Identifying objects...')

        for f in np.arange(z):
            # for f in np.arange(z-1):

            oridata = np.copy(data[f, :, :])
            indata = np.copy(data[f, :, :])
            meddata = np.copy(data[f, :, :])
            stddata = np.copy(data[f, :, :])
            meddata_src = np.copy(data[f, :, :])
            stddata_src = np.copy(data[f, :, :])

            # Set sources above objthresh  limit to NaN

            datamed, datastd = medclip(indata, clipsig=5, maxiter=5)

            src_piece = np.copy(indata[posy - half_src_pixels - 1:posy + half_src_pixels - 1,
                                posx - half_src_pixels - 1:posx + half_src_pixels - 1])
            meddata_src[:, :] = datamed
            stddata_src[:, :] = datastd
            stddata_src[
                posy - half_src_pixels - 1:posy + half_src_pixels - 1,
                posx - half_src_pixels - 1:posx + half_src_pixels - 1
            ] = 0.001  # take out of the image everything above the median sky level around the source position

            sourcepixels = np.where((indata - meddata_src) >= objthresh * stddata_src)
            # sourcepixels = np.where(abs(indata-datamed)>= objthresh*datastd)
            satpixels = np.where(indata >= satlevel)

            if len(sourcepixels[0]) > 0:
                indata[sourcepixels] = float('NaN')

            if len(satpixels[0]) > 0:
                indata[satpixels] = float('NaN')

            # # Keep sources as NaN Value
            if pipevar['verbose'] > 0:
                good = np.isfinite(indata)
                bad = ~good  # Opposite of boolean array good
                bad_ratio = np.sum(bad)/(bad.shape[0]*bad.shape[1])
                print("Bad pix ratio: {}".format(bad_ratio))

            # compute mean and std deviation in each part of the matrix considering a grid of nsq_piece*2 pixels squares
            nsq_piece = round(14 * (0.36 / pixsc))  # half square length (optimal 5 for pixscle 0.36)

            mi, mj = oridata.shape
            for mii in ((np.arange(mi / (nsq_piece * 2)) + 1) * 2 * nsq_piece) - nsq_piece:
                for mji in ((np.arange(mj / (nsq_piece * 2)) + 1) * 2 * nsq_piece) - nsq_piece:
                    # sq_piece = np.copy(
                    #     indata[
                    #         max(0, mii-nsq_piece):min(mii+nsq_piece+1, mi-1),
                    #         max(0, mji-nsq_piece):min(mji+nsq_piece + 1, mj-1)]
                    # )
                    sq_piece = np.copy(
                        indata[
                            int(mii - nsq_piece):int(mii + nsq_piece + 1), int(mji - nsq_piece):int(mji + nsq_piece + 1)
                        ]
                    )
                    datamedi, datastdi = medclip(sq_piece, clipsig=5, maxiter=5)
                    meddata[int(mii - nsq_piece):int(mii + nsq_piece + 1), int(mji - nsq_piece):int(mji + nsq_piece + 1)] = datamedi
                    stddata[int(mii - nsq_piece):int(mii + nsq_piece + 1), int(mji - nsq_piece):int(mji + nsq_piece + 1)] = datastdi

            # remove the sources above objthresh2 limitv in the squares

            sourcepixels_sq = np.where((indata - meddata) >= objthresh2 * stddata)
            # print("Length of bad pix: {}".format(len(sourcepixels_sq[0])))

            if len(sourcepixels_sq[0]) > 0:
                indata[sourcepixels_sq] = float('NaN')

            good = np.isfinite(meddata)
            bad = ~good  # Opposite of boolean array good
            indata[bad] = np.nan
            data[f, :, :] = indata
            if pipevar['debug'] != 0:
                good2 = np.isfinite(indata)
                out = usefiles[f].replace('fp', 'skymask2_fp')
                pf.writeto(out, good2.astype(np.int), head_m, overwrite=True)

    reflat = np.zeros((ny, nx)) + float('NaN')

    # If algorithm set to median, find 3 sigma clipped median of each pixel
    # excluding NaN values (which are eventually set to median)
    if algorithm == 'median':
        if pipevar['verbose'] > 0:
            print('  Median-combining...')

        for y in np.arange(ny):
            vector = data[:z - 1, y, :]
            temp = np.isfinite(vector)
            me, st = medclip2d(vector, clipsig=3, maxiter=5, overaxis=0)
            reflat[y, :] = me

        # reflat=np.nanmedian(data, axis=0)

        # # Keep sources as NaN Value

    # If algorithm set to mean, takes mean of trimmed sorted values. Default is to
    # trim 25% off top and bottom, if not enough good data, set trimming to 0
    if algorithm == 'mean':

        if pipevar['verbose'] > 0:
            print('  Combining via trimmed mean...')

        for y in np.arange(ny):
            for x in np.arange(nx):
                _slice = data[:, y, x]
                good = np.isfinite(_slice)

                cslice = _slice[good]
                ctgood = len(cslice)

                if ctgood == 0:
                    reflat[y, x] = 1

                itrimlo = trimlo
                itrimhi = trimhi

                while ctgood - itrimlo - itrimhi < 1:
                    itrimlo = max(itrimlo - 1, 0)
                    itrimhi = max(itrimhi - 1, 0)

                cslice = np.sort(cslice)
                cslice = cslice[itrimlo:ctgood - itrimhi]
                reflat[y, x] = np.nanmean(cslice)

    # Interpolates sky flat to remove NaN Values in case of source overlap
    if pipevar['verbose']:
        good = np.isfinite(reflat)
        bad = ~good  # Opposite of boolean array good
        bad_count = np.sum(bad)
        print("Local median filling {} NaN Values".format(bad_count))
    if pipevar['debug'] != 0:
        out = outfile.replace('.fits', 'skymask.fits')
        pf.writeto(out, good.astype(np.int), head_m, overwrite=True)

    t1 = time.perf_counter()
    skyflat = median_filter_masking(reflat)
    t2 = time.perf_counter()

    # skyflat = np.copy(reflat)
    # t1 = time.perf_counter()
    # mi, mj = reflat.shape
    # for mii in ((np.arange(mi / (nsq_piece * 2)) + 1) * 2 * nsq_piece) - nsq_piece:
    #     for mji in ((np.arange(mj / (nsq_piece * 2)) + 1) * 2 * nsq_piece) - nsq_piece:
    #         sq_piece = np.copy(reflat[int(mii - nsq_piece):int(mii + nsq_piece + 1), int(mji - nsq_piece):int(mji + nsq_piece + 1)])
    #         interp_flat = masked_interpolation(sq_piece, method='nearest')
    #         skyflat[int(mii - nsq_piece):int(mii + nsq_piece + 1),int(mji - nsq_piece):int(mji + nsq_piece + 1)] = interp_flat
    # t2 = time.perf_counter()

    if pipevar['debug'] != 0:
        print("Interpolation Time: {:0.2f} s".format(t2-t1))

    # Adds header information to signify what files we used
    for f in np.arange(z - 1):
        head_m['SKY' + str(f)] = usefiles[f]

    if type is not None:
        head_m['SKYTYPE'] = type

    date = datetime.datetime.now().isoformat()
    head_m.add_history('Processed by skypipecombine ' + date)

    if pipevar['verbose'] > 0:
        print('  Written to ' + outfile)

    pf.writeto(outfile, skyflat, head_m, overwrite=True)


def skypipeproc(filename, flatname, outfile, flatminval=None, flatmaxval=None):
    """
    NAME:
        skypipeproc
    PURPOSE:
        Subtracts sky flat from data and then subtracts median of that from remaining data.
    INPUTS:
        filename - file or list of files to be sky subtracted
        flatname - sky flat fits file
        outfile  - name of output file
    OPTIONAL KEYWORDS:
        flatminval - minimum required value in flat (default for skycts calculation is 0.1)
        flatmaxval - maximum required value in flat
    EXAMPLE:
        skypipeproc(filename, flatname, outfile)
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

    # Open flat
    flat = pf.getdata(flatname)

    # med = np.nanmedian(flat)

    # For each input file check if same size as flats (required). If there is a minimum
    # or maximum flat value set, forces values outside of that range to NaN. Use finite
    # values above 0.1 to determine skycounts, and subtract flat along with median of
    # flattened data. Saves to new fits file
    for file_ in files:
        f = pf.open(file_)
        data = f[0].data
        head = f[0].header
        f.close()

        if np.shape(data) != np.shape(flat):
            print(file_ + ' could not be flat subtracted because it is not the same' + \
                  ' size as the master flat, remove file to avoid confusion')
            return

        if flatmaxval is not None:
            w = np.where(flat > flatminval)
            if len(w[0]) != 0:
                flat[w] = float('NaN')

        if flatminval is not None:
            w = np.where(flat < flatminval)
            if len(w[0]) != 0:
                flat[w] = float('NaN')
            goodsignal = np.where((flat >= flatminval) & (np.isfinite(flat)))
        else:
            goodsignal = np.where((flat >= 0.1) & (np.isfinite(flat)))

            # Scale skyflat, subtract scaled skyflat, and subtract median of subsequent flat
        # subtracted data. Calculate skycounts from data (above minimum, or
        # by default above 0.1)
        flattmp = np.median(flat[np.isfinite(flat)])
        imgtmp = np.median(data[np.isfinite(data)])

        scalefr = imgtmp / flattmp
        fdata = data - scalefr * flat

        tmp = np.median(fdata[np.isfinite(fdata)])
        fdata = fdata - tmp

        skycts = np.nanmedian(fdata[goodsignal])

        # Adds header keywords to denote new median counts and file we used to flatfield
        head['SFLATFLD'] = flatname
        head['SKYCTS'] = (skycts, 'Sky counts')

        try:
            head['CTRATE'] = (skycts / head['EXPTIME'], 'Sky counts per second')
        except:
            print('No EXPTIME keyword')

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeproc ' + date)

        apd.write_fits(outfile, fdata, head)


def medclip(indata, clipsig=3.0, maxiter=5, verbose=0):
    """
    NAME:
        medclip
    PURPOSE:
        Median iterative sigma-clipping
    INPUTS:
        indata - array to be clipped
    OPTIONAL KEYWORDS:
        clipsig - sigma to clip around
        maxiter - maximum number of times to clip
        verbose - allow to print messages
    EXAMPLE:
        med, sigma = medclip(indata, sigma=5.0)
    """

    # Flatten array
    skpix = indata.reshape(indata.size, )
    keep = np.isfinite(skpix)
    skpix = skpix[keep]
    # print("Skypix: {}".format(len(skpix)))
    # print("Med: {}, Std: {}".format(np.nanmedian(skpix), np.nanstd(skpix)))
    if len(skpix) == 0:  # In case of all NaN slice
        return np.nan, np.nan
    ct = indata.size
    iteration = 0
    numrej = len(skpix)
    ndata = len(skpix)

    while (iteration < maxiter) and (numrej > min(ndata * 0.01, 50)):
        lastct = ct
        medval = np.nanmedian(skpix)
        sig = np.nanstd(skpix)
        wsm = np.where(abs(skpix - medval) <= clipsig * sig)
        ct = len(wsm[0])
        # print("Ct: {}".format(ct))
        if ct > 0:
            skpix = skpix[wsm]

        numrej = abs(ct - lastct)
        if ct <= 2:
            # print("Count case!!!!")
            return np.nan, np.nan
        iteration += 1

    med = np.nanmedian(skpix)
    sigma = np.nanstd(skpix)

    if verbose:
        print('%.1f-sigma clipped median' % clipsig)
        print('Mean computed in %i iterations' % iteration)
        # noinspection PyStringFormat
        print('Mean = %.6f, sigma = %.6f' % (med, sigma))

    return med, sigma


def medclip2d(indata, clipsig=3.0, maxiter=5, verbose=0, overaxis=0):
    """
    NAME:
        medclip2d
    PURPOSE:
        Median iterative sigma-clipping over 2d array
    INPUTS:
        indata - array to be clipped
    OPTIONAL KEYWORDS:
        clipsig - sigma to clip around
        maxiter - maximum number of times to clip
        verbose - allow to print messages
        overaxis - axis that we want to take median over
    EXAMPLE:
        med, sigma = medclip2d(indata, sigma=5.0, overaxis=0)
    """

    # Flatten array
    skpix = np.ma.masked_array(indata)
    skpix = np.ma.masked_invalid(skpix)
    iter = 0

    while (iter < maxiter):
        medval = np.nanmedian(skpix, axis=overaxis)
        sig = np.nanstd(skpix, axis=overaxis)

        if overaxis == 0:
            mask = (abs(skpix - medval) < clipsig * sig)
        else:
            mask = (abs(skpix.T - medval) < clipsig * sig)
        if (mask == skpix.mask).all:
            break
        skpix.mask = mask
        # if ct <=2: return 'Too few remaining'
        iter += 1

    med = np.nanmedian(skpix, axis=overaxis)
    sigma = np.nanstd(skpix, axis=overaxis)

    if verbose:
        print('%.1f-sigma clipped median' % clipsig)
        print('Mean computed in %i iterations' % iter)
        print('Mean = %.6f, sigma = %.6f' % (med, sigma))

    return med, sigma


def masked_interpolation(image, method='nearest'):
    """
    Interpolates over masked locations
    Parameters
    ----------

    image : np.array
    method : ‘linear’, ‘nearest’, ‘cubic’}, optional

    Returns np.array
    -------
    """
    bad_pixel_mask = ~np.isfinite(image)
    # print("Bad Pixel count {}".format(np.sum(bad_pixel_mask)))
    x = np.arange(0, image.shape[1])
    y = np.arange(0, image.shape[0])
    # interpolated_image = np.copy(image)
    # interpolated_image[bad_pixel_mask] = np.nan
    assert isinstance(image, np.ndarray)
    interpolated_image = np.ma.masked_invalid(image.copy())
    xx, yy = np.meshgrid(x, y)
    x1 = xx[~interpolated_image.mask]
    y1 = yy[~interpolated_image.mask]
    newarr = interpolated_image[~interpolated_image.mask]
    assert isinstance(x1, np.ndarray)
    assert isinstance(y1, np.ndarray)
    interpolated_image = scipy.interpolate.griddata((x1, y1), newarr.ravel(), (xx, yy),method=method)
    return interpolated_image
