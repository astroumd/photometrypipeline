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
    # skyflat for filter, sky subtract if it exists using skypipeproc
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
        print("Bad pix: {}".format(np.sum(bad)))

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

        date = datetime.datetime.now().isoformat()
        head.add_history('Processed by skypipeprocmed ' + date)

        pf.writeto(outfile, data, head, overwrite=pipevar['overwrite'])

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')


def autopipemakesky(pipevar=inpipevar):
    """
    NAME:
        autopipemakesky
    PURPOSE:
        Combine sky flats based on filter type (sigma clipping for sources)
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
    for file in files:
        head = pf.getheader(file)
        filter = head['FILTER']
        filters += [filter]
    filters = np.array(filters)

    # Unique list of filters
    filterlist = set(filters)

    # For each unique filter, combine sky files using skycombine if more than 2 files
    # Otherwise return list of unprocessed files
    for filt in filterlist:
        skyflats = np.where(filters == filt)
        outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '.fits'

        if len(skyflats[0]) >= 2:

            if os.path.isfile(outflatname) and pipevar['overwrite'] == 0:
                print('Skipping makesky for ' + filt + '. File already exists')
                continue

            files = np.array(files)
            if pipevar['verbose']:
                print(filt, '-band sky flats.')
                print(files[skyflats])

                apd.skypipecombine_new(files[skyflats], outflatname, filt,
                                       #                apd.skypipecombine(files[skyflats], outflatname, file,
                                       pipevar, removeobjects=True, type='sky')
        else:
            print('Unable to produce a flat field for this setting: ' + filt)
            print('Will not be able to further process ' + str(len(skyflats)) + \
            ' image(s) without a flat from another source:')

            #for i in np.arange(len(skyflats[0])):
                #print(skyflats[0])
                #print('    ' + str(files[skyflats[i]]))

                # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')

def autopipemakesky_targets(pipevar=inpipevar):
    """
    NAME:
        autopipemakesky_targets
    PURPOSE:
        Combine sky flats based on filter type and source target (sigma clipping for sources)
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipemakesky_targets(pipevar=inpipevar)
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
        filter = head['FILTER']
        target = re.sub(r'\s+', '', head['TARGNAME'])
        filters += [filter]
        targets += [target]
    filters = np.array(filters)
    targets = np.array(targets)

    # Unique list of targets
    targetlist = set(targets)
    print(targetlist)

    # For each unique target and filter combination, combine sky files using skycombine_new if more than 2 files
    # Otherwise return list of unprocessed files
    for targ in targetlist:

        thistarget = np.where(targets == targ)
        if len(thistarget) == 0:
            print("issue!!!")
            continue

        thistargetfilts = set(filters[thistarget])
        print("{} Filters: {}".format(targ,thistargetfilts))

        for filt in thistargetfilts:
            skyflats = np.where(np.logical_and(filters == filt, targets == targ))
            outflatname = pipevar['imworkingdir'] + 'sky-' + filt + '_' + targ + '.fits'

            if len(skyflats[0]) >= 2:

                if os.path.isfile(outflatname) and pipevar['overwrite'] == 0:
                    print('Skipping makesky for ' + filt + '. File already exists')
                    continue

                files = np.array(files)
                if pipevar['verbose']:
                    print(filt, '-band sky flats.')
                    print(files[skyflats])

                    apd.skypipecombine_new(files[skyflats], outflatname, filt, pipevar, removeobjects=True, type='sky')
            else:
                print('Unable to produce a flat field for this setting: ' + filt)
                print('Will not be able to further process ' + str(len(skyflats)) + \
                ' image(s) without a flat from another source:')

                #for i in np.arange(len(skyflats[0])):
                    #print(skyflats[0])
                    #print('    ' + str(files[skyflats[i]]))

                    # If remove intermediate files keyword set, delete p(PREFIX)*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')


def autopipeskysub(pipevar=inpipevar):
    """
    NAME:
        autopipeskysub
    PURPOSE:
        Subtracts master sky flat from data and subtracts median.
    OPTIONAL KEYWORDS:
        pipevar  - input pipeline parameters (typically set in ratautoproc.pro,
                   but can be set to default)
    EXAMPLE:
        autopipeskysub(pipevar=inpipevar)
    DEPENDENCIES:
        autoproc_depend.skypipeproc
    """

    print('SKY-SUBTRACT')

    # Find data that needs to be sky subtracted
    files = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    sfiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    skys = glob.glob(pipevar['imworkingdir'] + '*sky-*.fits')

    if len(skys) == 0:
        print('No master sky files found, cannot sky subtract')
        return

    # Find the associated filter of each master skyflat
    skyfilts = []
    for sky in skys:
        head = pf.getheader(sky)
        filter = head['FILTER']
        skyfilts += [filter]

    # For each file if output files don't exist or override set check if we have master
    # skyflat for filter, sky subtract if it exists using skypipeproc
    for file in files:

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping sky subtraction for ' + file + '. File already exists')
            continue

        head = pf.getheader(file)
        filter = head['FILTER']

        # Find corresponding master skyflat
        try:
            skyloc = skyfilts.index(filter)
            skyfile = skys[skyloc]
        except:
            print('Sky field not found for ', file)
            pipevar['flatfail'] += ' ' + file
            continue

        if pipevar['verbose'] > 0:
            print('Sky subtracting', file, 'using', skyfile)

        apd.skypipeproc(file, skyfile, outfile)

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')

def autopipeskysub_targets(pipevar=inpipevar):
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
        autoproc_depend.skypipeproc
    """

    print('SKY-SUBTRACT')

    # Find data that needs to be sky subtracted
    files = glob.glob(pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
    sfiles = glob.glob(pipevar['imworkingdir'] + 'sfp' + pipevar['prefix'] + '*.fits')

    if len(files) == 0:
        print('Did not find any files! Check your data directory path!')
        return

    skys = glob.glob(pipevar['imworkingdir'] + '*sky-*.fits')

    if len(skys) == 0:
        print('No master sky files found, cannot sky subtract')
        return

    # Find the associated filter of each master skyflat
    skyset = []
    for sky in skys:
        head = pf.getheader(sky)
        filter = head['FILTER']
        target = [re.sub(r'\s+', '', head['TARGNAME'])]
        skyset += [(filter,target)]


    # For each file if output files don't exist or override set check if we have master
    # skyflat for filter, sky subtract if it exists using skypipeproc
    for file in files:

        fileroot = os.path.basename(file)
        outfile = pipevar['imworkingdir'] + 's' + fileroot

        if os.path.isfile(outfile) and pipevar['overwrite'] == 0:
            print('Skipping sky subtraction for ' + file + '. File already exists')
            continue

        head = pf.getheader(file)
        filter = head['FILTER']
        target = [re.sub(r'\s+', '', head['TARGNAME'])]

        # Find corresponding master skyflat
        try:
            skyloc = skyset.index((filter,target))
            skyfile = skys[skyloc]
        except:
            print('Sky field not found for ', file)
            pipevar['flatfail'] += ' ' + file
            continue

        if pipevar['verbose'] > 0:
            print('Sky subtracting', file, 'using', skyfile)

        apd.skypipeproc(file, skyfile, outfile)

    # If remove intermediate files keyword set, delete p(PREFIX)*.fits, fp(PREFIX)*.fits,
    # and sky-*.fits files
    if pipevar['rmifiles'] != 0:
        os.system('rm -f ' + pipevar['imworkingdir'] + 'p' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + 'fp' + pipevar['prefix'] + '*.fits')
        os.system('rm -f ' + pipevar['imworkingdir'] + '*sky-*.fits')