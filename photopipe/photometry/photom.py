"""
NAME:
	photom

PURPOSE:
	Creates a stacked image with all filters (saved to multicolor.fits and multicolor.weight.fits). 
	Crops all of the filter images and combined image to the same size and coordinates (file).ref.multi.fits.
	Then finds sources using the master combined image, and calculates the magnitude based on
	just the filtered images (with weight file using sextractor.  Saves new sextractor values to 
	"fluxes_(FILTER).txt'. Calculates absolute magnitude errors based on fluxes_(FILTER).txt and keyword 
	in file for absolute zeropoint RMS. Saves final magnitudes to 'finalphot(FILTER).am'

OUTPUTS:
	multicolor.[weights.]fits - files with all filter images stacked
	(file).ref.[weights.]fits - files resampled from swarp using all files
	(file).ref.multi.[weights.]fits - files cropped so that they include all filters
	coords(FILTER) - RA and DEC coordinates from sextractor from cropped images
	fluxes_(FILTER).txt - sextractor output from cropped images
	finalphot(FILTER).am - file containing pixel and coordinate location, mag and 
						   corrected mag error, flux and flux error of sextractor sources found (fluxes_(FILTER).txt)
	
Translated from icoords.pro by John Capone (jicapone@astro.umd.edu).
Modified by Vicki Toy (vtoy@astro.umd.edu) 5/21/2014
"""

import numpy as np
import os
import astropy.io.fits as pf
from photopipe.photometry.dependencies import photprocesslibrary as pplib
from photopipe.reduction.auto import autoproc_steps as ap
#from string import index
from numpy import shape
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
# Disable interactive mode
pl.ioff()


def photom(prefchar='coadd'):

	pipevar = {'getsedcommand':'get_SEDs', 'sexcommand':'sex' , 'swarpcommand':'swarp', 'refdatapath':'', 'defaultspath':'', 'imworkingdir':'' }
	ap.autopipedefaults(pipevar=pipevar)

	files = pplib.choosefiles(prefchar + '*.fits')
	filters = []
	for file in files:
		head = pf.getheader(file)
		filter = head['FILTER']
		filters += [filter]
	filters = np.array(filters)
	filters = set(filters)

	# Identify files (must have same number of images files as weight files)
	zffiles = []
	weightfiles = []
	for filter in filters:
		zffiles.extend(pplib.choosefiles(prefchar + '*_{}.fits'.format(filter)))
		#print(zffiles)
		weightfiles.extend(pplib.choosefiles(prefchar + '*_{}.weight.fits'.format(filter)))
		#print(weightfiles)
	coaddfiles = np.copy(zffiles)
	
	if len(zffiles) > len(weightfiles):
		print('Must have matching weight file to each image file to run automatic crop.')
		print('To use manual crop user manualcrop keyword and change crop values by hand')
		return -1
	
	numfiles = len(zffiles)
	print("Number of files {}".format(numfiles))

	# Resample all images using SWarp to a reference image called multicolor using weight files
	swarpstr = ''
	for i in range(numfiles):		
		swarpstr = swarpstr + zffiles[i] + ' '

	stackcmd = 'swarp ' + swarpstr + '-DELETE_TMPFILES N -WRITE_XML N -SUBTRACT_BACK N -WEIGHT_TYPE MAP_WEIGHT -IMAGEOUT_NAME multicolor.fits -WEIGHTOUT_NAME multicolor.weight.fits'
	stackcmd = stackcmd + ' -COPY_KEYWORDS OBJECT,TARGNAME,TELESCOP,FILTER,INSTRUME,OBSERVAT,PIXSCALE,ORIGIN,CCD_TYPE,JD,DATE-OBS,DATE1,DATEN,AIRMASS,TOTALEXP,FLATFLD,FLATTYPE,SEEPIX,ABSZPT,ABSZPTSC,ABSZPRMS'
	print(stackcmd)
	os.system(stackcmd)

	# Rename all the resampled files to crop files
	for i in range(numfiles):
		tmp = zffiles[i].split('.')[0]
		ifile = tmp + '.resamp.fits'
		ofile = tmp + '.ref.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		print(mvcmd)
		os.system(mvcmd)
		
		ifile = tmp + '.resamp.weight.fits'
		ofile = tmp + '.ref.weight.fits'
		mvcmd = 'mv -f ' + ifile + ' ' + ofile
		os.system(mvcmd)

	coaddreffiles = pplib.choosefiles(prefchar+'*.ref.fits')

	ra1arr  = []
	dec1arr = []
	ra2arr  = []
	dec2arr = []

	# Finds the RA and DEC of the first and the last pixel of each cropped coadded file
	for files in coaddreffiles:

		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		imSize = shape(data)
	
		# Converts pixel value to RA and DEC using header information (AstroPy function)
		w = wcs.WCS(fitsheader)
		pixcrd = [[0.,0.], [imSize[1]-1.0, imSize[0]-1.0]]
		[[ra1, dec1], [ra2, dec2]] = w.wcs_pix2world(pixcrd, 0)
	
		# Stores information into arrays
		ra1arr.append(ra1)
		dec1arr.append(dec1)
		ra2arr.append(ra2)
		dec2arr.append(dec2)

	# Finds the coordinates that fit all of the data
	raleft  = min(ra1arr)
	raright = max(ra2arr)
	decbot  = max(dec1arr)
	dectop  = min(dec2arr)

	# Crops data so the size of every filter image matches and saves to file 'coadd*.multi.fits'
	# Same for weight file
	for files in coaddreffiles:
		newfile = files[:-4]+'multi.fits'
		fitsfile = pf.open(files)
		fitsheader = fitsfile[0].header
		data = fitsfile[0].data
	
		w = wcs.WCS(fitsheader)
		[[x1,y1],[x2,y2]] = w.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
		
		pplib.hextractlite(newfile, data, fitsheader, x1+1, x2, y1+1, y2)
		
		wnewfile = files[:-4]+'multi.weight.fits'
		wfitsfile = pf.open(files[:-4]+'weight.fits')
		wfitsheader = wfitsfile[0].header
		wdata = wfitsfile[0].data
	
		w = wcs.WCS(wfitsheader)
		[[x1,y1],[x2,y2]] = w.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
		
		pplib.hextractlite(wnewfile, wdata, wfitsheader, x1+1, x2, y1+1, y2)

	# Crops the multicolor (data file and weight) fits files to match the same coordinates
	mixfile = 'multicolor.fits'
	mixfitsfile = pf.open(mixfile)
	mixfitsheader = mixfitsfile[0].header
	mixdata = mixfitsfile[0].data

	mixw = wcs.WCS(mixfitsheader)
	[[mx1,my1],[mx2,my2]] = mixw.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
	pplib.hextractlite(mixfile, mixdata, mixfitsheader, mx1+1, mx2, my1+1, my2)

	wmixfile = mixfile[:-4]+'weight.fits'
	wmixfitsfile = pf.open(wmixfile)
	wmixfitsheader = wmixfitsfile[0].header
	wmixdata = wmixfitsfile[0].data

	wmixw = wcs.WCS(wmixfitsheader)
	[[wmx1,wmy1],[wmx2,wmy2]] = mixw.wcs_world2pix([[raleft, decbot], [raright, dectop]], 0)
	pplib.hextractlite(wmixfile, wmixdata, wmixfitsheader, wmx1+1, wmx2, wmy1+1, wmy2)	

	# Find directory where this python code is located
	propath = os.path.dirname(os.path.realpath(__file__))
	
	# Make sure configuration file is in current working directory, if not copy it from
	# location where this code is stored
	if not os.path.exists('ratir_weighted.sex'): 
		os.system('cp '+propath+'/defaults/ratir_weighted.sex .')
	
	if not os.path.exists('weightedtemp.param'): 
		os.system('cp '+propath+'/defaults/weightedtemp.param .')
		
	if not os.path.exists('ratir.conv'): 
		os.system('cp '+propath+'/defaults/ratir.conv .')
				
	if not os.path.exists('ratir_nir.nnw'): 
		os.system('cp '+propath+'/defaults/ratir_nir.nnw .')

	#coaddfiles = pplib.choosefiles(prefchar+'*_?.fits')

	# Uses sextractor to find the magnitude and location of sources for each file
	# Saves this information into 'fluxes_*.txt' files
	for files in coaddfiles:

		hdr = pf.getheader(files)

		try:
			filter = hdr['FILTER']
			abszpt = hdr['ABSZPT']
			abszprms = hdr['ABSZPRMS']
			pixscale = hdr['PIXSCALE']
		except KeyError as error:
			print(error)
			continue

		# Finds filter name and makes sure it is capitalized correctly
		filter = files.split('_')[-1].split('.')[0]

		if filter in ('SDSS-U', 'SDSS-G', 'SDSS-R', 'SDSS-I', 'SDSS-Z'):
			filter = filter[-1].lower()
		elif filter.lower() in ('j','h','k'):
			filter = filter.upper()
		else:
			filter = filter.lower()
		#print(filter)
		compfile = files

		# Use individual image, not multi image now for individual detections
		# Call to sextractor in double image mode (image1 used for detection of sources, image2 only for measurements - must be same size)
		os.system('sex -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE '+compfile[:-4]+'weight.fits' + \
			' -c ratir_weighted.sex -SEEING_FWHM 1.5 -PIXEL_SCALE '+str(pixscale)+' -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -PHOT_APERTURES ' + \
			str(hdr['SPIX'])+ ' -MAG_ZEROPOINT ' + str(hdr['ABSZPT'])+' ' + compfile)
		os.system('mv -f temp.cat fluxes_'+filter+'.txt')
		
		# Columns unpacked for fluxes*.txt are: (x,y,ra,dec,mag,magerr,flux,fluxerr,e,fwhm,flags)
		sexout = np.loadtxt('fluxes_'+filter+'.txt', unpack=True)
		
		magcol = 4
		magerrcol = 5

		# SEDs
		# Create catalog star file
		xim  = sexout[0,:]
		yim  = sexout[1,:]

		w = wcs.WCS(hdr)
		wrd = w.all_pix2world(np.transpose([xim, yim]), 0)
		imfile  = files + '.seds.im'
		catfile = files + '.seds.cat'

		# Save stars from image
		np.savetxt(imfile, np.transpose([wrd[:,0],wrd[:,1],sexout[magcol,:]]))

		sedcmd = 'python ' + pipevar['getsedcommand'] + ' ' + imfile + ' ' +\
			filter + ' ' + catfile + " 15 True True"
		os.system(sedcmd)
		#
		# tout = np.transpose(sexout[0:8,:]) #Only include through fluxerr
		indexes = np.array([np.arange(len(sexout[0,:]))])
		tout = np.hstack((sexout[0:8,:].T, indexes.T))

		for i in np.arange(len(tout[:,magerrcol])):
			tout[i,magerrcol] = max(tout[i,magerrcol], 0.01)
		
		tout[:,magerrcol] = np.sqrt(tout[:,magerrcol]**2 + abszprms**2)

		tsorted =  tout[np.argsort(tout[:,magcol])]
			
		# Creates Absolute Magnitude file with coordinates
		amfile = 'finalphot'+filter+'.am'
		np.savetxt(
			amfile, tsorted, fmt='%15.6f ' * 8 + '%i',
			header='X\t Y\t RA\t DEC\t CAL_MAG\t CAL_MAG_ERR\t CAL_FLUX\t CAL_FLUX_ERR\t CAT_INDEX\t'
		)

		# Creates Region File using RA, DEC, and FWMH
		print("Creating region file for filter : {}".format(filter))
		fwmh = sexout[9,:] #FWMH each source in pixels
		region_file = open("region_" + filter +".REG", "w")
		a = '# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
		region_file.write(a)
		for i in range(len(wrd[:,0])):
			c = SkyCoord(ra=wrd[i,0] * u.degree, dec=wrd[i,1] * u.degree)
			ra, dec = c.to_string('hmsdms').split()
			ra = ra.replace("h", ":").replace("m", ":").replace("s", "")
			dec = dec.replace("d", ":").replace("m", ":").replace("s", "")
			region_file.write('\ncircle({},{},{:.4f}")'.format(ra, dec, fwmh[i]*pixscale))
		region_file.close()

		# Plot Astrometry Error
		print("Creating astronomy error file for filter : {}".format(filter))
		cvars = np.loadtxt(catfile, unpack=True)
		ra_cat = cvars[0, :]
		dec_cat = cvars[1, :]
		ra_im = sexout[2, :]
		dec_im = sexout[3, :]
		nsources = len(ra_cat)
		ra_err = np.zeros(nsources)
		dec_err = np.zeros(nsources)
		for i in range(nsources):
			ra_err[i] = ra_im[i] - ra_cat[i]
			#print("Im:{:.4f}, Cat:{:.4f}, Err:{:.4f}".format(ra_im[i],ra_cat[i],ra_err[i]))
			dec_err[i] = dec_im[i] - dec_cat[i]
			#print("Im:{:.4f}, Cat:{:.4f}, Err:{:.4f}".format(dec_im[i], dec_cat[i], dec_err[i]))

		fig = pl.figure(1, figsize=(8, 8))
		ax1 = fig.add_subplot(211)
		#ax1.set_xlim([min(ra_im), max(ra_im)])
		#ax1.set_ylim([min(ra_err),max(ra_err)])
		ax1.plot(ra_im, ra_err, marker='o', linestyle='None')
		ax1.set_xlabel('RA (deg)')
		ax1.set_ylabel(r"$\Delta$RA (deg)")
		ax2 = fig.add_subplot(212)
		# ax2.set_xlim([min(ra_im), max(ra_im)])
		# ax2.set_ylim([min(ra_err),max(ra_err)])
		ax2.plot(dec_im, dec_err, marker='o', linestyle='None')
		ax2.set_xlabel('DEC (deg)')
		ax2.set_ylabel(r"$\Delta$DEC (deg)")
		pl.savefig('astronomy_stats_' + filter + '.png', bbox_inches='tight')
		pl.clf()