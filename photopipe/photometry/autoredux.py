from photopipe.photometry import photom
from photopipe.photometry import plotphotom
import glob
import astropy.io.fits as pf
import re
import numpy as np
import os
from shutil import move


def autoredux(data_path, prefchar='coadd', noplot=False):

	# Must be in directory with Coadd Frames
	os.chdir(data_path)
	print(data_path)

	# Find Targets
	files = glob.glob(prefchar + '*.fits')
	if len(files) == 0:
		print('Did not find any files! Check your data directory path!')
		return
	targets = []
	for file in files:
		head = pf.getheader(file)
		target = re.sub(r'\s+', '', head['TARGNAME'])
		targets += [target]
	targets = np.array(targets)
	targetlist = set(targets)
	print("Making Target Directories:")
	print(targetlist)

	# Creat directory for each target and organize
	for target in targetlist:
		targ_dir = os.path.join(data_path, target)
		try: os.mkdir(targ_dir)
		except FileExistsError: pass
		targ_files = glob.glob(prefchar + target + '*')
		for file in targ_files:
			f_base = os.path.join(data_path, file)
			f_new = os.path.join(targ_dir, file)
			move(f_base,f_new)

	# Begin Photometry
	for target in targetlist:
		targ_dir = os.path.join(data_path, target)
		os.chdir(targ_dir)
		print('photom: {}'.format(target))
		photom.photom()
		if noplot is False:
			print('plotphotom: {}'.format(target))
			plotphotom.plotphotom()
		print('Processing complete for {}. Results saved to photom.html and photom.json'.format(target))
