Below is a full list of dependencies required for the pipeline to run.

* [Installation Tools](#installation-tools)
* [Python Packages](#python-packages)
* [AstrOmatic Software](#astromatic-software)
* [Additional Dependencies](#additional-dependencies)
* [Platform-specific Packages (Linux/Debian)](#platform-specific-packages-linuxdebian)

## Installation Tools

* [macports]() (macOS only) - Package manager
  * Required to install [AstrOmatic software](http://www.astromatic.net/software)
  * Requires also [Xcode and Xcode Developer Tools](https://guide.macports.org/#installing.xcode) 
* [pip](https://pip.pypa.io/en/stable/installing/) - Python package manager
  * Required to install `astropy` and other Python modules

## Python Packages

* [numpy](http://www.numpy.org/) - Fundamental package for scientific computing with Python
* [scipy](https://www.scipy.org/) - Open Source Library of Scientific Tools
* [matplotlib](http://matplotlib.org/) - 2D plotting library 
* [pyfits](http://www.stsci.edu/institute/software_hardware/pyfits) - Convenience functions that allow the user to work with FITS data
* [astropy](http://www.astropy.org/)-  A core package for astronomy

## AstrOmatic Software

* [SExtractor](http://www.astromatic.net/software/sextractor) - Source extraction
* [SCAMP](http://www.astromatic.net/software/scamp) - Astrometric calibration and photometric homogenisation
* [SWarp](http://www.astromatic.net/software/swarp) - Image regridding and co-addition
* [missFITS](http://www.astromatic.net/software/missfits) - FITS file management

## Additional Dependencies

* [CDSClient](http://cdsarc.u-strasbg.fr/doc/cdsclient.html) -  Query astronomical catalogues from the command line
* [CLAPACK](http://www.netlib.org/clapack/) - f2c'ed version of LAPACK
  * Required by SCAMP

## Platform-specific Packages (Linux/Debian)

* [gfortran](https://packages.debian.org/jessie/gfortran), [libatlas-base-dev](https://packages.debian.org/jessie/libatlas-base-dev)
  * Required by numpy, scipy, AstrOmatic
* [libatlas3gf-base](https://packages.debian.org/jessie/libatlas3gf-base), libatlas-headers
  * Required by SExtractor, SCAMP
* [liblapack-dev](https://packages.debian.org/jessie/liblapack-dev)
  * Required by CLAPACK
* [libshp-dev](https://packages.debian.org/jessie/libshp-dev), [libfftw3-dev](https://packages.debian.org/jessie/libfftw3-dev), [libplplot-dev](https://packages.debian.org/jessie/libplplot-dev)
  * Required by SCAMP
* [libpng12-dev](https://packages.debian.org/jessie/libpng12-dev), [libfreetype6-dev](https://packages.debian.org/jessie/libfreetype6-dev), [libxft-dev](https://packages.debian.org/jessie/libxft-dev)
  * Required by matplotlib
* [python-tk](https://packages.debian.org/jessie/python-tk)
  * Required by pyfits