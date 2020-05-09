
# PhotoPipe

### Photometry Pipeline for RATIR and RIMAS

[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![PyPI version](https://badge.fury.io/py/photopipe.svg)](https://badge.fury.io/py/photopipe)

PhotoPipe is a pipeline for automated reduction, photometry and astrometry, written in Python (2.7), and designed to process imaging data from the following instruments:

* [RATIR](http://butler.lab.asu.edu/RATIR/): the Reionization and Transients Infrared/Optical Project.
* [RIMAS](https://lowell.edu/research/research-facilities/4-3-meter-ldt/rimas/): the Rapid infrared IMAger Spectrometer.
* [LMI](https://lowell.edu/research/research-facilities/4-3-meter-ldt/the-large-monolithic-imager/): the Large Monolithic Imager.


This project is based on the previous versions of the pipeline by [cenko](https://github.com/cenko/RATIR-GSFC) and [vickitoy](https://github.com/vickitoy/photometry_pipeline).  
Built at the NASA Goddard Space Flight Center, in collaboration with the University of Maryland.
<br><br><br>
![NASA GSFC - RATIR - UMD](https://github.com/maxperry/photometrypipeline/raw/master/docs/readme-logos.jpg)


## Full Documentation

See [this] (https://github.com/RIMAS-RATIR-DCT/photometrypipeline/tree/master/Documentation) for full documentation, examples, operational details and other information.


## Prerequisites

The pipeline can be installed with very little additional software. However, depending on the installation, different libraries and software may be required. Please see the [Installation](#installation) section below, and the dedicated [file](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Prerequisites.md) for a full list of prerequisites.


## Installation

The easiest way to get up and running with the pipeline is to run the ready-to-use **docker container**.

Alternatively, it can be installed on Linux or macOS either with `pip`, or with the provided installation scripts.
_(**WARNING:** compiling the dependencies can take more than 6 hours.)_

#### 1) Download a pre-configured Virtual Machine (PhotoPipe-VM)

The docker container can be used if you have docker installed by running:

    $ docker run -it joedurbak/photopipe bash
    
More detailed instructions on installing, configuring, and rerunning the same container are available for both [Mac](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Docker-Instructions-on-Mac.md) and [Windows](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Docker-Instructions-on-Windows.md) operating systems

#### 2) Install on your machine from PyPI or git

* Run `sudo -H pip install photopipe` to install the latest stable version from [PyPI](https://pypi.python.org/pypi/photopipe). 

* Or, clone from `git`:

 ```
 $ git clone https://github.com/RIMAS-RATIR-DCT/photometrypipeline.git
 $ cd photometrypipeline
 $ sudo python setup.py install
 ```
 
**NOTE (macOS)**: If the installation fails with `sudo: port: command not found` make sure that [MacPorts](https://guide.macports.org/#installing) is installed and `/opt/local/bin` is in the $PATH (e.g. `export PATH=/opt/local/bin:/opt/local/sbin:$PATH`).

#### 3) Manual Installation
If you to run the pipeline from the Python environment rather than using the `photopipe` command as described in the [Usage](#usage) section, please follow the step by step instructions in the wiki pages below to install all the dependencies manually.

* [Instructions for macOS](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Mac-Installation-Instructions.md)
* [Instructions for Linux-Debian](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Linux-Debian-Installation-Instructions.md)
* [Instructions for Windows-WSL-Debian](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Windows-Installation-Instructions.md)


## Usage

The following steps can be reproduced using the test data downloadable [here](https://drive.google.com/file/d/0BzMOBEOpFL9LaHpkWnFXc0IzRmM/view?usp=sharing), either from your host machine or from the virtual machine.

#### 1) Create a new directory with the following structure:
 
 **NOTE for PhotoPipe-VM**: Create the `imdata` directory in the VM's shared data dir (should be `./photopipe-vm/data`) 
 
 ```
imdata  
│
└───bias
│   │   20160628T032914C0b.fits
│   │   20160628T032914C1b.fits
│   │   ...
│   
└───dark
│   │   20160628T040211C0d.fits
│   │   20160628T040211C1d.fits
│   │   ...   
│
└───flat
│   │   20160628T024207C0f.fits
│   │   20160628T024207C1f.fits
│   │   ...    
│
└───science
│   │   20160628T043940C0o.fits
│   │   20160628T043940C1o.fits
│   │   ...    
│
└───science_selected
│
└───reduced
```

**If running on your Host Machine**: Open the terminal and `cd ./imdata/reduced`

**If running on PhotoPipe-VM**: Launch the VM first ([see instructions](https://github.com/maxperry/photometrypipeline-vm#usage)), and `cd /vagrant_data/imdata/reduced`

**NOTE**: It's important to perform the steps below from the `reduced` folder, otherwise make sure to move the master frames there after running `mkmaster` (it saves to the current dir!). 

#### 2) Preprocessing
 1. Enter Python environment: `$ python`
 2. Run the following script:
 
  ```python
  from photopipe.reduction import preproc
  
  # Bias frames calibration 
  
  bias_calib = preproc.choose_calib( 'ratir', 
                                     'bias', 
                                     workdir='/vagrant_data/imdata/bias/', 
                                     cams=[0,1], 
                                     auto=True, 
                                     amin=0.0, amax=1.0, 
                                     reject_sat=False, 
                                     save_select=True, 
                                     noplot=False )
 
 
  # Dark frames calibration
  
  dark_calib = preproc.choose_calib( 'ratir', 
                                     'dark', 
                                     workdir='/vagrant_data/imdata/dark/', 
                                     cams=[0,1], 
                                     auto=True, 
                                     amin=0.0, amax=1.0, 
                                     reject_sat=False, 
                                     save_select=True, 
                                     noplot=False )

  # Flat frames calibration 
  
  flat_calib = preproc.choose_calib( 'ratir', 
                                     'flat', 
                                     workdir='/vagrant_data/imdata/flat/', 
                                     cams=[0,1,2,3], 
                                     auto=True, 
                                     amin=0.2, amax=0.8, 
                                     reject_sat=False, 
                                     save_select=True, 
                                     noplot=False )
                                    
  # WARNING: RATIR flats are often bad even when 
  # the median value is in the acceptable range. 
  # Auto mode is only recommended for bias frame 
  # selection.     
  
  
  # Select science frames
  # (selected frames will be copied to target_dir)
  
  science_dict = preproc.choose_science( 'ratir', 
                                         workdir='/vagrant_data/imdata/science, 
                                         targetdir='/vagrant_data/imdata/science_selected', 
                                         cams=[0,1,2,3], 
                                         auto=True, 
                                         save_select=True, 
                                         calibrate=False, 
                                         noplot=False ) 

  # WARNING: When auto is True, all science frames
  # are selected. Since the telescope occasionally
  # has tracking issues, it is recommended to check
  # all frames.  
  
  
  # Make master frames
  # (saves to the current dir)
  
  preproc.mkmaster('ratir', bias_calib, 'bias')

  preproc.mkmaster('ratir', dark_calib, 'dark')
 
  preproc.mkmaster('ratir', flat_calib, 'flat')  
 ```
 
 See [documentation](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/preproc.py.md) for `preproc` functions reference.

 **WARNING**:
  1. Always use absolute paths terminated by a `/` when setting `dir` parameters. 
  2. Make sure to set the correct number of **cameras** where required (e.g. `cams[0,1,2,3`).
  3. `amin/amax` set the min and max **saturation** values. Make sure to select the appropriate values for each type of frame, only frames with median values in this range will be selected.
  
#### 3) Reduction
 1. Start a new python environment
 2. Execute the script below:
 
  ```python
   from photopipe.reduction.auto.autoproc import autoproc

   autoproc( datadir='/vagrant_data/imdata/science_selected/', 
              imdir='/vagrant_data/imdata/reduced/', 
              redo=1 )
  ```
 3. Reduced frames will be saved to `imdir` using `coadd` as prefix in the filename.
 
 See [documentation](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/autoproc.py.md) for `autoproc` function reference.

#### 4) Photometry

 1. Move the `coadd` files to a new dir called `photometry`, and cd into it
 2. Start a new python environment, and execute the script below:
 
  ```python
   from photopipe.photometry.autoredux import autoredux
   autoredux()
  ```
 **NOTE**: It's important to start python from the `photometry` dir containing the `coadd` files.
 
 **[autoredux](https://github.com/RIMAS_RATIR_DCT/photometrypipeline/blob/master/photopipe/photometry/autoredux.py)** will run the following scripts:

  1.  **[photom.py](https://github.com/RIMAS_RATIR_DCT/photometrypipeline/blob/master/photopipe/photometry/photom.py)**: Samples and crops all files, creates a multicolor image used to find all sources, then finds the *aperture photometry* for each of those sources (resampled) using sextractor (with corrected zeropoint).  
        - Output: *.am (absolute magnitude) files for w/ RA and DEC identified by sextractor
  2. **[plotphotom.py](https://github.com/RIMAS_RATIR_DCT/photometrypipeline/blob/master/photopipe/photometry/plotphotom.py)**: Plots photometry results to a HTML page.
 
## Bugs and Feedback

For bugs, questions and discussions please use the [Github Issues](https://github.com/maxperry/photometrypipeline/issues).


## License
This project is licensed under the GNU GPLv3 License - see the [LICENSE](https://github.com/scrapy/scrapy/blob/master/LICENSE) file for details.
