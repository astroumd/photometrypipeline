Below are the step-by-step instructions for installing the pipeline on macOS (OS X), as an alternative to the package installation with `pip` or `python setup.py`.

_The instructions below are also provided as an automated [bash script](), however preliminary steps 1. and 2. below are still **required**._

1. [Installation Tools](#installation-tools)
2. [Python Environment](#python-environment)
3. [Python Packages](#python-packages)
4. [AstrOmatic Software](#astromatic-software)
5. [Download and Configure](#final-step-download-the-repository-and-create-the-configuration-file) 

___WARNING: compiling some dependencies will take several hours, see step 4.2.___
  

### 1. Installation tools

1. Download and Install [Xcode](https://itunes.apple.com/us/app/xcode/id497799835?mt=12) from the AppStore

2. Once you have Xcode installed, open a terminal, run `xcode-select --install`

 See [here](https://guide.macports.org/#installing.xcode) for more information.

2. [Install MacPorts](https://www.macports.org/install.php)


### 2. Python Environment

While macOS already ships with Python, we recommend to get the current stable production version. This can be done with MacPorts. Otherwise, jump to step 3.

1. In the terminal, run `sudo port selfupdate` to update MacPorts's list of Portfiles
2. To install Python, run `sudo port install python27`

3. Finally, we need `pip` to install `astropy` and other Python packages. In the terminal, run `sudo easy_install pip`

**Note:** If `easy_install` is not found, run `curl https://bootstrap.pypa.io/ez_setup.py -o - | sudo python` and repeat the last step.


### 3. Python Packages

1. Before installing the Python packages, run the following command to upgrade `pip`

 `python -m pip install --upgrade pip`

2. Run the following list of commands to install all the required packages:

 ```
 pip install --user numpy scipy matplotlib
 pip install --user pyfits
 pip install --user --no-deps astropy
 ```

**Note 1**: It's recommended to set the `--user` flag to install the packages for the local user only. Alternatively, you can perform a [systemwide installation](https://www.scipy.org/install.html#install-systemwide-via-a-mac-package-manager) with MacPorts.

**Note 2**: The `--no-deps` flag is set to tell pip to not upgrade numpy when installing astropy.

### 4. AstrOmatic Software

1. If you haven't yet, ensure that MacPorts is up-to-date 

 `sudo port selfupdate`

2. Install the [ATLAS](http://math-atlas.sourceforge.net/) library

 `sudo port install atlas +nofortran`
 
 ___WARNING: this will take several hours, and is best done overnight.___

2. Run the following list of commands to install the AstrOmatic software:
 ```sh
 sudo port install sextractor 
 sudo port install swarp
 sudo port install scamp
 sudo port install missfits
 sudo port install cdsclient
 ```

### Final Step: Download the repository and create the configuration file

1. [Download](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/archive/master.zip) and unzip the repository into your desired location, then rename the unzipped directory to `photometrypipeline`

2. Open the terminal and `cd` into the unzipped directory

3. Use the provided make_config script to create a configuration file and set the paths to the AstrOmatic software, as well as other default parameters

 ```
 sh ./make_config.sh
 ```

 The configuration file will be saved to `photometrypipeline/photometrypipeline/reduction/auto/pipeautoproc.par` where the first component of the path is the repository's root directory. 
