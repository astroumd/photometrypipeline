Below are the step-by-step instructions for installing the pipeline on Debian-based distributions of Linux, as an alternative to the package installation with `pip` or `python setup.py`.

_The instructions below are also provided as an automated [bash script]()_

1. [Python Environment](#python-environment)
2. [Python Packages](#python-packages)
3. [AstrOmatic Software](#astromatic-software)
4. [Download and Configure](#final-step-download-the-repository-and-create-the-configuration-file) 

___WARNING: compiling some dependencies will take several hours, see step 4.2.___
  

### 1. Python Environment

1. In the terminal, run `sudo apt-get update` to download the updated lists of packages
2. Install Python 2.7 (if you already have Python, make sure also `python2.7-dev` is installed)

 ```
 sudo apt-get install python2.7 python2.7-dev
 sudo ln -s /usr/bin/python2.7 /usr/bin/python
 ```
3. Finally, we need `pip` to install `astropy` and other Python packages

 ```
 sudo apt-get install python-pip
 ```

### 2. Python Packages

1. Before installing the Python packages, run the following command to upgrade `pip`

 `sudo pip install -U pip`

2. Run the following list of commands to install all the required packages

 **Numpy and Scipy**

 ```
 sudo apt-get install gfortran libatlas-base-dev
 pip install --user numpy scipy
 ```
 **Matplotlib**

 ```
 sudo apt-get install libpng12-dev libfreetype6-dev libxft-dev
 pip install --user matplotlib
 ```

 **Astropy**

 ```
 pip install --user --no-deps astropy
 ```

 **PyFITS**

 ```
 sudo apt-get install python-tk
 pip install --user pyfits
 ```

**Note 1**: It's recommended to set the `--user` flag to install the packages for the local user only. Alternatively, you can perform a [systemwide installation](https://www.scipy.org/install.html#install-systemwide-via-a-linux-package-manager) with `apt-get`.

**Note 2**: The `--no-deps` flag is set to tell pip to not upgrade numpy when installing astropy.

### 4. AstrOmatic Software

We recommend to compile the software directly from source, in order to get the latest release.

1. **Install dependencies**
 ```
 sudo apt-get install libatlas3gf-base libatlas-headers libfftw3-dev liblapack-dev libplplot-dev libshp-dev
 ```

2. **SExtractor** 

 ```
 wget http://www.astromatic.net/download/sextractor/sextractor-2.19.5.tar.gz
 tar -zxvf sextractor-2.19.5.tar.gz
 cd sextractor-2.19.5
 ./configure
 make
 sudo make install
 ```

3. **SWarp** 

 ```
 wget http://www.astromatic.net/download/swarp/swarp-2.38.0.tar.gz
 tar -zxvf swarp-2.38.0.tar.gz
 cd swarp-2.38.0
 ./configure
 make
 sudo make install
 ```

4. **CDSClient**

 ```
 wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz 
 tar -zxvf cdsclient.tar.gz 
 cd cdsclient-3.83/ 
 ./configure 
 make 
 sudo make install
 ```

5. **CLAPACK** (required by SCAMP)

 1. Download and decompress the source

   ```
   wget http://www.netlib.org/clapack/clapack.tgz
   tar -zxvf clapack.tgz
   cd CLAPACK-3.2.1
   ```
 
 2. Rename make.inc.example to make.inc

   `mv make.inc.example make.inc`

 3. Open the file and replace 
   `BLASLIB      = ../../blas$(PLAT).a`

    with `BLASLIB      = ../../libcblaswr.a -lcblas -latlas`

    or, run `sed -i 's/blas$(PLAT).a/libcblaswr.a -lcblas -latlas/' make.inc`

 4. Run make and copy libraries to system folders	
 
   ```
   make cblaswrap
   sudo make lapacklib
   mv lapack_LINUX.a liblapack.a
   sudo cp libcblaswr.a /usr/local/lib/
   sudo cp liblapack.a /usr/local/lib/
   ```

6. **SCAMP**

 ```
 wget http://www.astromatic.net/download/scamp/scamp-2.0.4.tar.gz
 tar -zxvf scamp-2.0.4.tar.gz
 cd scamp-2.0.4
 ./configure --with-atlas-incdir=/usr/include/atlas/
 sudo make
 sudo make install
 ```

7. **MissFITS**

 ```
 wget http://www.astromatic.net/download/missfits/missfits-2.8.0.tar.gz
 tar -zxvf missfits-2.8.0.tar.gz
 cd missfits-2.8.0
 ./configure
 sudo make
 sudo make install
 ```

### Final Step: Download the repository and create the configuration file

1. [Download](https://github.com/maxperry/photometrypipeline/archive/master.zip) and unzip the repository into your desired location, then rename the unzipped directory to `photometrypipeline`

2. Open the terminal and `cd` into the unzipped directory

3. Use the provided [make_config]() script to create a configuration file and set the paths to the AstrOmatic software, as well as other default parameters

 ```
 sh ./make_config.sh
 ```

 The configuration file will be saved to `photometrypipeline/photometrypipeline/reduction/auto/pipeautoproc.par` where the first component of the path is the repository's root directory. 