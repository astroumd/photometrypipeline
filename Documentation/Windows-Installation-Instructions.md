This installation uses Windows Subsystem for Linux with a Debian terminal. I recommend using Debian as it has ATLAS, and multiple AstrOmatic libraries as part of its package manager. ATLAS is notoriously difficult to install.

# Install Windows Dependencies

* Enable Windows Subsystem for Linux

![windows_features](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/windows_features.jpg)

* Install the Debian terminal from the [Microsoft Store](https://www.microsoft.com/en-us/p/debian/9msvkqc78pk6?activetab=pivot:overviewtab)

* [Install VcXsrv](https://sourceforge.net/projects/vcxsrv/)

# Open and Configure VcXsrv

Open XLaunch, and configure as shown.

<img src="https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_xlaunch_start.jpg" alt="windows_xlaunch_start" width="500"/>

![xlaunch settings](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_xlaunch_settings_1.jpg)

![xlaunch settings](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_xlaunch_settings_2.jpg)

![xlaunch settings](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_xlaunch_settings_3.jpg)

![xlaunch settings](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_xlaunch_settings_4.jpg)

### Recommended steps

* Save the configuration

* Place your config.xlaunch file in the Startup folder so that your X Server starts automatically

    * C:\Users\\**your_username**\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Startup

# Install Debian Environment Dependencies

* Open Debian Terminal

* Install libraries from the Debian package manager

      $ sudo apt-get update -y
      $ sudo apt-get install -y apt-utils git sextractor swarp scamp missfits gfortran libatlas-base-dev libpng-dev libfreetype6-dev libxft-dev tcl-dev tk-dev python-tk libsm-dev x11-apps saods9 make wget

* Link dependencies to match photopipe's expected commands

      $ sudo ln -s /usr/bin/sextractor /usr/bin/sex
      $ sudo ln -s /usr/bin/SWarp /usr/bin/swarp

* Install miniconda 2

      $ mkdir ~/installation_files
      $ cd ~/installation_files
      $ wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
      $ bash Miniconda2-latest-Linux-x86_64.sh

* Restart terminal to make conda commands available

* Install python dependencies

      $ conda create -n photopipe python=2.7 numpy scipy matplotlib astropy jsonschema ipython conda-build
      $ conda activate photopipe

* Fix ds9 installation bug

      $ cd ~/installation_files
      $ wget http://ds9.si.edu/download/debian9/ds9.debian9.8.1.tar.gz
      $ tar xzvf ds9.debian9.8.1.tar.gz
      $ sudo mv ds9 /usr/local/bin
      $ sudo rm /usr/bin/ds9

* Install remaining AstrOmatic dependency

      $ cd ~/installation_files
      $ wget http://cdsarc.u-strasbg.fr/ftp/pub/sw/cdsclient.tar.gz
      $ tar -zxvf cdsclient.tar.gz
      $ cd cdsclient-3.84/
      $ ./configure
      $ make
      $ sudo make install

# Install photopipe

* clone and build photopipe

      $ cd ~  # or wherever you would like this folder to live
      $ conda activate photopipe
      $ git clone https://github.com/astroumd/photometrypipeline.git
      $ conda develop photometrypipeline
      $ cd photometrypipeline
      $ python setup.py build

* Switch matplotlib interactive backend to get rid of QT error

      $ cd ~/installation_files
      $ conda activate photopipe
      $ pltrc=$(python -c "import matplotlib; print(matplotlib.matplotlib_fname())")
      $ sed 's/#backend : qt5agg/backend : TKagg/g' $pltrc > ./temprc
      $ cp $pltrc matplotlibrc_old
      $ cp ./temprc $pltrc

# Test the installation

* Test that X forwarding is working

      $ xlogo

    You should see the X server logo appear

* Test that matplotlib plots can be shown

    Start an ipython console

      $ ipython

    Create a plot

      from matplotlib import pyplot as plt
      
      x = [1, 2, 3]
      y = [4, 5, 6]

      plt.plot(x,y)
      plt.show()

* Test that the pipeline is working

      $ cd ~/photometrypipeline
      $ conda activate photopipe
      $ python test.py  # respond 'y' to the prompts

* Open the co-added image

      $ ds9 test/copy/reduced/coaddGRB190829A_20200129T021951840_r.fits
