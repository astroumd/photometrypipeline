## Installation

The easiest way to get up and running with the pipeline is to run the ready-to-use **docker container**.

Alternatively, it can be installed on Linux or macOS either with `pip`, or with the provided installation scripts.
_(**WARNING:** compiling the dependencies can take more than 6 hours.)_

#### 1) Download a pre-configured docker container (joedurbak/photopipe)

The docker container can be used if you have docker installed by running:

    $ docker run -it joedurbak/photopipe bash
    
More detailed instructions on installing, configuring, and rerunning the same container are available for both [Mac](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Docker-Instructions-on-Mac.md) and [Windows](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Docker-Instructions-on-Windows.md) operating systems

#### 2) Install on your machine from PyPI or git

* Run `sudo -H pip install photopipe` to install the latest stable version from [PyPI](https://pypi.python.org/pypi/photopipe). 

* Or, clone from `git`:

 ```
 $ git clone https://github.com/astroumd/photometrypipeline.git
 $ cd photometrypipeline
 $ sudo python setup.py install
 ```
 
**NOTE (macOS)**: If the installation fails with `sudo: port: command not found` make sure that [MacPorts](https://guide.macports.org/#installing) is installed and `/opt/local/bin` is in the $PATH (e.g. `export PATH=/opt/local/bin:/opt/local/sbin:$PATH`).

#### 3) Manual Installation
If you to run the pipeline from the Python environment rather than using the `photopipe` command as described in the [Usage](#usage) section, please follow the step by step instructions in the wiki pages below to install all the dependencies manually.

* [Instructions for macOS](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Mac-Installation-Instructions.md)
* [Instructions for Linux-Debian](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Linux-Debian-Installation-Instructions.md)
* [Instructions for Windows-WSL-Debian](https://github.com/RIMAS-RATIR-DCT/photometrypipeline/blob/master/Documentation/Windows-Installation-Instructions.md)

