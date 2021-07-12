These instructions assume that you have little, to no experience with docker. 

It's worth noting that configuring Windows to work with Linux Docker containers is a bit more complicated than if you were doing this on a Unix or Linux OS. This will likely improve when Windows Subsystem for Linux 2 (WSL2) is released widely. I would also note that there are a lot of useful tools available through docker, so it's worth going through this process.

# Install dependencies

* Make sure you have Windows 10 Professional edition. Docker doesn't work on Windows 10 Home edition. You can check this by following these [instructions](https://support.microsoft.com/en-us/help/13443/windows-which-version-am-i-running).

* Enable necessary Windows Features

    * Containers

    * Hyper-V

    * SMB 1.0/CIFS File Sharing Support

    * SMB Direct

    * Windows PowerShell 2.0

![windows_features](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/windows_features.jpg)

* Create a Virtual Switch with Hyper-V using these [instructions](https://docs.docker.com/machine/drivers/hyper-v/#example)

* [Install docker](https://docs.docker.com/docker-for-windows/install/)

    * Choose Linux containers

* [Install VcXsrv](https://sourceforge.net/projects/vcxsrv/)

# Start the Necessary Programs

## Start the X Server

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

## Start the Docker Daemon

* Open Docker Desktop in the Start menu

<img src="https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/window_docker_desktop_start.jpg" alt="windows_docker_desktop_start" width="500"/>

* Check to see that the docker whale icon is with your notification icons, and check for the X server icon while you're at it

![docker icon](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/windows_docker_icon.jpg)

* Open docker settings

![docker icon](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/windows_docker_icon_settings.jpg)

* Make your drive containing your data available for file sharing

![docker settings](https://github.com/astroumd/photometrypipeline/blob/master/wiki_files/windows_docker_settings.jpg)

# Fetch and run the docker container

* Open PowerShell

* Find out what your IP Address is by running `ipconfig`

    * Note: If you have multiple IP Addresses, you may need to go through some trial and error to determine which one your X Server is bound to.

* Now we're ready to run the docker container in the PowerShell by replacing `/path/to/your/data/directory` and `IPADDRESS` with the necessary values and running:

        PS C:\> docker run -it -v /path/to/your/data/directory:/mnt/data -e DISPLAY="IPADDRESS:0.0" --name photopipe1 joedurbak/photopipe

    This will start a bash shell within the container.
    
    * Note: on Windows, you will still need to use "/" instead "\\" for your file paths.

        * Example: If I put the data in a directory called Data in my Documents folder, I would replace '/path/to/your/data/directory' as follows

                PS C:\> docker run -it -v c:/Users/durba/Documents/Data:/mnt/data -e DISPLAY="IPADDRESS:0.0" --name photopipe1 joedurbak/photopipe

Now that this is working, let's run some tests.

* Check that your files are available within the container by running
      
        # ls /mnt/data

* Test that x-forwarding is working by running `xlogo`, and seeing if the logo appears

* Test that the pipeline is working by running, the following and replying 'y' to all of the prompts:

        # cd /work/photometrypipeline

        # python test.py


* Check the co-added output by running:

      # ds9 /work/photometrypipeline/test/copy/reduced/coaddGRB190829A_20200129T021951840_r.fits

* Exit the container by running `exit`

# Using the Container Again

Now that you've downloaded and run this container, you can use it again by running:

    PS C:\> docker start -i photopipe1

Also, if you want to create a new container, you just need to do the same docker run command as before, but with a new name. For example:

    PS C:\> docker run -it -v /path/to/your/data/directory:/mnt/data -e DISPLAY="IPADDRESS:0.0" --name photopipe2 joedurbak/photopipe
