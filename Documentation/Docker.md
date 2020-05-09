If you already have docker installed, and configured for X forwarding. You can use the container as follows.

# Fetch and run the docker container

Now we're ready to run the docker container in the terminal by running:

    $ docker run -it -e DISPLAY=host.docker.internal:0 -v /path/to/your/data/directory:/mnt/data --name photopipe1 joedurbak/photopipe

This will start a bash shell within the container.

Now that this is working, let's run some tests.

* Test that x-forwarding is working by running `xlogo`, and seeing if the logo appears

* Test that the pipeline is working by running, the following and replying 'y' to all of the prompts:

      # cd /work/photometrypipeline
      # python test.py


* Check the co-added output by running:

      # ds9 /work/photometrypipeline/test/copy/reduced/coaddGRB190829A_20200129T021951840_r.fits

* Check that your files are available within the container by running

      # ls /mnt/data

* Exit the container by running `exit`

# Using the Container Again

Now that you've downloaded and run this container, you can use it again by running:

    $ docker start -i photopipe1

Also, if you want to create a new container, you just need to do the same docker run command as before, but with a new name. For example:

    $ docker run -it -e DISPLAY=host.docker.internal:0 -v /path/to/your/data/directory:/mnt/data --name photopipe2 joedurbak/photopipe
