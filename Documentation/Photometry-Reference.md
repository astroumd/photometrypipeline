Below is a list of all functions that are part of the photometry module.

* [autoredux](#autoredux) -  **automated photometry script** 
  - Runs `photom` and `plotphotom`.
* [photom](#photom) -  **main photometry script**
  - Finds sources with `sextractor`, and calculates magnitudes and corrected magnitude errors of each source.
* [plotphotom](#plotphotom) -  **plotting script**
  - Plots overlay of all filters, images of each filter with annotated sources, and HTML report with all information.
* [printhtml](#printhtml) - **HTML report script**

##photom
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/photometry/photom.py "Source") **photom**(<i>[prefchar='coadd']</i>)

1. Creates a stacked image with all filters (saved to `multicolor.fits` and `multicolor.weight.fits`). 
2. Crops all of the filter images and combined image to the same size and coordinates `(file).ref.multi.fits`. 
3. Then finds sources using the master combined image, and calculates the magnitude based on just the filtered images (with weight file using `sextractor`.  Saves new `sextractor` values to `fluxes_(FILTER).txt`.
4. Calculates absolute magnitude errors based on `fluxes_(FILTER).txt` and keyword in file for absolute zeropoint RMS. Saves final magnitudes to `finalphot(FILTER).am`.

####Params

* `prefchar` **{string}**: Filename prefix of output files from reduction module.

####Output
* `multicolor.[weights.]fits` - files with all filter images stacked.
* `(file).ref.[weights.]fits` - files resampled from `swarp` using all files.
* `(file).ref.multi.[weights.]fits` - files cropped so that they include all filters.
* `coords(FILTER)` - RA and DEC coordinates from `sextractor` from cropped images.
* `fluxes_(FILTER).txt` - `sextractor` output from cropped images.
* `finalphot(FILTER).am` - file containing pixel and coordinate location, mag and 						   corrected mag error of `sextractor` sources found (fluxes_(FILTER).txt)

##plotphotom
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/photometry/plotphotom.py "Source") **plotphotom**(<i>[prefchar='coadd']</i>)

1. Reads in `prefchar*(FILTER).multi.ref.fits` and `finalphot(FILTER).am`. 
2. Finds number of unique stars and saves each filter's magnitude and error to the same star. Saves this to `finalmags.txt`.
3. Creates color plot with all of the images overlaid (red = J/H, green = z/y, blue = r/i) and creates plot of each filter field with green circles around each object. Saves as `coadd*(FILTER).png`. 
4. Calls `printhtml` which create HTML to view all information. 

*Supports multiple filters.*

####Params
* `prefchar` **{string}**: Filename prefix of output files from reduction module.

####Output
* `finalmags.txt` - text file with all the magnitudes for each filter of each source.
* `color.png` - overlay of all filters (red = J/H, green = z/y, blue = r/i).
* `coadd*(FILTER).png` - images of each filter field with green circles over source.
* `ratir.html`- html showing all information.

###Dependencies
[printhtml.py](#printhtml)

##printhtml
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/photometry/printhtml.py "Source") **printhtml**(<i>filters, colnames</i>)

Create `photcomp.png` comparing magnitude and errors as well as create HTML page that has all of the data easily displayed for up to 9 filters.

####Params
* `filters` **{list of str}**: List of filter types (e.g. `['r','i','z','y','J','H']`).
* `colnames` **{list of str}**: List columns names to use in table header.

####Output
* `photcomp.png` - shows magnitude vs. error for each filter
* `photom.html` - html page showing information about sources

####Dependencies
Images and `finalmags.txt` file created from [plotphotom](#plotphotom).

