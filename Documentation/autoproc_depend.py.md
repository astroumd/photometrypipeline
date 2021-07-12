## Functions
* [pipeprepare](#piprepare) - Normalize header keywords in FITS files.
* [flatpipeproc](#flatpipeproc) - Check if flat is same size as data, then divide for correct filter.
* [skypipecombine](#skypipecombine) - Create sigma clipped median sky flat.
* [skypipeproc](#skypipeproc) - Subtract sky flat from data, then subtract median of that from remaining data. 
* [cosmiczap](#cosmiczap) - Remove cosmic rays using Laplacian cosmic ray identification written in `cosmics.py`.
* [astrometry](#astrometry) - Run `sextractor` and `scamp` to refine astrometric solution.
* [findsexobj](#findsexobj) - Find `sextractor` objects with optional inputs. Estimates seeing from stars found.  
* [calc_zpt](#calc_zpt) - Find zeropoint using robust scatter.
* [robust_scat](#robust_scat) - Calculate robust scatter and set the weight of those above this limit to 0.
* [medclip](#medclip) - Median iterative sigma-clipping.
* [medclip2d](#medclip2d) - Median iterative sigma-clipping over 2d array.
* [identify_matches](#identify_matches) - Use a kd-tree (3d) to match two lists of stars, using full spherical coordinate distances.

##piprepare
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L12 "Source") **piprepare**(<i>filename, [outname=None, biasfile=None, darkfile=None, verbose=1]</i>)

Adds additional header keywords needed later on in the pipeline and removes unnecessary header keywords by looking through a list of mandatory keywords. Also runs bias and dark subtraction for filters with an existing master bias/dark (CCDs). The prepared images are written to disk with `outname`.

####Params

* `filename` **{str or list of str}**: Absolute path to FITS file, or list of paths, or file w/list of paths.
* `outname` **{str, optional}**: Specify output file to write to disk.
* `biasfile` **{str, optional}**: Absolute path to bias file.
* `darkfile` **{str, optional}**: Absolute path to dark file.
* `verbose` **{bool, optional}**: Print out comments.

##flatpipeproc
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L131 "Source") **flatpipeproc**(<i>filename, flatname, [flatminval=0, flatmaxval=0]</i>)

Checks if flat is same size as data, then divides for correct filter.

####Params

* `filename` **{str or list of str}**: Absolute path to FITS file, or array of paths, or file w/list of paths.
* `flatname` **{str}**: Absolute path to FITS master flat file.
* `flatminval` **{float, optional}**: If not set to 0 below this value will set to NaNs.
* `flatmaxval` **{float, optional}**: If not set to 0 below this value will set to NaNs.

##skypipecombine
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L220 "Source") **skypipecombine**(<i>filelist, outfile, filt, pipevar, [removeobjects=None, objthresh=6, algorithm='median', trimlo=None, trimhi=None, mincounts=1, maxcounts=55000, satlevel=30000, type=None]</i>)

Create sigma clipped median sky flat. Scales each file based on the overall sigma clipped median, then removes objects selected with `sextractor` (uses flux fraction radius) in each file. Removes saturated pixels. Calculates sigma clipped median of each pixel and saves anything with non-finite values (saturated or source) to the median of the entire frame.  Save with `outfile` name.

####Params

* `filename` **{str or list of str}**: Absolute path to FITS file, or list of paths, or file w/list of paths.
* `outfile` **{str}**: Name for output fits file.
* `filt` **{str}**: Absolute path to filter file.
* `pipevar` **{dict}**: (Input pipeline parameters)[https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#pipeline-input-parameters].
* `removeobjects` **{bool, optional}**: Specifies if objects should be removed.
* `objthresh` **{float, optional}**: Sets sigma in [removeobjects](#removeobjects) (default is 6).
* `algorithm` **{str, optional}**: Algorithm to solve (mean or median, default is `median`).
* `trimlo` **{bool, optional}**: Trim off bottom of data in mean algorithm mode (default is 25%).
* `trimhi` **{bool, optional}**: Trim off top of data in mean algorithm mode (default is 25%).
* `mincounts` **{int, optional}**: Sets minimum counts allowed (default is 1).
* `maxcounts` **{int, optional}**: Sets minimum counts allowed (default is 55000).
* `satlevel` **{int, optional}**: Sets saturation level (default is 30000).
* `type` **{str, optional}**: Sets `SKYTYPE` keyword in header of `outfile` to this string.

#### Dependencies
[medclip](#medclip), `sextractor`.

#### Example & Usage
See autoproc_steps.py at [line 277](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L277).

```python
skypipecombine( files[skyflats], 
                outflatname, 
                file,
                pipevar, 
                removeobjects=True, 
                type='sky' )
```

##skypipeproc
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L428 "Source") **skypipeproc**(<i>filename, flatname, outfile, [flatminval=None, flatmaxval=None]</i>)

Subtracts sky flat from data and then subtracts median of that from remaining data. 

####Params

* `filename` **{str or list of str}**: Absolute path to FITS file, or list of paths to be sky subtracted.
* `flatname` **{str}**: Absolute path to sky flat FITS file.
* `outfile` **{str}**: Name for output fits file.
* `flatminval` **{float, optional}**: Minimum required value in flat (default for `SKYCTS` calculation is 0.1).
* `flatmaxval` **{float, optional}**: Maximum required value in flat.

##cosmiczap
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L534 "Source") **cosmiczap**(<i>filename, outname, [sigclip=6.0, maxiter=3, verbose=True]</i>)

Removes cosmic rays using Laplacian cosmic ray identification written in [cosmics.py](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/cosmics.py).  

####Params

* `filename` **{str or list of str}**: Absolute path to FITS file, or list of paths to be cosmic ray zapped.
* `outfile` **{str}**: Name for output fits file.
* `sigclip` **{float, optional}**: Sigma to clip.
* `maxiter` **{int, optional}**: Maximum number of times to iterate loop.
* `verbose` **{bool, optional}**: Print out comments.

#### Example & Usage
See autoproc_steps.py at [line 484](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L484).

```python
cosmiczap( file, 
           outfile, 
           sigclip=6.0, 
           maxiter=1, 
           verbose=pipevar['verbose'] )
```

##astrometry
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L573 "Source") **astrometry**(<i>atfimages, [scamprun=1, pipevar=None]</i>)

Run `sextractor` and `scamp` to refine astrometric solution.  

####Params

* `atfimages` **{list of str}**: List of absolute paths of images to run through `scamp`.
* `scamprun` **{int, optional}**: The first run does a LOOSE run with distortion degree 1, any other run will look for high distortion parameters, if it finds it will use distortion degree 7, otherwise 3 (will also cut out `FLXSCALE` on runs after 1).
* `pipevar` **{dict}**: (Input pipeline parameters)[https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#pipeline-input-parameters].

#### Example & Usage
See autoproc_steps.py at [line 619](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L619) and [line 622](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L622).

##findsexobj
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L673 "Source") **findsexobj**(<i>file, sigma, pipevar, [masksfx=None, zeropt=25.0, maptype='MAP_WEIGHT', wtimage=None, fwhm=1.5, pix=0.3787, aperture=5.0, elong_cut=1.5, quite=0]</i>)

Finds `sextractor` objects with optional inputs. Estimates seeing from stars found.   

####Params

* `file` **{str}**: Absolute path to FITS file to run `sextractor` on.
* `sigma` **{float}**: Detection threshold and analysis threshold for `sextractor`.
* `pipevar` **{dict}**: (Input pipeline parameters)[https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#pipeline-input-parameters].
* `masksfx` **{str, optional}**: Identifier for `sextractor` `CHECKIMAGE_NAME`.
* `zeropt` **{float, optional}**: Input value for `sextractor` `MAG_ZEROPOINT`.
* `wtimage` **{str, optional}**: Absolute file for input for `sextractor` `WEIGHT_IMAGE`.
* `fwhm` **{float, optional}**: Input value for `sextractor` `SEEING_FWHM`.
* `pix` **{float, optional}**: Input value for `sextractor` `PIXEL_SCALE`.
* `aperture` **{float, optional}**: Input value for `sextractor` `PHOT_APERTURES`.
* `elong_cut` **{float, optional}**: Cutoff limit for `FWHM` calculation of elongation to eliminate non-stars.
* `quiet` **{bool, optional}**: Disable output from `sextractor`.

#### Example & Usage
See autoproc_steps.py (`autopipestack`) at [line 745](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L745), [line 935](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L935) and  [line 944](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L944).

##calc_zpt
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L785 "Source") **calc_zpt**(<i>catmag, obsmag, wts, [sigma=3.0, plotter=None]</i>)

Find zeropoint using robust scatter.  

####Params

* `catmag` **{list}**: 2d array with catalog magnitudes `catmag[nobs,nstar]`.
* `obsmag` **{list}**: 2d array with catalog magnitudes `obsmag[nobs,nstar]`.
* `wts` **{list}**: 2d array with weights `wts[nobs,nstar]`.
* `sigma` **{float, optional}**: Sigma value for how far values can be from robust scatter value.
* `plotter` **{str, optional}**: Absolute path to save zeropoint plot.

####Output
* `z2` - Zeropoint correction.
* `scats` - Robust scatter of each observation.
* `rmss` - Standard deviation (without bad weight points) of each observation.

#### Example & Usage
See autoproc_steps.py (`autopipestack`) at [line 844](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L844), and [line 1039](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L1039).

##robust_scat
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L856 "Source") **robust_scat**(<i>diff, wts, nobs, nstars, sigma</i>)

Calculate robust scatter and set the weight of those above this limit to 0.  

####Params

* `diff` **{list}**: Values to calculate robust scatter over.
* `wts` **{list}**: Weights (0 is bad).
* `nobs` **{int}**: Number of observations to iterate over.
* `nstars` **{int, optional}**: Number of stars to iterate over.
* `sigma` **{float, optional}**: Sigma*robust scatter that is acceptable.

####Output
* `scats` - robust scatter of each observation.
* `rmss` - standard deviation (without bad weight points) of each observation.

#### Example & Usage
See autoproc_depend.py (`calc_zpt`) at [line 829](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L829).

##medclip
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L892 "Source") **medclip**(<i>indata, [clipsig=3.0, maxiter=5, verbose=0]</i>)

Median iterative sigma-clipping.

####Params

* `indata` **{list}**:  Array to be clipped.
* `clipsig` **{float, optional}**: Sigma to clip around.
* `maxiter` **{int, optional}**: Maximum number of times to clip.
* `verbose` **{bool, optional}**: Print out comments.

#### Example & Usage
See autoproc_depend.py (`skypipecombine`) at [line 305](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L305), [line 352](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L352).

##medclip2d
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L892 "Source") **medclip**(<i>indata, [clipsig=3.0, maxiter=5, verbose=0, overaxis=0]</i>)

Median iterative sigma-clipping.

####Params

* `indata` **{list}**:  Array to be clipped.
* `clipsig` **{float, optional}**: Sigma to clip around.
* `maxiter` **{int, optional}**: Maximum number of times to clip.
* `verbose` **{bool, optional}**: Print out comments.
* `overaxis` **{int, optional}**: Axis that we want to take median over.

#### Example & Usage
See autoproc_depend.py (`skypipecombine`) at [line 377](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L377).

##identify_matches
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L990 "Source") **identify_matches**(<i>queried_stars, found_stars, [match_radius=3.0]</i>)

Use a kd-tree (3d) to match two lists of stars, using full spherical coordinate distances.

####Params

* `queried_stars` **{list}**: Numpy arrays of `[ [ra,dec],[ra,dec], ... ]` (all in decimal degrees).
* `found_stars` **{list}**: Numpy arrays of `[ [ra,dec],[ra,dec], ... ]` (all in decimal degrees).
* `match_radius` **{float, optional}**: Max distance (in arcseconds) allowed to identify a match between two stars.

#### Example & Usage
See autoproc_steps.py (`autopipestack`) at [line 790](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L790), and [line 988](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L988).