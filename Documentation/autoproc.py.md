Automated image reduction, performed in a sequence of [steps](#steps).

## Functions

* [autoproc](#autoproc) - main function.

##autoproc
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc.py#L4 "Source") **autoproc**(<i>datadir='./', imdir='./imredux/', [start, stop, only, step, redo, nocrclean, nomastersky, quiet, rmifiles, customcat]</i>)

Fully-automated image reduction. Returns processed data in `imdir` folder.

####Params

* `datadir` **{str}**: Location of raw data (`current directory` if unspecified), followed by a slash.
* `imdir` **{str}**: Location of processed data (`imredux` if unspecified), followed by a slash.
* `start` **{str, optional}**: Start with this step, skipping previous ones.
* `stop` **{str, optional}**: End with this step, skipping subsequent ones.
* `only` **{str, optional}**: Do only this step.
* `step` **{str, optional}**: Completely identical to `only`, takes precedence.
* `nocrclean` **{bool, optional}**: Do not zap cosmic rays.
* `nomastersky` **{bool, optional}**: Do not create master sky, only subtract median of sky.
* `redo` **{bool, optional}**: Repeat step(s), overwriting any existing files.
* `quiet` **{bool, optional}**: (mainly) silent output unless errors.
* `rmifiles` **{bool, optional}**: Removes intermediate files.
* `customcat` **{str, optional}**: Custom catalog (txt file) to determine instrumental zeropoint corrections. Must be in same format as what `get_SEDs.py` produces:

	```
	ra(deg)	dec(deg)	u	g	r	i	z	y	B	V	R
	I	J	H	K	u_err	g_err	r_err	i_err y_err B_err
	V_err	R_err	I_err	J_err	H_err	K_err Mode
	```
First line will be skipped (use for headers) everything but JHK are expected to be in AB magnitudes, JHK should be in same units as 2MASS (Vega mags).
* `customcatfilt` **{list of str, optional}**: Filters relevant to custom catalog file (all other filters will use `get_SEDs.py` to calculate catalog from 2MASS + (SDSS or APASS or USNOB1) in that order.

**Params for `start`, `stop`, `only`, `step`** 
 listed in order of execution: `'prepare'`, `'flatten'`, `'makesky'`, `'skysub'`, `'crclean'`, `'astrometry'`, `'stack'`.

####Example

```python
autoproc(datadir='./raw/', imdir='./reduced/', redo=1)
```

####Additional Options

If any of the following files are found in the directory where autoproc is run, it will change the default behavior.
       
* [pipeautoproc.par](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/pipeautoproc.par): Contains various defaults (mainly directory paths).

####Comments

This code is meant to be fully automated, producing science-quality
output with a single command line with zero intervention.  Reduces
RATIR and LMI data.

The production of images is quite high-level and includes photometric 
calibration of the field (although the accuracy of this has not been
robustly tested).

The program works in a series of steps following standard CCD reduction
techniques, automatically recognizing calibration files, matching
calibrations together and science frames with appropriate calibrations,
etc.  Users looking for more control over the reductions can run each
step individually with the `step` command.

If the reduction is interrupted, the program will resume without
incident (although any failed steps may be repeated); each task is 
run independently of all others.
   
The program tries very hard to correct for observer mistakes. But it's not perfect.
If problems occur, generally the easiest fix is to delete any 
offending files and rerun. More significant problems generally require direct 
modification of the code. 

Filenames for the input raw images are expected to be in 
`2*.fits` format. They can be either in the working directory or 
in a different directory specified by datadir.

####Steps
The code runs in a series of steps in the following order:

1. **Prepare** ([autopipeprepare](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L62 "Source"))
 the data by converting the multi-header extension FITS files
to a standard single frame, and adding extra information to the header. 

 Output: `p*.fits` (written to `imredux/` by default.)

2. **Flatten** ([autopipeimflatten](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L142 "Source")) data.  Divide each image by the flatfield.  A more
refined cropping is also done at this stage, depending on the placement
of the image area during the run (which is variable.) 

 Assumes master flat (eg. `flat_H.fits` in `imredux/` folder.)

3. **Makesky** ([autopipemakesky](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L214 "Source")) makes a master sky (eg. `sky-H.fits` in `imredux/` folder.)

4. **Skysubtract** ([autopipeskysub](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L351 "Source"), or [autopipeskysubmed](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L291 "Source") if `nomastersky=True`) subtracts out master sky.

 Output: `sp*.fits` (written to `imredux/` by default.)

5. **Remove cosmic rays** ([autopipecrcleanim](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L426 "Source")) using the independent routines [cosmics.py](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/cosmics.py) that uses Laplacian cosmic ray identification. See the source for more information.  

 *This can be a time-consuming process*. 
 
  Output: `zfp*.fits`

6. **Solve astrometry** ([autopipeastrometry](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L496 "Source")) of the field against the best available online catalog
(SDSS/2MASS/APASS/USNO-B1.0), using the independent [vlt_autoastrometry.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/astrom/vlt_autoastrometry.py) code.
(Also requires `sextractor` to be installed.) 

 Uses two passes of `Scamp` for a 
secondary correction.  Scamp accounts for distortion.  Uses higher distortion
parameters if already supplied distortion keywords.

 Output: `azfp*.fits`

7. **Stack exposures** ([autopipestack](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L635 "Source"))  A weighted, masked median is performed for each field. Requires `swarp` to be installed and properly linked.

 Output: `coadd[object].[filter].fits`

#### Troubleshooting
If the pipeline crashes during operation:
    
* ***Did an external call (autoastrometry, swarp, sex) fail?***
    Check that you have specified paths correctly and have the required 
    software installed for astrometry/coadding
    
* ***Did it encounter some other error while processing a file?***  
    Check which file the program was working on when it crashed. If the
    file is not essential, try deleting it and re-running the pipeline
    starting with the current step (or delete ALL of the file's precursors
    with the same file number and rerun the pipeline.) 
    
If processing completed, but the results are problematic:

* ***Did it return without actually processing any data?*** 
    Make sure that it is in the current working directory or that you 
    have correctly pointed to the directory containing raw data with the 
    "datadir" keyword. If you are re-doing a step, it will not overwrite 
    existing files by default. Set the redo keyword to overwrite files 
    from a previously-attempted step (be sure to set "start" or "step" 
    unless you want to restart the pipeline from the beginning.)

* ***Were some files skipped?***
    If the pipeline encountered a non-fatal problem processing an
    individual image (such as an inability to flatfield) then it will 
    not process that file any further. For most cases a summary of the 
    problems will be printed out at the end of processing. If a file is 
    not being processed and you do not see it in the final summary,
    you can simply rerun the pipeline (without deleting any files and 
    without setting the redo flag) and it will try to repeat any failed 
    steps of this nature. 
    