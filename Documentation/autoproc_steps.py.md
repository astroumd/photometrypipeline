Reduction steps called by [autoproc.py](https://github.com/astroumd/photometrypipeline/wiki/autoproc.py)

## Pipeline Input Parameters
The following structure is used to pass commonly used pipeline variables throughout each steps. See autoproc.py at [line 153](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc.py#L153) and [autopipedefaults](#autopipedefaults) for usage and definition.

```python
pipevar = { 'autoastrocommand':'autoastrometry', 
            'getsedcommand':'get_SEDs', 
            'sexcommand':'sex', 
            'swarpcommand':'swarp', 
            'rmifiles':0, 
            'prefix':'', 
            'datadir':'' , 
            'imworkingdir':'', 
            'overwrite':0 , 
            'verbose':1,
            'flatfail':'', 
            'fullastrofail':'',
            'pipeautopath':'',
            'refdatapath':'', 
            'defaultspath':'' 
          }
```

## Functions
* [autopipedefaults](#autopipedefaults) - set commonly used variables to use throughout each step.
* [autopipeprepare](#autopipeprepare) - update image headers and performs bias/dark subtraction.
* [autopipeimflatten](#autopipeimflatten) - flatten data using flat with matching filter name.
* [autopipemakesky](#autopipemakesky) - combine sky flats based on filter type (sigma clipping for sources).
* [autopipeskysub](#autopipeskysub) - subtract both master sky and median.
* [autopipeskysubmed](#autopipeskysubmed) - subtract median, does NOT use master sky.
* [autopipecrcleanim](#autopipecrcleanim) - remove cosmic rays.
* [autopipeastrometry](#autopipeastrometry) - calculate astrometry of image files to fix WCS coordinates.
* [autopipestack](#autopipestack) - create flux scale and stack images with SWarp.

Functions are called in the same order as above, see autoproc.py at [line 240](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L240) (snippet below).

```python
for step in steps:        
  if step == 'prepare': ap.autopipeprepare(pipevar=pipevar)
  if step == 'flatten': ap.autopipeimflatten(pipevar=pipevar)
  if step == 'makesky' and nomastersky == False: ap.autopipemakesky(pipevar=pipevar)
  if step == 'skysub' and nomastersky == False:  ap.autopipeskysub(pipevar=pipevar)
  if step == 'skysub' and nomastersky == True:  ap.autopipeskysubmed(pipevar=pipevar)    
  if step == 'crclean' and nocrclean == False: ap.autopipecrcleanim(pipevar=pipevar)
  if step == 'astrometry': ap.autopipeastrometry(pipevar=pipevar),
  if step == 'stack'     : ap.autopipestack(pipevar=pipevar, customcat=customcat, customcatfilt=customcatfilt)
```

##autopipedefaults
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L18 "Source") **autopipedefaults**(<i>pipevar</i>)

Sets commonly used variables for pipeautoproc to use throughout each step. Uses [pipeautoproc.par](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/pipeautoproc.par) to set variables, otherwise set to default values. Saves in a dictionary.

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

##autopipeprepare
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L62 "Source") **autopipeprepare**(<i>pipevar</i>)

Runs [pipeprepare](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#piprepare) on every valid file and saves files with prefix `p`. Changes header with more manageable keywords and does bias/dark subtraction if bias/dark master exists (compares header keywords in files and bias/dark master).

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Saved in `pipevar['datadir']` with prefix defined in `pipevar['prefix']`.

####Output files

Saved to `pipevar['imworkingdir']` with `p` prefix.

####Dependencies

[autoproc_depend.pipeprepare](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#piprepare) 

##autopipeimflatten
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L142 "Source") **autopipeimflatten**(<i>pipevar</i>)

Flattens data using flat with matching filter name, and Saves files with prefix `f`. Flat files that don't match any filter are stored in `pipevar['flatfail']`.

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

* `prepared` files: In `pipevar['imworkingdir']` with `p` prefix.
* `flat` files: In `pipevar['imworkingdir']` containing `flat` in filename.

####Output files

Saved to `pipevar['imworkingdir']` with `f` prefix.

####Dependencies

[autoproc_depend.flatpipeproc](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#flatpipeproc).

##autopipemakesky
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L214 "Source") **autopipemakesky**(<i>pipevar</i>)

Combine sky flats based on filter type (sigma clipping for sources, at least 2 files per filter). See [skypipecombine](https://github.com/maxperry/photometrypipeline/wiki/autoproc_depend.py#skypipecombine).

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Stored in `pipevar['imworkingdir']` with `fp` prefix.

####Output files

Saved to `pipevar['imworkingdir']` with `sky-` prefix.

##autopipeskysub
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L351 "Source") **autopipeskysub**(<i>pipevar</i>)

Subtracts master sky flat from data and subtracts median.

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

* `flattened` files: Stored in `pipevar['imworkingdir']` with `fp` prefix.
* `sky` files: Stored in `pipevar['imworkingdir']` containing `sky-` in filename.

####Output files

Saved to `pipevar['imworkingdir']` with `s` prefix.

####Dependencies

[autoproc_depend.skypipeproc](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#skypipeproc).

##autopipeskysubmed
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L291 "Source") **autopipeskysubmed**(<i>pipevar</i>)

Subtracts median ONLY, does NOT use master sky.

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Stored in `pipevar['imworkingdir']` with `fp` prefix.

####Output files

Saved to `pipevar['imworkingdir']` with `s` prefix.

##autopipecrcleanim
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L426 "Source") **autopipecrcleanim**(<i>pipevar</i>)

Removes cosmic rays 

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Stored in `pipevar['imworkingdir']` with `sfp` prefix.

####Output files

Saved to `pipevar['imworkingdir']` with `z` prefix.

####Dependencies

[autoproc_depend.cosmiczap](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#cosmiczap).

##autopipeastrometry
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L496 "Source") **autopipeastrometry**(<i>pipevar</i>)

Calculates astrometry of image files to fix WCS coordinates (shift and rotation) in header. Using fast astrometry solver ([vlt_autoastrometry.py](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/astrom/vlt_autoastrometry.py)) with pair-distance matching and asterism matching.  Returns file with corrected WCS coordinates saved with `a` prefix. Runs `Scamp` for additional astrometry corrections, twice, once for basic individual LOOSE correction, second correct all together. Uses distortion of 3 as default, but uses 7 if distortion parameters high (i.e. RATIR H2RG).

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Stored in `pipevar['imworkingdir']` with `zsfp` prefix.

####Output files

Saved to `pipevar['imworkingdir']` with `a` prefix.

####Dependencies

[autoproc_depend.astrometry](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#astrometry), [vlt_autoastrometry.py](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/astrom/vlt_autoastrometry.py), SCAMP, SExtractor.

####Usage of [vlt_autoastrometry.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/astrom/vlt_autoastrometry.py)

See autoproc_steps.py at [line 541](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L541), `pipevar['autoastrocommand']` is `vlt_autoastrometry.py` as defined in [pipeautoproc.par](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/pipeautoproc.par).

```python
head = pf.getheader(file)
        
targ = head['TARGNAME']
sat  = head['SATURATE']
        
if 'flat' in targ: continue
  cmd = 'python ' + pipevar['autoastrocommand'] + ' ' + file + ' -l ' + str(sat)
```

##autopipestack
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L635 "Source") **autopipestack**(<i>pipevar</i>)

Does zeropoint correction on each individual frame using `sextractor` and `get_SEDs`. Creates flux scale (newflxsc) from how close to median of zeropoint values. Uses flux scale to stack images in `swarp` (has moved bad zeropoint values and bad newflxsc values to marked folders - `badzptfit/` and `badflxsc/`) and calculates absolute zeropoint correction of coadd. Saves zeropoint plot as `zpt_(FILTER).ps`.

####Params

* `pipevar` **{dict}**: (Input pipeline parameters)[#pipeline-input-parameters].

####Input files

Stored in `pipevar['imworkingdir']` with `a` prefix.

####Output files

Saved to `pipevar['imworkingdir']` with `coadd` prefix as `.fits` and `.weight.fits`.

####Dependencies

SWarp, [get_SEDs.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/photometry/dependencies/get_SEDs.py), [autoproc_depend.calc_zpt](https://github.com/maxperry/photometrypipeline/wiki/autoproc_depend.py#calc_zpt), [autoproc_depend.findsexobj](https://github.com/maxperry/photometrypipeline/wiki/autoproc_depend.py#findsexobj) (SExtractor), [autoproc_depend.identify_matches](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_depend.py#L989).

####Usage of [autoproc_depend.findsexobj](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#findsexobj)

1. See autoproc_steps.py at [line 741](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L741).

 ```python
for sfile in stacklist:
  head = pf.getheader(sfile)
  ipixscl = head['PIXSCALE']
  apd.findsexobj(sfile, 3.0, pipevar,pix=ipixscl,aperture=20.0, quiet=quiet)
```

 **NOTE**: `stacklist` is a `list` of files that meet the following conditions:

 * same filter (`head['FILTER']`)
 * same target (`head['TARGNAME']`)
 * `head['ASTRRMS1']` in (2.0e-4, 5.0e-6)
 * `head['ASTRRMS2']` in (2.0e-4, 5.0e-6)

2. See autoproc_steps.py at [line 931](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L931).

 ```python
head   = pf.getheader(outfl)
pixscl = head['PIXSCALE']
try:
  apd.findsexobj(outfl, 10.0, pipevar, pix=pixscl, aperture=20.0, wtimage=outwt, quiet=quiet)
```

 **NOTE**: `outfl` is the fits file coadded with the flux scale using SWarp (see [line 903](https://github.com/maxperry/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L903)).

3. See autoproc_steps.py at [line 940](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L940).

 ```python                        
head   = pf.getheader(outfl)
cpsfdi = 1.34 * float(head['SEEPIX'])
            
# Run sextractor again on new coadd file
apd.findsexobj(outfl, 3.0, pipevar, pix=pixscl, aperture=cpsfdi, wtimage=outwt, quiet=quiet)
```

####Usage of [get_SEDs.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/photometry/dependencies/get_SEDs.py)

See autoproc_steps.py at [line 805](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L805).

```python
# If custom catalog not provided, catalog doesn't include band_filter, or 
# no objects from catalog found in image then
# use get_SEDs.py to make catalog using 2MASS + (SDSS or APASS or USNOB1)
if nocustomcat:
  # Create catalog star file 
  # (python get_SEDs.py imfile band_filter catfile USNOB_THRESH alloptstars)
  sedcmd = 'python ' + pipevar['getsedcommand'] + ' ' + imfile + ' ' +\
                         filter + ' ' + catfile + " 15 True "+ qtcmd
```

####Usage of SWarp

See autoproc_steps.py at [line 907](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L907).

```python
stackcmd = pipevar['swarpcommand']
# Keywords to carry through
stackcmd += ' -COPY_KEYWORDS OBJECT,TARGNAME,FILTER,' +\
                   'INSTRUME,PIXSCALE,WAVELENG,DATE-OBS,AIRMASS,FLATFLD,FLATTYPE '
# Create output variables that will be used by SWarp
outfl = pipevar['imworkingdir'] + 'coadd' + targ + '_'+ re.sub(r'[^\w]', '', medtime)+'_'+ filter + '.fits'
outwt = pipevar['imworkingdir'] + 'coadd' + targ + '_'+ re.sub(r'[^\w]', '', medtime)+'_'+ filter + '.weight.fits'
            
if pipevar['verbose'] > 0:
  stackcmd += ' -VERBOSE_TYPE NORMAL '
else:
  stackcmd = stackcmd + ' -VERBOSE_TYPE QUIET '
            
# Coadd with flux scale
stackcmd = stackcmd + ' -SUBTRACT_BACK N -WRITE_XML N -IMAGEOUT_NAME ' +\
                       outfl + ' -WEIGHTOUT_NAME ' + outwt +\
                       ' -FSCALE_KEYWORD NEWFLXSC ' + newtextslist
```

####Usage of [autoproc_depend.calc_zpt](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#calc_zpt)

1. See autoproc_steps.py at [line 844](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L844).

 ```python
zpt, scats, rmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]), np.array([obswts]), sigma=3.0)
```

2. See autoproc_steps.py at [line 1039](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/auto/autoproc_steps.py#L1039).

 ```python
czpts, cscats, crmss = apd.calc_zpt(np.array([refmag]), np.array([obskpm]), 
                                    np.array([obswts]), sigma=1.0,
                                    plotter=pipevar['imworkingdir']+'zpt_'+filter+'.ps')
```