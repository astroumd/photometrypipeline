Collection of pre-processing functions for general use.

## Functions
* [choose_calib](#choose_calib) - calibration images visualization and selection.
* [choose_science](#choose_science) - display science images for verification by user.
* [mkmaster](#mkmaster) - make master calibration frames (bias, dark, flat).

##choose_calib
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/preproc.py#L26 "Source") **choose\_calib**(<i>instrument, ftype, [workdir='./', cams=[0,1,2,3], auto=False, reject\_sat=True, amin=0.2, amax=0.8, save\_select=True, figsize=(8,5), noplot=False]</i>)

Either auto-select or display calibration images for user verification.

####Params

* `instrument` **{str}**: Instrument name defined in [specific\_instruments.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/instruments/specific_instruments.py) (e.g. `ratir`)
* `ftype` **{str}**: Type of calibration frames (e.g. `bias`, `dark`, `flat`).
* `workdir` **{str, optional}**: Location of calibration frames (current directory if unspecified), followed by a slash.
* `cams` **{list of str, optional}**: Camera numbers, default is all. (e.g. `[0,1,2,3]`)
* `auto` **{bool, optional}**: Automated frame selection. If `bias`, will select all, if `flat` will select non-saturated frames with sufficient counts.
* `reject_sat` **{bool, optional}**: Reject frames with saturated pixels.
* `amin` **{float, optional}**: Minimum fraction of saturation value for median (automated).
* `amax` **{float, optional}**: Maximum fraction of saturation value for median (automated).
* `save_select` **{bool, optional}**: Save dictionary of selected frames to python pickle file.
* `figsize` **{tuple of int}**: Dimensions of figure used to display frames for selection.
* `noplot` **{bool, optional}**: Don't display calibration frames.

####Example

```python
dict = choose_calib( instrument='ratir', 
                     ftype='bias',
                     workdir='./fits/bias/',
                     cams =[0,1,2,3] )
```
**Note**: Call [mkmaster](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/preproc.py#L813) using resulting dictionary, or path to pickle file (if `save_select=True`).  


##choose_science 
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/preproc.py#L473 "Source") **choose\_science**(<i>instrument, [workdir='./', targetdir='./', cams=[0,1,2,3], auto=False, save\_select=True, figsize=(10,10), calibrate=False, noplot=False]</i>)

Display science images for verification by user.

####Params

* `instrument` **{str}**: Instrument name defined in [specific\_instruments.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/instruments/specific_instruments.py) (e.g. `ratir`)
* `workdir` **{str, optional}**: Location of all science frames (`current directory` if unspecified), followed by a slash.
* `targetdir` **{str, optional}**: Location of selected science frames (`current directory` if unspecified), followed by a slash.
* `cams` **{list of str, optional}**: Camera numbers, default is all. (e.g. `[0,1,2,3]`)
* `auto` **{str, optional}**: Select all science frames.
* `save_select` **{bool, optional}**: Save dictionary of selected frames to python pickle file.
* `figsize` **{tuple of int}**: Dimensions of figure used to display frames for selection.
* `calibrate` **{bool, optional}**: Subtract master bias, master dark to science frames, and divide by master flat before displaying.
* `noplot` **{bool, optional}**: Don't display science frames.

####Example

```python
dict = choose_science( instrument='ratir', 
                       workdir='./fits/science',
                       targetdir='./fits/science/selected', 
                       cams=[0,1,2,3], 
                       calibrate=True )
```

##mkmaster 
[<>](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/preproc.py#L813 "Source") **mkmaster**(<i>instrument, fn_dict, mtype, [fmin=5]</i>)

Make master calibration frames (bias, dark, flat). Currently no outlier rejection other than median combine.
        
####Params

* `instrument` **{str}**: Instrument name defined in [specific\_instruments.py](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/instruments/specific_instruments.py) (e.g. `ratir`)
* `fn_dict` **{dict or str}**: Dictionary output by [choose_calib](https://github.com/astroumd/photometrypipeline/blob/master/photopipe/reduction/preproc.py#L26) containing organized FITS file names. Can also provide filename of pickled dictionary.
* `mtype` **{str}**: Type of master frame. Should be either `bias`, `dark` or `flat`.
* `fmin` **{str, optional}**: Minimum number of files needed to make a master frame.

####Example

```python
mkmaster( instrument='ratir', 
          fn_dict='./output_from_choose_calib'
          mtype='bias' )
```