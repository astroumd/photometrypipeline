Below is a list of all functions that are part of the reduction module, grouped by filename.

* [preproc.py](https://github.com/astroumd/photometrypipeline/wiki/preproc.py) -  **preprocessing functions**.
  * [choose_calib](https://github.com/astroumd/photometrypipeline/wiki/preproc.py#choose_calib) - calibration images visualization and selection.
  * [choose_science](https://github.com/astroumd/photometrypipeline/wiki/preproc.py#choose_science) - display science images for verification by user.
  * [mkmaster](https://github.com/astroumd/photometrypipeline/wiki/preproc.py#mkmaster) - make master calibration frames (bias, dark, flat).

* [autoproc.py](https://github.com/astroumd/photometrypipeline/wiki/autoproc.py) - **automated pipeline script**.
  * [autoproc](https://github.com/astroumd/photometrypipeline/wiki/autoproc.py#autoproc) - main function.
* [autoproc_steps.py](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py) - **reduction processing steps**.

  * [autopipedefaults](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipedefaults) - set commonly used variables to use throughout each step.
  * [autopipeprepare](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipeprepare) - update image headers and performs bias/dark subtraction.
  * [autopipeimflatten](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipeimflatten) - flatten data using flat with matching filter name.
  * [autopipemakesky](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipemakesky) - combine sky flats based on filter type (sigma clipping for sources).
  * [autopipeskysub](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipeskysub) - subtract both master sky and median.
  * [autopipeskysubmed](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipeskysubmed) - subtract median, does NOT use master sky.
  * [autopipecrcleanim](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipecrcleanim) - remove cosmic rays.
  * [autopipeastrometry](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipeastrometry) - calculate astrometry of image files to fix WCS coordinates.
  * [autopipestack](https://github.com/astroumd/photometrypipeline/wiki/autoproc_steps.py#autopipestack) - create flux scale and stack images with SWarp.

* [autoproc_depend.py](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py) - **reduction processing dependencies**.
  * [pipeprepare](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#piprepare) - Normalize header keywords in FITS files.
  * [flatpipeproc](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#flatpipeproc) - Check if flat is same size as data, then divide for correct filter.
  * [skypipecombine](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#skypipecombine) - Create sigma clipped median sky flat.
  * [skypipeproc](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#skypipeproc) - Subtract sky flat from data, then subtract median of that from remaining data. 
  * [cosmiczap](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#cosmiczap) - Remove cosmic rays using Laplacian cosmic ray identification written in `cosmics.py`.
  * [astrometry](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#astrometry) - Run `sextractor` and `scamp` to refine astrometric solution.
  * [findsexobj](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#findsexobj) - Find `sextractor` objects with optional inputs. Estimates seeing from stars found.  
  * [calc_zpt](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#calc_zpt) - Find zeropoint using robust scatter.
  * [robust_scat](#robust_scat) - Calculate robust scatter and set the weight of those above this limit to 0.
  * [medclip](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#medclip) - Median iterative sigma-clipping.
  * [medclip2d](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#medclip2d) - Median iterative sigma-clipping over 2d array.
  * [identify_matches](https://github.com/astroumd/photometrypipeline/wiki/autoproc_depend.py#identify_matches) - Use a kd-tree (3d) to match two lists of stars, using full spherical coordinate distances.