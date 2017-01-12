## LaSRC Version 1.0.0 Release Notes
Release Date: {sometime} 2017

### Downloads
LaSRC (Landsat Surface Reflectance Code) source code

    git clone https://github.com/USGS-EROS/espa-surface-reflectance.git

LaSRC auxiliary files

    http://edclpdsftp.cr.usgs.gov/downloads/auxiliaries/l8sr_auxiliary/l8sr_auxiliary.tar.gz

See git tag [lasrc-version_1.0.0]

### Installation
  * Install dependent libraries - ESPA product formatter (https://github.com/USGS-EROS/espa-product-formatter)
  * Set up environment variables.  Can create an environment shell file or add the following to your bash shell.  For C shell, use 'setenv VAR "directory"'.
```
    export PREFIX="path_to_directory_for_lasrc_build_data"
```

  * Install baseline auxiliary files and set up the environment variables.
```
    tar -xvzf l8sr_auxiliary.tar.gz
    export L8_AUX_DIR="directory_saved_auxiliary_files"
    (or in c shell use 
    setenv L8_AUX_DIR "directory_saved_auxiliary_files")
```

  * Download (from Github USGS-EROS surface-reflectance project) and install source files. The following build will create a list of executable files under $PREFIX/bin (tested in Linux with the gcc compiler). It will also copy various scripts from the scripts directory up the the $PREFIX/bin directory.
```
    cd lasrc\c_version\src
    make
    make install
    make clean

    cd lasrc\landsat_aux\src
    make
    make install
    make clean
```

  * Test - Download Landsat Level 1 files.  Run the do\_lasrc Python script in the PREFIX bin directory to run the applications.  Use do\_lasrc.py --help for the usage information.  This script requires that your LaSRC binaries are in your $PATH or that you have a $BIN environment variable set up to point to the PREFIX bin directory.
```
    convert_lpgs_to_espa --mtl <Landsat_MTL_file>
    do_lasrc.py --xml <Landsat_ESPA_XML_file>
```

  * Check output
```
    {scene_name}_toa_*: top-of-atmosphere (TOA) reflectance (or brightness temperature for thermal bands) in internal ESPA file format
    {scene_name}_sr_*: surface reflectance in internal ESPA file format
```

  * Note that the FORTRAN version contains the code as delivered from Eric Vermote and team at NASA Goddard Space Flight Center.  The C version contains the converted FORTRAN code into C to work in the ESPA environment.  It also contains any bug fixes, agreed upon by Eric's team, along with performance enhancements.  The FORTRAN code contains debugging and validation code, which is not needed for production processing.

### Dependencies
  * ESPA raw binary and ESPA common libraries from ESPA product formatter and associated dependencies
  * XML2 library
  * Auxiliary data products
    1. LAADS Terra and Aqua CMG and CMA data
    2. CMGDEM HDF file
    3. Various input files and model information provided with the LaSRC auxiliary .tar.gz file

### Auxiliary Data Updates
The baseline auxiliary files provided don't include the daily climate data.  In order to generate or update the auxiliary files to the most recent day of year (actually the most current auxiliary files available will be 2-3 days prior to the current day of year do to the latency of the underlying LAADS products) the user will want to run the updatelads.py script available in $PREFIX/bin.  This script can be run with the "--help" argument to print the usage information.  In general the --quarterly argument will reprocess/update all the LAADS data back to 2013.  This is good to do every once in a while to make sure any updates to the LAADS data products are captured.  The --today command-line argument will process the LAADS data for the most recent year.  In general, it is suggested to run the script with --quarterly once a quarter.  Then run the script with --today on a nightly basis.  This should provide an up-to-date version of the auxiliary input data for LaSRC.  The easiest way to accomplish this is to set up a nightly and quarterly cron job.

The updatelads script requires a username/password to access the ladssci.nascom.nasa.gov FTP site.  The user will need to contact USGS EROS Customer Services to obtain a username/password for the LAADS FTP site.  In your email explain that you will be using this ftp access to obtain LAADS data for processing Landsat 8 products using the LaSRC application provided by the USGS EROS.  For questions regarding this information, please contact the Landsat Contact Us page and specify USGS CDR/ECV in the "Regarding" section. https://landsat.usgs.gov/contactus.php

The provided username and password should be used in the --username and --password command-line arguments for the updatelads.py script.  If not specified the source code will try to use the ESPA_LAADS_CONFIG http service to automatically determine the username/password, which is only available to the USGS LSRD systems.

### Data Preprocessing
This version of the LaSRC application requires the input Landsat products to be in the ESPA internal file format.  After compiling the product formatter raw\_binary libraries and tools, the convert\_lpgs\_to\_espa command-line tool can be used to create the ESPA internal file format for input to the LaSRC application.

### Data Postprocessing
After compiling the product-formatter raw\_binary libraries and tools, the convert\_espa\_to\_gtif and convert\_espa\_to\_hdf command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Verification Data

### User Manual

### Product Guide

## Release Notes
Changes to the C version were made based on the new version 3.5.2 of the
FORTRAN code, which is also included in this release.  Changes to that FORTRAN
code include the following:

* subaeroret.f (and new subaeroretwat.f)
1. Modifications to how aerosol retrieval is handled.  Modified the algorithm
   and added a new algorithm for retrieval over water.  This routine is very
   similar to the non-water retrieval, but the residual computations are
   different for the water pixels.  There are also a couple of other minor
   differences in the two aerosol retrieval functions.  (subaeroret.c contains
   subaeroret and subaeroretwat)

* LUTldcm_subr_v4.f
2. The atmospheric correction now uses the angstrom coefficient to modify the
   passed in raot550nm.  This modified value is used in the atmospheric
   corrections versus the original value.  (Affects atmcorlamb2 in lut_subr.c)

* possolp.f
3. Added a routine to compute the solar azimuth and zenith angles.  These are
   computed for each pixel based on lat/long.  (New possolp.c file)

* LDCMSR-v3.5.2.f
4. The overall algorithm now computes the lat/long for each pixel and then
   determines the solar azimuth and solar zenith angles by calling possolp
   for each pixel (sza, saa).  The viewzenith and viewazimuth are also
   approximated.  These per-pixel angles are used in the subaeroret
   computations.  Note: the view zenith/azimuth corrections will be minor
   in the results due to the small view angle window of L8 (range of
   approximately 8 degrees).

5. The ratiob1, ratiob2, slpratiob1, slpratiob2, intratiob1, and intratiob2
   values read from the ratiomapndwiexp.hdf file have changed which SDSs are
   used.  In general, the b1 values have switched from using b3 to b10 in the
   HDF file.  The b2 values have switched from using b8 to b9 in the HDF file.
   Eric indicated this change in SDSs provides a better match between the MODIS
   data in the HDF file as the Landsat imagery.

6. The four non-fill image corners are found and the scene center is computed
   from those four corners.  The image rotation angle is computed using these
   corners as well.  These corners and rotation angle are used in the
   approximation of the per-pixel view zenith angle and view azimuth angle.
   They are used to produce an index into the l8geom.hdf file which has a
   table of view zenith and view azimuth values.

7. The TOA reflectance corrections use the approximated solar zenith angle for
   the current pixel versus the solar zenith angle for the scene center.

8. Previously cirrus pixels were flagged as clouds if the value in band 9
   was greater than 100.0 / tp[curr_pix] / 1013.  These pixels were not run
   through aerosol inversion.  The new code does not flag these pixels and
   therefore aerosol inversion is run on all pixels.

9. If/then statements were changed for the slpratio and intratio processing.
   The values of intratiob1 and intratiob2 were changed from 327, 482 to
   550, 600. 

10. The actual call to subaeroret for the aerosol retrievals now uses a
    pressure, ozone, and water vapor value that was interpolated for the
    current pixel versus a constant value for each variable.

11. subaeroret is called three times with eps = 1.0, 1.75, and 2.5 for each
    pixel.  The algorithm looks for the eps that minimizes the residual.  Then
    that eps, residual, and raot are carried forward.  Previously subaeroret
    was only called once, and the eps was not used.

12. The if-check on the model residual in relation to the aerosol impact has
    been modified.

13. The model residual is checked with several checks, including the NDVI check.
    Pixels failing are flagged as water.  These pixels are then reprocessed for
    aerosol retrieval, using the new aerosol retrieval algorithm for water.

14. The aerosol interpolation is handled differently.  The average taero and
    teps values in a surrounding 11x11 window are computed and used.  In the
    event that not enough valid values are available, the taeros is set to
    0.05 and teps is set to 1.5 as the defaults.

15. Cloud QA processing is now handled at the end of the algorithm, after the
    aerosol retrieval/interpoloation and the surface reflectance corrections
    have been completed.  Therefore, the algorithm is no longer reliant upon
    cloud QA information, and I confirmed that with Jean-Claude.  I asked him
    if we could simply skip doing the cloud QA processing and he reconfirmed
    that the cloud QA is no longer needed.  The C code will not run the cloud
    QA processing, since there are way better cloud algorithms in existence.
    This also means that the algorithm is no longer reliant upon the TIRS bands.

16. The aerosol QA information has changed slightly and that QA band will now
    also include the aerosol QA confidence that used to be stored in the cloud
    QA band.  This new band will be the "aerosol" band and it replaces the
    "cloud" and "ipflag" bands.

* ESPA Enhancements
17. For collection products, the per-pixel approximations for the solar and
    view angles have been replaced with the angle values from the Landsat
    per-pixel angle bands, which are generated from the Level-1 angle
    coefficient file. This makes the newly delivered l8geom.hdf file obsolete.
    There is a src directory and src_pre_collection directory under the c_code
    for LaSRC.  The src_pre_collection source code still uses the original
    angle approximations, since pre-collection scenes don't have the angle
    coefficient file.  The src directory contains the source code to be used
    for the collection products.  do_lasrc.py has been modified to call the
    different versions of C-code, depending on whether the data is
    pre-collection or collection.
