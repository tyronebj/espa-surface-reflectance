#include "lasrc.h"
#include "time.h"
#include "aero_interp.h"

/******************************************************************************
MODULE:  compute_toa_refl

PURPOSE:  Computes the TOA reflectance and at-sensor brightness temps for all
the bands except the pan band. Uses a per-pixel solar zenith angle for the
TOA corrections. Also determines radiometric saturation for each band, as
available.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
  1. These TOA and BT algorithms match those as published by the USGS Landsat
     team in http://landsat.usgs.gov/Landsat8_Using_Product.php
******************************************************************************/
int compute_toa_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    char *instrument,   /* I: instrument to be processed (OLI, TIRS) */
    int16 *sza,         /* I: scaled per-pixel solar zenith angles (degrees),
                              nlines x nsamps */
    int16 **sband,      /* O: output TOA reflectance and brightness temp
                              values (scaled) */
    uint16 *radsat      /* O: radiometric saturation QA band, nlines x nsamps;
                              array should be all zeros on input to this
                              routine*/
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_toa_refl";   /* function name */
    int i;               /* looping variable for pixels */
    int line, samp;      /* looping variables for lines and samples */
    int ib;              /* looping variable for input bands */
    int sband_ib;        /* looping variable for output bands */
    int iband;           /* current band */
    float rotoa;         /* top of atmosphere reflectance */
    float tmpf;          /* temporary floating point value */
    float refl_mult;     /* reflectance multiplier for bands 1-9 */
    float refl_add;      /* reflectance additive for bands 1-9 */
    float xcals;         /* radiance multiplier for bands 10 and 11 */
    float xcalo;         /* radiance additive for bands 10 and 11 */
    float k1b10;         /* K1 temperature constant for band 10 */
    float k1b11;         /* K1 temperature constant for band 11 */
    float k2b10;         /* K2 temperature constant for band 10 */
    float k2b11;         /* K2 temperature constant for band 11 */
    float xmus;          /* cosine of solar zenith angle (per-pixel) */
    uint16 *uband = NULL;  /* array for input image data for a single band,
                              nlines x nsamps */

time_t mytime;
mytime = time(NULL);
printf ("DEBUG: Start compute_toa_refl: %s\n", ctime(&mytime));
    /* Allocate memory for band data */
    uband = calloc (nlines*nsamps, sizeof (uint16));
    if (uband == NULL)
    {
        sprintf (errmsg, "Error allocating memory for uband");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the bands (except the pan band) and compute the TOA
       reflectance and at-sensor brightness temp */
    for (ib = DN_BAND1; ib <= DN_BAND11; ib++)
    {
        /* Don't process the pan band */
        if (ib == DN_BAND8)
            continue;
        printf ("%d ... ", ib+1);

        /* Read the current band and calibrate bands 1-9 (except pan) to
           obtain TOA reflectance. Bands are corrected for the sun angle at
           the center of the scene. */
        if (ib <= DN_BAND9)
        {
            if (ib <= DN_BAND7)
            {
                iband = ib;
                sband_ib = ib;
            }
            else
            {  /* don't count the pan band */
                iband = ib - 1;
                sband_ib = ib - 1;
            }

            if (get_input_refl_lines (input, iband, 0, nlines, uband) !=
                SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get TOA reflectance coefficients for this reflectance band from
               XML file */
            refl_mult = input->meta.gain[iband];
            refl_add = input->meta.bias[iband];

#ifdef _OPENMP
            #pragma omp parallel for private (line, samp, i, xmus, rotoa)
#endif
            for (line = 0; line < nlines; line++)
            {
                i = line * nsamps;
                for (samp = 0; samp < nsamps; samp++, i++)
                {
                    /* If this pixel is not fill */
                    if (!level1_qa_is_fill (qaband[i]))
                    {
                        /* Compute the TOA reflectance based on the per-pixel
                           sun angle (need to unscale). Scale the TOA value for
                           output. */
                        xmus = cos(sza[i] * 0.01 * DEG2RAD);
                        rotoa = (uband[i] * refl_mult) + refl_add;
                        rotoa = rotoa * MULT_FACTOR / xmus;
    
                        /* Save the scaled TOA reflectance value, but make
                           sure it falls within the defined valid range. */
                        if (rotoa < MIN_VALID)
                            sband[sband_ib][i] = MIN_VALID;
                        else if (rotoa > MAX_VALID)
                            sband[sband_ib][i] = MAX_VALID;
                        else
                            sband[sband_ib][i] = (int) (roundf (rotoa));

                        /* Check for saturation. Saturation is when the pixel
                           reaches the max allowed value. */
                        if (uband[i] == L1_SATURATED)
                            radsat[i] |= 1 << (ib+1);
                    }
                    else
                    {
                        sband[sband_ib][i] = FILL_VALUE;
                        radsat[i] = RADSAT_FILL_VALUE;
                    }
                }  /* for samp */
            }  /* for line */
        }  /* end if band <= band 9 */

        /* Read the current band and calibrate thermal bands.  Not available
           for OLI-only scenes. */
        else if (ib == DN_BAND10 && strcmp (instrument, "OLI"))
        {
            if (get_input_th_lines (input, 0, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[0];
            xcalo = input->meta.bias_th[0];
            k1b10 = input->meta.k1_const[0];
            k2b10 = input->meta.k2_const[0];

            /* Compute brightness temp for band 10.  Make sure it falls
               within the min/max range for the thermal bands. */
#ifdef _OPENMP
            #pragma omp parallel for private (i, tmpf)
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill */
                if (!level1_qa_is_fill (qaband[i]))
                {
                    /* Compute the TOA spectral radiance */
                    tmpf = xcals * uband[i] + xcalo;

                    /* Compute the at-satellite brightness temp (K) and
                       scale for output */
                    tmpf = k2b10 / log (k1b10 / tmpf + 1.0);
                    tmpf = tmpf * MULT_FACTOR_TH;  /* scale the value */

                    /* Make sure the brightness temp falls within the specified
                       range */
                    if (tmpf < MIN_VALID_TH)
                        sband[SR_BAND10][i] = MIN_VALID_TH;
                    else if (tmpf > MAX_VALID_TH)
                        sband[SR_BAND10][i] = MAX_VALID_TH;
                    else
                        sband[SR_BAND10][i] = (int) (roundf (tmpf));

                    /* Check for saturation */
                    if (uband[i] == L1_SATURATED)
                        radsat[i] |= 1 << (ib+1);
                }
                else
                {
                    sband[SR_BAND10][i] = FILL_VALUE;
                    radsat[i] = RADSAT_FILL_VALUE;
                }
            }
        }  /* end if band 10 */

        else if (ib == DN_BAND11 && strcmp (instrument, "OLI"))
        {
            if (get_input_th_lines (input, 1, 0, nlines, uband) != SUCCESS)
            {
                sprintf (errmsg, "Reading band %d", ib+1);
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            /* Get brightness temp coefficients for this band from XML file */
            xcals = input->meta.gain_th[1];
            xcalo = input->meta.bias_th[1];
            k1b11 = input->meta.k1_const[1];
            k2b11 = input->meta.k2_const[1];

            /* Compute brightness temp for band 11.  Make sure it falls
               within the min/max range for the thermal bands. */
#ifdef _OPENMP
            #pragma omp parallel for private (i, tmpf)
#endif
            for (i = 0; i < nlines*nsamps; i++)
            {
                /* If this pixel is not fill */
                if (!level1_qa_is_fill (qaband[i]))
                {
                    /* Compute the TOA spectral radiance */
                    tmpf = xcals * uband[i] + xcalo;

                    /* Compute the at-satellite brightness temp (K) and
                       scale for output */
                    tmpf = k2b11 / log (k1b11 / tmpf + 1.0);
                    tmpf = tmpf * MULT_FACTOR_TH;  /* scale the value */

                    /* Make sure the brightness temp falls within the specified
                       range */
                    if (tmpf < MIN_VALID_TH)
                        sband[SR_BAND11][i] = MIN_VALID_TH;
                    else if (tmpf > MAX_VALID_TH)
                        sband[SR_BAND11][i] = MAX_VALID_TH;
                    else
                        sband[SR_BAND11][i] = (int) (roundf (tmpf));

                    /* Check for saturation only */
                    if (uband[i] == L1_SATURATED)
                        radsat[i] |= 1 << (ib+1);
                }
                else
                {
                    sband[SR_BAND11][i] = FILL_VALUE;
                    radsat[i] = RADSAT_FILL_VALUE;
                }
            }
        }  /* end if band 11 */
    }  /* end for ib */
    printf ("\n");

    /* The input data has been read and calibrated. The memory can be freed. */
    free (uband);

mytime = time(NULL);
printf ("DEBUG: End compute_toa_refl: %s\n", ctime(&mytime));

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  compute_sr_refl

PURPOSE:  Computes the surfance reflectance for all the reflectance bands.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Initializes the variables and data arrays from the lookup table and
   auxiliary files.
2. The tauray array was originally read in from a static ASCII file, but it is
   now hardcoded to save time from reading the file each time.  This file was
   generated (like many of the other auxiliary input tables) by running 6S and
   storing the coefficients.
4. Aerosols are retrieved for all non-fill pixels.  If the aerosol fails the
   model residual or NDVI test, then the pixel is flagged as water.  All water
   pixels are run through a water-specific aerosol retrieval.  If the model
   residual fails, then that pixel is marked as failed aerosol retrieval.  Any
   pixel that failed retrieval is then interpolated using an average of the
   clear (valid land pixel aerosols) and water (valid water pixel aerosols).
   Those final aerosol values are used for the surface reflectance corrections.
5. Cloud-based QA information is not processed in this algorithm.
******************************************************************************/
int compute_sr_refl
(
    Input_t *input,     /* I: input structure for the Landsat product */
    Espa_internal_meta_t *xml_metadata,
                        /* I: XML metadata structure */
    char *xml_infile,   /* I: input XML filename */
    uint16 *qaband,     /* I: QA band for the input image, nlines x nsamps */
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    float pixsize,      /* I: pixel size for the reflectance bands */
    int16 **sband,      /* I/O: input TOA and output surface reflectance */
    int16 *sza,         /* I: scaled per-pixel solar zenith angles (degrees),
                              nlines x nsamps */
    int16 *saa,         /* I: scaled per-pixel solar azimuth angles (degrees),
                              nlines x nsamps */
    int16 *vza,         /* I: scaled per-pixel view zenith angles (degrees),
                              nlines x nsamps */
    int16 *vaa,         /* I: scaled per-pixel view azimuth angles (degrees),
                              nlines x nsamps */
    float xts,          /* I: solar zenith angle (deg) */
    float xfs,          /* I: solar azimuth angle (deg) */
    float xmus,         /* I: cosine of solar zenith angle */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm         /* I: auxiliary filename for ozone and water vapor */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "compute_sr_refl";   /* function name */
    int retval;          /* return status */
    int i, j, k, l;      /* looping variable for pixels */
    int ib;              /* looping variable for input bands */
    int iband;           /* current band */
    int curr_pix;        /* current pixel in 1D arrays of nlines * nsamps */
    int center_pix;      /* current pixel in 1D arrays of nlines * nsamps for
                            the center of the aerosol window */
    int win_pix;         /* current pixel in the line,sample window */
    int center_line;     /* line for the center of the aerosol window */
    int center_samp;     /* sample for the center of the aerosol window */
    int nearest_line;    /* line for nearest non-fill/cloud pixel in the
                            aerosol window */
    int nearest_samp;    /* samp for nearest non-fill/cloud pixel in the
                            aerosol window */
    float tmpf;          /* temporary floating point value */
    float rotoa;         /* top of atmosphere reflectance */
    float roslamb;       /* lambertian surface reflectance */
    float tgo;           /* other gaseous transmittance (tgog * tgoz) */
    float roatm;         /* atmospheric intrinsic reflectance */
    float ttatmg;        /* total atmospheric transmission */
    float satm;          /* atmosphere spherical albedo */
    float xrorayp;       /* reflectance of the atmosphere due to molecular
                            (Rayleigh) scattering */
    float next;
    float erelc[NSR_BANDS];    /* band ratio variable for bands 1-7 */
    float troatm[NSR_BANDS];   /* atmospheric reflectance table for bands 1-7 */
    float btgo[NSR_BANDS];     /* other gaseous transmittance for bands 1-7 */
    float broatm[NSR_BANDS];   /* atmospheric reflectance for bands 1-7 */
    float bttatmg[NSR_BANDS];  /* ttatmg for bands 1-7 */
    float bsatm[NSR_BANDS];    /* atmosphere spherical albedo for bands 1-7 */

    int iband1, iband3; /* band indices (zero-based) */
    float raot;         /* AOT reflectance */
    float sraot1, sraot2, sraot3;
                        /* raot values for three different eps values */
    float residual;     /* model residual */
    float residual1, residual2, residual3;
                        /* residuals for 3 different eps values */
    float rsurf;        /* surface reflectance */
    float corf;         /* aerosol impact (higher values represent high
                           aerosol) */

    float ros1,ros4,ros5; /* surface reflectance for bands 1, 4, and 5 */
    int tmp_percent;      /* current percentage for printing status */
#ifndef _OPENMP
    int curr_tmp_percent; /* percentage for current line */
#endif

    float lat, lon;       /* pixel lat, long location */
    int lcmg, scmg;       /* line/sample index for the CMG */
    int lcmg1, scmg1;     /* line+1/sample+1 index for the CMG */
    float u, v;           /* line/sample index for the CMG */
    float one_minus_u;    /* 1.0 - u */
    float one_minus_v;    /* 1.0 - v */
    float one_minus_u_x_one_minus_v;  /* (1.0 - u) * (1.0 - v) */
    float one_minus_u_x_v;  /* (1.0 - u) * v */
    float u_x_one_minus_v;  /* u * (1.0 - v) */
    float u_x_v;          /* u * v */
    float ndwi_th1, ndwi_th2; /* values for NDWI calculations */
    float xcmg, ycmg;     /* x/y location for CMG */
    float xndwi;          /* calculated NDWI value */
    int uoz11, uoz21, uoz12, uoz22;  /* ozone at line,samp; line, samp+1;
                                        line+1, samp; and line+1, samp+1 */
    float pres11, pres12, pres21, pres22;  /* pressure at line,samp;
                             line, samp+1; line+1, samp; and line+1, samp+1 */
    float wv11, wv12, wv21, wv22;  /* water vapor at line,samp;
                             line, samp+1; line+1, samp; and line+1, samp+1 */
    uint8 *ipflag = NULL; /* QA flag to assist with aerosol interpolation,
                             nlines x nsamps */
    float *twvi = NULL;   /* interpolated water vapor value,
                             nlines x nsamps */
    float *tozi = NULL;   /* interpolated ozone value, nlines x nsamps */
    float *tp = NULL;     /* interpolated pressure value, nlines x nsamps */
    float *taero = NULL;  /* aerosol values for each pixel, nlines x nsamps */
    float *teps = NULL;   /* angstrom coeff for each pixel, nlines x nsamps */
    int16 *aerob1 = NULL; /* atmospherically corrected band 1 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob2 = NULL; /* atmospherically corrected band 2 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob4 = NULL; /* atmospherically corrected band 4 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob5 = NULL; /* atmospherically corrected band 5 data
                             (TOA refl), nlines x nsamps */
    int16 *aerob7 = NULL; /* atmospherically corrected band 7 data
                             (TOA refl), nlines x nsamps */

    /* Vars for forward/inverse mapping space */
    Geoloc_t *space = NULL;       /* structure for geolocation information */
    Space_def_t space_def;        /* structure to define the space mapping */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */

    /* Lookup table variables */
    float eps;           /* angstrom coefficient */
    float eps1, eps2, eps3;  /* eps values for three runs */
    float xtv;           /* observation zenith angle (deg) */
    float xmuv;          /* cosine of observation zenith angle */
    float xfi;           /* azimuthal difference between the sun and
                            observation angle (deg) */
    float cosxfi;        /* cosine of azimuthal difference */
    float xtsstep;       /* solar zenith step value */
    float xtsmin;        /* minimum solar zenith value */
    float xtvstep;       /* observation step value */
    float xtvmin;        /* minimum observation value */
    float *rolutt = NULL;  /* intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt = NULL;  /* transmission table
                       [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUNANGLE_VALS] */
    float *sphalbt = NULL; /* spherical albedo table 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext = NULL; /* aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *tsmax = NULL;   /* maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin = NULL;   /* minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi = NULL;    /* number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfic = NULL;   /* communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *ttv = NULL;     /* view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22];         /* sun angle table */
    int32 indts[22];       /* index for sun angle table */
    int iaots;             /* index for AOTs */

    /* Auxiliary file variables */
    int16 *dem = NULL;        /* CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi = NULL;      /* avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi = NULL;      /* standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1 = NULL;    /* mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2 = NULL;    /* mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7 = NULL;    /* mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1 = NULL;   /* intercept band1 ratio,
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob2 = NULL;   /* intercept band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *intratiob7 = NULL;   /* intercept band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob1 = NULL;   /* slope band1 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob2 = NULL;   /* slope band2 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    int16 *slpratiob7 = NULL;   /* slope band7 ratio
                                   RATIO_NBLAT x RATIO_NBLON */
    uint16 *wv = NULL;       /* water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz = NULL;        /* ozone values [CMG_NBLAT x CMG_NBLON] */
    float raot550nm;    /* nearest input value of AOT */
    float uoz;          /* total column ozone */
    float uwv;          /* total column water vapor (precipital water vapor) */
    float pres;         /* surface pressure */
    float rb1;          /* band ratio 1 (unscaled) */
    float rb2;          /* band ratio 2 (unscaled) */
    float slpr11, slpr12, slpr21, slpr22;  /* band ratio slope at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float intr11, intr12, intr21, intr22;  /* band ratio intercept at line,samp;
                           line, samp+1; line+1, samp; and line+1, samp+1 */
    float slprb1, slprb2, slprb7;  /* interpolated band ratio slope values for
                                      band ratios 1, 2, 7 */
    float intrb1, intrb2, intrb7;  /* interpolated band ratio intercept values
                                      for band ratios 1, 2, 7 */
    int ratio_pix11;  /* pixel location for ratio products [lcmg][scmg] */
    int ratio_pix12;  /* pixel location for ratio products [lcmg][scmg+1] */
    int ratio_pix21;  /* pixel location for ratio products [lcmg+1][scmg] */
    int ratio_pix22;  /* pixel location for ratio products [lcmg+1][scmg+1] */
    int cmg_pix11;    /* pixel location for CMG/DEM products [lcmg][scmg] */
    int cmg_pix12;    /* pixel location for CMG/DEM products [lcmg][scmg+1] */
    int cmg_pix21;    /* pixel location for CMG/DEM products [lcmg+1][scmg] */
    int cmg_pix22;    /* pixel location for CMG/DEM products [lcmg+1][scmg+1] */

    /* Variables for the aerosol interpolation */
// TODO GAIL smflag, taeros, and tepss are not used unless the aerosol
// interpolation is being done after the fact.  And, if the NxN window is in
// place this won't be needed.
    bool *smflag = NULL;  /* flag for whether or not the window average was
                             computed and is valid for this pixel */
    int nbpixnf;          /* number of non-filled aerosol pixels */
    int prev_nbpixnf;     /* number of non-filled aerosol pixels in previous
                             loop */
    int nbpixtot;         /* total number of pixels in the window */
    float taeroavg;       /* average of the taero values in the window */
    float tepsavg;        /* average of the teps values in the window */
    int nbaeroavg;        /* number of pixels in the window used for computing
                             the taeroavg and tepsavg */
    float *taeros = NULL; /* array of the average taero values in the local
                             window, nlines x nsamps */
    float *tepss = NULL;  /* array of the average teps values in the local
                             window, nlines x nsamps */

    /* Variables for finding the eps that minimizes the residual */
    double xa, xb, xc, xd, xe, xf;  /* coefficients */
    double coefa, coefb;            /* coefficients */
    float epsmin;                  /* eps which minimizes the residual */

    /* Output file info */
    Output_t *sr_output = NULL;  /* output structure and metadata for the SR
                                    product */
    Envi_header_t envi_hdr;      /* output ENVI header information */
    char envi_file[STR_SIZE];    /* ENVI filename */
    char *cptr = NULL;           /* pointer to the file extension */

    /* Table constants */
    float aot550nm[NAOT_VALS] =  /* AOT look-up table */
        {0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.20,
         1.40, 1.60, 1.80, 2.00, 2.30, 2.60, 3.00, 3.50, 4.00, 4.50, 5.00};
    float tpres[NPRES_VALS] =      /* surface pressure table */
        {1050.0, 1013.0, 900.0, 800.0, 700.0, 600.0, 500.0};

    /* Atmospheric correction variables */
    /* Look up table for atmospheric and geometric quantities */
    float tauray[NSR_BANDS] =  /* molecular optical thickness coefficients --
        produced by running 6S */
        {0.23638, 0.16933, 0.09070, 0.04827, 0.01563, 0.00129, 0.00037,
         0.07984};
    double oztransa[NSR_BANDS] =   /* ozone transmission coeff */
        {-0.00255649, -0.0177861, -0.0969872, -0.0611428, 0.0001, 0.0001,
          0.0001, -0.0834061};
    double wvtransa[NSR_BANDS] =   /* water vapor transmission coeff */
        {2.29849e-27, 2.29849e-27, 0.00194772, 0.00404159, 0.000729136,
         0.00067324, 0.0177533, 0.00279738};
    double wvtransb[NSR_BANDS] =   /* water vapor transmission coeff */
        {0.999742, 0.999742, 0.775024, 0.774482, 0.893085, 0.939669, 0.65094,
         0.759952};
    double ogtransa1[NSR_BANDS] =  /* other gases transmission coeff */
        {4.91586e-20, 4.91586e-20, 4.91586e-20, 1.04801e-05, 1.35216e-05,
         0.0205425, 0.0256526, 0.000214329};
    double ogtransb0[NSR_BANDS] =  /* other gases transmission coeff */
        {0.000197019, 0.000197019, 0.000197019, 0.640215, -0.195998, 0.326577,
         0.243961, 0.396322};
    double ogtransb1[NSR_BANDS] =  /* other gases transmission coeff */
        {9.57011e-16, 9.57011e-16, 9.57011e-16, -0.348785, 0.275239, 0.0117192,
         0.0616101, 0.04728};

#define WRITE_TAERO 1
#ifdef WRITE_TAERO
    FILE *aero_fptr=NULL;
#endif

time_t mytime;
mytime = time(NULL);
printf ("DEBUG: Start compute_sr_refl: %s\n", ctime(&mytime));

    /* Allocate memory for the many arrays needed to do the surface reflectance
       computations */
    retval = memory_allocation_sr (nlines, nsamps, &aerob1, &aerob2, &aerob4,
        &aerob5, &aerob7, &ipflag, &twvi, &tozi, &tp, &taero, &taeros, &teps,
        &tepss, &smflag, &dem, &andwi, &sndwi, &ratiob1, &ratiob2, &ratiob7,
        &intratiob1, &intratiob2, &intratiob7, &slpratiob1, &slpratiob2,
        &slpratiob7, &wv, &oz, &rolutt, &transt, &sphalbt, &normext, &tsmax,
        &tsmin, &nbfic, &nbfi, &ttv);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error allocating memory for the data arrays needed "
            "for surface reflectance calculations.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the geolocation space applications */
    if (!get_geoloc_info (xml_metadata, &space_def))
    {
        sprintf (errmsg, "Getting the space definition from the XML file");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    space = setup_mapping (&space_def);
    if (space == NULL)
    {
        sprintf (errmsg, "Setting up the geolocation mapping");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Initialize the look up tables and atmospheric correction variables */
    retval = init_sr_refl (nlines, nsamps, input, space, anglehdf, intrefnm,
        transmnm, spheranm, cmgdemnm, rationm, auxnm, &eps, &iaots, &xtv,
        &xmuv, &xfi, &cosxfi, &raot550nm, &pres, &uoz, &uwv, &xtsstep,
        &xtsmin, &xtvstep, &xtvmin, tsmax, tsmin, tts, ttv, indts, rolutt,
        transt, sphalbt, normext, nbfic, nbfi, dem, andwi, sndwi, ratiob1,
        ratiob2, ratiob7, intratiob1, intratiob2, intratiob7, slpratiob1,
        slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Error initializing the lookup tables and "
            "atmospheric correction variables.");
        error_handler (false, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through all the reflectance bands and perform atmospheric
       corrections based on climatology */
mytime = time(NULL);
printf ("DEBUG: Start atm corrections: %s\n", ctime(&mytime));
    printf ("Performing atmospheric corrections for each reflectance "
        "band ...");
    for (ib = 0; ib <= SR_BAND7; ib++)
    {
        printf (" %d ...", ib+1);

        /* Get the parameters for the atmospheric correction */
        /* rotoa is not defined for this call, which is ok, but the
           roslamb value is not valid upon output. Just set it to 0.0 to
           be consistent. */
        rotoa = 0.0;
        raot550nm = aot550nm[1];
        eps = 2.5;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
            raot550nm, ib, pres, tpres, aot550nm, rolutt, transt, xtsstep,
            xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic,
            nbfi, tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
            ogtransb1, wvtransa, wvtransb, oztransa, rotoa, &roslamb,
            &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Save these band-related parameters for later */
        btgo[ib] = tgo;
        broatm[ib] = roatm;
        bttatmg[ib] = ttatmg;
        bsatm[ib] = satm;

        /* Perform atmospheric corrections for bands 1-7 */
#ifdef _OPENMP
        #pragma omp parallel for private (i, j, curr_pix, rotoa, roslamb)
#endif
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If this pixel is not fill.  Otherwise fill pixels have
                   already been marked in the TOA calculations. */
                if (!level1_qa_is_fill (qaband[curr_pix]))
                {
                    /* Store the TOA scaled TOA reflectance values for later
                       use before completing atmospheric corrections */
                    if (ib == DN_BAND1)
                        aerob1[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND2)
                        aerob2[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND4)
                        aerob4[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND5)
                        aerob5[curr_pix] = sband[ib][curr_pix];
                    else if (ib == DN_BAND7)
                        aerob7[curr_pix] = sband[ib][curr_pix];
    
                    /* Apply the atmospheric corrections (ignoring the Rayleigh
                       scattering component and water vapor), and store the
                       scaled value for further corrections.  (NOTE: the full
                       computations are in atmcorlamb2) */
                    rotoa = sband[ib][curr_pix] * SCALE_FACTOR;
                    roslamb = rotoa / tgo;
                    roslamb = roslamb - roatm;
                    roslamb = roslamb / ttatmg;
                    roslamb = roslamb / (1.0 + satm * roslamb);
                    sband[ib][curr_pix] = (int) (roslamb * MULT_FACTOR);
                }
            }  /* end for j */
        }  /* end for i */
    }  /* for ib */
    printf ("\n");

    /* Interpolate the auxiliary data for each pixel location */
mytime = time(NULL);
printf ("DEBUG: Start interpolating the auxiliary data: %s\n", ctime(&mytime));
    printf ("Interpolating the auxiliary data ...\n");
    tmp_percent = 0;
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, curr_pix, img, geo, lat, lon, xcmg, ycmg, lcmg, scmg, lcmg1, scmg1, u, v, one_minus_u, one_minus_v, one_minus_u_x_one_minus_v, one_minus_u_x_v, u_x_one_minus_v, u_x_v, cmg_pix11, cmg_pix12, cmg_pix21, cmg_pix22, wv11, wv12, wv21, wv22, uoz11, uoz12, uoz21, uoz22, pres11, pres12, pres21, pres22)
#endif

    for (i = 0; i < nlines; i++)
    {
#ifndef _OPENMP
        /* update status, but not if multi-threaded */
        curr_tmp_percent = 100 * i / nlines;
        if (curr_tmp_percent > tmp_percent)
        {
            tmp_percent = curr_tmp_percent;
            if (tmp_percent % 10 == 0)
            {
                printf ("%d%% ", tmp_percent);
                fflush (stdout);
            }
        }
#endif

        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If this pixel is fill, then don't process */
            if (qaband[curr_pix] == 1)
            {
                ipflag[curr_pix] |= (1 << IPFLAG_FILL);
                continue;
            }

            /* Get the lat/long for the current pixel */
            img.l = i - 0.5;
            img.s = j + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                    "geolocation coords", i, j);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /*** Handle all the variables related to the current pixel in the
                 auxiliary products ***/
            /* Use that lat/long to determine the line/sample in the
               CMG-related lookup tables, using the center of the UL
               pixel. Note, we are basically making sure the line/sample
               combination falls within -90, 90 and -180, 180 global climate
               data boundaries.  However, the source code below uses lcmg+1
               and scmg+1, which for some scenes may wrap around the dateline
               or the poles.  Thus we need to wrap the CMG data around to the
               beginning of the array. */
            /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
               pixel for each calculation.  Negative latitude values should be
               the largest line values in the CMG grid.  Negative longitude
               values should be the smallest sample values in the CMG grid. */
            /* The line/sample calculation from the x/ycmg values are not
               rounded.  The interpolation of the value using line+1 and
               sample+1 are based on the truncated numbers, therefore rounding
               up is not appropriate. */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) ycmg;
            scmg = (int) xcmg;

            /* Handle the edges of the lat/long values in the CMG grid */
            if (lcmg < 0)
                lcmg = 0;
            else if (lcmg >= CMG_NBLAT)
                lcmg = CMG_NBLAT;

            if (scmg < 0)
                scmg = 0;
            else if (scmg >= CMG_NBLON)
                scmg = CMG_NBLON;

            /* If the current CMG pixel is at the edge of the CMG array,
               then allow the next pixel for interpolation to wrap around
               the array */
            if (scmg >= CMG_NBLON-1)  /* 180 degrees so wrap around */
                scmg1 = 0;
            else
                scmg1 = scmg + 1;

            if (lcmg >= CMG_NBLAT-1)  /* -90 degrees so wrap around */
                lcmg1 = 0;
            else
                lcmg1 = lcmg + 1;

            /* Determine the four CMG pixels to be used for the current
               Landsat pixel */
            cmg_pix11 = lcmg * CMG_NBLON + scmg;
            cmg_pix12 = lcmg * CMG_NBLON + scmg1;
            cmg_pix21 = lcmg1 * CMG_NBLON + scmg;
            cmg_pix22 = lcmg1 * CMG_NBLON + scmg1;

            /* Get the water vapor pixels. If the water vapor value is
               fill (=0), then use it as-is. */
            wv11 = wv[cmg_pix11];
            wv12 = wv[cmg_pix12];
            wv21 = wv[cmg_pix21];
            wv22 = wv[cmg_pix22];

            /* Get the ozone pixels. If the ozone value is fill (=0), then use
               a default value of 120. */
            uoz11 = oz[cmg_pix11];
            if (uoz11 == 0)
                uoz11 = 120;

            uoz12 = oz[cmg_pix12];
            if (uoz12 == 0)
                uoz12 = 120;

            uoz21 = oz[cmg_pix21];
            if (uoz21 == 0)
                uoz21 = 120;

            uoz22 = oz[cmg_pix22];
            if (uoz22 == 0)
                uoz22 = 120;

            /* Get the surface pressure from the global DEM.  Set to 1013.0
               (sea level) if the DEM is fill (= -9999), which is likely ocean.
               The dimensions on the DEM array is the same as that of the CMG
               arrays. Use the current pixel locations already calculated. */
            if (dem[cmg_pix11] != -9999)
                pres11 = 1013.0 * exp (-dem[cmg_pix11] * ONE_DIV_8500);
            else
                pres11 = 1013.0;

            if (dem[cmg_pix12] != -9999)
                pres12 = 1013.0 * exp (-dem[cmg_pix12] * ONE_DIV_8500);
            else
                pres12 = 1013.0;

            if (dem[cmg_pix21] != -9999)
                pres21 = 1013.0 * exp (-dem[cmg_pix21] * ONE_DIV_8500);
            else
                pres21 = 1013.0;

            if (dem[cmg_pix22] != -9999)
                pres22 = 1013.0 * exp (-dem[cmg_pix22] * ONE_DIV_8500);
            else
                pres22 = 1013.0;

            /*** Handle all the variables related to the current pixel in the
                 Landsat scene, which means interpolating the global-level
                 variables ***/
            /* Determine the fractional difference between the integer location
               and floating point pixel location to be used for interpolation */
            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            one_minus_u = 1.0 - u;
            one_minus_v = 1.0 - v;
            one_minus_u_x_one_minus_v = one_minus_u * one_minus_v;
            one_minus_u_x_v = one_minus_u * v;
            u_x_one_minus_v = u * one_minus_v;
            u_x_v = u * v;

            /* Interpolate water vapor, and unscale */
            twvi[curr_pix] = wv11 * one_minus_u_x_one_minus_v +
                             wv12 * one_minus_u_x_v +
                             wv21 * u_x_one_minus_v +
                             wv22 * u_x_v;
            twvi[curr_pix] = twvi[curr_pix] * 0.01;   /* vs / 100 */

            /* Interpolate ozone, and unscale */
            tozi[curr_pix] = uoz11 * one_minus_u_x_one_minus_v +
                             uoz12 * one_minus_u_x_v +
                             uoz21 * u_x_one_minus_v +
                             uoz22 * u_x_v;
            tozi[curr_pix] = tozi[curr_pix] * 0.0025;   /* vs / 400 */


            /* Interpolate surface pressure */
            tp[curr_pix] = pres11 * one_minus_u_x_one_minus_v +
                           pres12 * one_minus_u_x_v +
                           pres21 * u_x_one_minus_v +
                           pres22 * u_x_v;
        }  /* end for j */
    }  /* end for i */

#ifndef _OPENMP
    /* update status */
    printf ("100%%\n");
    fflush (stdout);
#endif

    /* Start the aerosol inversion */
mytime = time(NULL);
printf ("DEBUG: Start Aerosol Inversion: %s\n", ctime(&mytime));
    printf ("Aerosol Inversion using %d x %d aerosol window ...\n", AERO_WINDOW,
        AERO_WINDOW);
    tmp_percent = 0;
#ifdef _OPENMP
    #pragma omp parallel for private (i, j, center_line, center_samp, nearest_line, nearest_samp, curr_pix, center_pix, img, geo, lat, lon, xcmg, ycmg, lcmg, scmg, lcmg1, scmg1, u, v, one_minus_u, one_minus_v, one_minus_u_x_one_minus_v, one_minus_u_x_v, u_x_one_minus_v, u_x_v, ratio_pix11, ratio_pix12, ratio_pix21, ratio_pix22, rb1, rb2, slpr11, slpr12, slpr21, slpr22, intr11, intr12, intr21, intr22, slprb1, slprb2, slprb7, intrb1, intrb2, intrb7, xndwi, ndwi_th1, ndwi_th2, xtv, xts, xmus, xmuv, xfi, cosxfi, iband, iband1, iband3, pres, uoz, uwv, iaots, retval, eps, eps1, eps2, eps3, residual, residual1, residual2, residual3, raot, sraot1, sraot2, sraot3, xa, xb, xc, xd, xe, xf, coefa, coefb, epsmin, corf, next, rotoa, raot550nm, roslamb, tgo, roatm, ttatmg, satm, xrorayp, ros5, ros4, ros1, erelc, troatm)
#endif
    for (i = HALF_AERO_WINDOW; i < nlines; i += AERO_WINDOW)
    {
#ifndef _OPENMP
        /* update status, but not if multi-threaded */
        curr_tmp_percent = 100 * i / nlines;
        if (curr_tmp_percent > tmp_percent)
        {
            tmp_percent = curr_tmp_percent;
            if (tmp_percent % 10 == 0)
            {
                printf ("%d%% ", tmp_percent);
                fflush (stdout);
            }
        }
#endif

        curr_pix = i * nsamps + HALF_AERO_WINDOW;
        for (j = HALF_AERO_WINDOW; j < nsamps;
             j += AERO_WINDOW, curr_pix += AERO_WINDOW)
        {
            /* Keep track of the center pixel for the current aerosol window;
               may need to return here if this is fill, cloudy or water */
            center_line = i;
            center_samp = j;
            center_pix = curr_pix;
//printf ("DEBUG: Processing aerosol window line, sample: %d, %d\n", center_line, center_samp);

            /* If this pixel is fill */
            if (level1_qa_is_fill (qaband[curr_pix]))
            {
//printf ("DEBUG:    Found fill pixel at line, sample: %d, %d\n", i, j);
                /* Look for other non-fill pixels in the window */
                if (find_closest_non_fill (qaband, nlines, nsamps, center_line,
                    center_samp, &nearest_line, &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill pixel for
                       further processing of aerosols. However we will still
                       write to the center of the aerosol window for the
                       current window. */
                    i = nearest_line;
                    j = nearest_samp;
                    curr_pix = i * nsamps + j;
//printf ("DEBUG:    Closest non-fill pixel at line, sample: %d, %d\n", i, j);
                }
                else
                {
                    /* No other non-fill pixels found.  Pixel is already
                       flagged as fill. Move to next aerosol window. */
                    continue;
                }
            }

            /* If this non-fill pixel is cloud, then look for a pixel which
               is not cloudy or fill.  If none are found, then just use this
               pixel. */
            if (is_cloud (qaband[curr_pix]))
            {
//printf ("DEBUG:    Found cloud pixel at line, sample: %d, %d\n", i, j);
                /* Look for other non-fill/non-cloud pixels in the window.
                   Start with the center of the window and search outward. */
                if (find_closest_non_cloud (qaband, nlines, nsamps,
                    center_line, center_samp, &nearest_line, &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill/non-cloud
                       pixel for further processing */
                    i = nearest_line;
                    j = nearest_samp;
                    curr_pix = i * nsamps + j;
//printf ("DEBUG:    Closest non-cloud pixel at line, sample: %d, %d\n", i, j);
                }
            }

            /* If this non-fill/non-cloud pixel is water, then look for a pixel
               which is not cloudy, not fill, and not water.  If none are
               found, then just use this pixel. */
            if (is_water (sband[SR_BAND4][curr_pix],
                          sband[SR_BAND5][curr_pix]))
            {
//printf ("DEBUG:    Found water pixel at line, sample: %d, %d\n", i, j);
                /* Look for other non-fill/non-cloud/non-water pixels in the
                   window.  Start with the center of the window and search
                   outward. */
                if (find_closest_non_water (qaband, sband, nlines, nsamps,
                    center_line, center_samp, &nearest_line, &nearest_samp))
                {
                    /* Use the line/sample location of the non-fill/non-cloud/
                       non-water pixel for further processing */
                    i = nearest_line;
                    j = nearest_samp;
                    curr_pix = i * nsamps + j;
//printf ("DEBUG:    Closest non-water pixel at line, sample: %d, %d\n", i, j);
                }
            }

            /* Get the lat/long for the current pixel (which may not be the
               center of the aerosol window), for the center of that pixel */
            img.l = i - 0.5;
            img.s = j + 0.5;
            img.is_fill = false;
            if (!from_space (space, &img, &geo))
            {
                sprintf (errmsg, "Mapping line/sample (%d, %d) to "
                    "geolocation coords", i, j);
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }
            lat = geo.lat * RAD2DEG;
            lon = geo.lon * RAD2DEG;

            /* Use that lat/long to determine the line/sample in the
               CMG-related lookup tables, using the center of the UL
               pixel. Note, we are basically making sure the line/sample
               combination falls within -90, 90 and -180, 180 global climate
               data boundaries.  However, the source code below uses lcmg+1
               and scmg+1, which for some scenes may wrap around the
               dateline or the poles.  Thus we need to wrap the CMG data
               around to the beginning of the array. */
            /* Each CMG pixel is 0.05 x 0.05 degrees.  Use the center of the
               pixel for each calculation.  Negative latitude values should
               be the largest line values in the CMG grid.  Negative
               longitude values should be the smallest sample values in the
               CMG grid. */
            /* The line/sample calculation from the x/ycmg values are not
               rounded.  The interpolation of the value using line+1 and
               sample+1 are based on the truncated numbers, therefore
               rounding up is not appropriate. */
            ycmg = (89.975 - lat) * 20.0;   /* vs / 0.05 */
            xcmg = (179.975 + lon) * 20.0;  /* vs / 0.05 */
            lcmg = (int) ycmg;
            scmg = (int) xcmg;

            /* Handle the edges of the lat/long values in the CMG grid */
            if (lcmg < 0)
                lcmg = 0;
            else if (lcmg >= CMG_NBLAT)
                lcmg = CMG_NBLAT;

            if (scmg < 0)
                scmg = 0;
            else if (scmg >= CMG_NBLON)
                scmg = CMG_NBLON;

            /* If the current CMG pixel is at the edge of the CMG array, then
               allow the next pixel for interpolation to wrap around the
               array */
            if (scmg >= CMG_NBLON-1)  /* 180 degrees so wrap around */
                scmg1 = 0;
            else
                scmg1 = scmg + 1;

            if (lcmg >= CMG_NBLAT-1)  /* -90 degrees so wrap around */
                lcmg1 = 0;
            else
                lcmg1 = lcmg + 1;

            /* Determine the fractional difference between the integer location
               and floating point pixel location to be used for interpolation */
            u = (ycmg - lcmg);
            v = (xcmg - scmg);
            one_minus_u = 1.0 - u;
            one_minus_v = 1.0 - v;
            one_minus_u_x_one_minus_v = one_minus_u * one_minus_v;
            one_minus_u_x_v = one_minus_u * v;
            u_x_one_minus_v = u * one_minus_v;
            u_x_v = u * v;

            /* Determine the band ratios and slope/intercept */
            ratio_pix11 = lcmg * RATIO_NBLON + scmg;
            ratio_pix12 = ratio_pix11 + 1;
            ratio_pix21 = lcmg1 * RATIO_NBLON + scmg;
            ratio_pix22 = ratio_pix21 + 1;

            rb1 = ratiob1[ratio_pix11] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix11] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix11] = 0;
                slpratiob2[ratio_pix11] = 0;
                slpratiob7[ratio_pix11] = 0;
                intratiob1[ratio_pix11] = 550;
                intratiob2[ratio_pix11] = 600;
                intratiob7[ratio_pix11] = 2000;
            }
            else if (sndwi[ratio_pix11] < 200)
            {
                slpratiob1[ratio_pix11] = 0;
                slpratiob2[ratio_pix11] = 0;
                slpratiob7[ratio_pix11] = 0;
                intratiob1[ratio_pix11] = ratiob1[ratio_pix11];
                intratiob2[ratio_pix11] = ratiob2[ratio_pix11];
                intratiob7[ratio_pix11] = ratiob7[ratio_pix11];
            }

            rb1 = ratiob1[ratio_pix12] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix12] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix12] = 0;
                slpratiob2[ratio_pix12] = 0;
                slpratiob7[ratio_pix12] = 0;
                intratiob1[ratio_pix12] = 550;
                intratiob2[ratio_pix12] = 600;
                intratiob7[ratio_pix12] = 2000;
            }
            else if (sndwi[ratio_pix12] < 200)
            {
                slpratiob1[ratio_pix12] = 0;
                slpratiob2[ratio_pix12] = 0;
                slpratiob7[ratio_pix12] = 0;
                intratiob1[ratio_pix12] = ratiob1[ratio_pix12];
                intratiob2[ratio_pix12] = ratiob2[ratio_pix12];
                intratiob7[ratio_pix12] = ratiob7[ratio_pix12];
            }

            rb1 = ratiob1[ratio_pix21] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix21] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix21] = 0;
                slpratiob2[ratio_pix21] = 0;
                slpratiob7[ratio_pix21] = 0;
                intratiob1[ratio_pix21] = 550;
                intratiob2[ratio_pix21] = 600;
                intratiob7[ratio_pix21] = 2000;
            }
            else if (sndwi[ratio_pix21] < 200)
            {
                slpratiob1[ratio_pix21] = 0;
                slpratiob2[ratio_pix21] = 0;
                slpratiob7[ratio_pix21] = 0;
                intratiob1[ratio_pix21] = ratiob1[ratio_pix21];
                intratiob2[ratio_pix21] = ratiob2[ratio_pix21];
                intratiob7[ratio_pix21] = ratiob7[ratio_pix21];
            }

            rb1 = ratiob1[ratio_pix22] * 0.001;  /* vs. / 1000. */
            rb2 = ratiob2[ratio_pix22] * 0.001;  /* vs. / 1000. */
            if (rb2 > 1.0 || rb1 > 1.0 || rb2 < 0.1 || rb1 < 0.1)
            {
                slpratiob1[ratio_pix22] = 0;
                slpratiob2[ratio_pix22] = 0;
                slpratiob7[ratio_pix22] = 0;
                intratiob1[ratio_pix22] = 550;
                intratiob2[ratio_pix22] = 600;
                intratiob7[ratio_pix22] = 2000;
            }
            else if (sndwi[ratio_pix22] < 200)
            {
                slpratiob1[ratio_pix22] = 0;
                slpratiob2[ratio_pix22] = 0;
                slpratiob7[ratio_pix22] = 0;
                intratiob1[ratio_pix22] = ratiob1[ratio_pix22];
                intratiob2[ratio_pix22] = ratiob2[ratio_pix22];
                intratiob7[ratio_pix22] = ratiob7[ratio_pix22];
            }

            /* Compute the NDWI variables */
            ndwi_th1 = (andwi[ratio_pix11] + 2.0 *
                        sndwi[ratio_pix11]) * 0.001;
            ndwi_th2 = (andwi[ratio_pix11] - 2.0 *
                        sndwi[ratio_pix11]) * 0.001;

            /* Interpolate the slope/intercept for each band, and unscale */
            slpr11 = slpratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob1[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob1[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob1[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob1[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb1 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb1 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            slpr11 = slpratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob2[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob2[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob2[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob2[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb2 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb2 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            slpr11 = slpratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
            intr11 = intratiob7[ratio_pix11] * 0.001;  /* vs / 1000 */
            slpr12 = slpratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
            intr12 = intratiob7[ratio_pix12] * 0.001;  /* vs / 1000 */
            slpr21 = slpratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
            intr21 = intratiob7[ratio_pix21] * 0.001;  /* vs / 1000 */
            slpr22 = slpratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
            intr22 = intratiob7[ratio_pix22] * 0.001;  /* vs / 1000 */
            slprb7 = slpr11 * one_minus_u_x_one_minus_v +
                     slpr12 * one_minus_u_x_v +
                     slpr21 * u_x_one_minus_v +
                     slpr22 * u_x_v;
            intrb7 = intr11 * one_minus_u_x_one_minus_v +
                     intr12 * one_minus_u_x_v +
                     intr21 * u_x_one_minus_v +
                     intr22 * u_x_v;

            /* Calculate NDWI variables for the band ratios */
            xndwi = ((double) sband[SR_BAND5][curr_pix] -
                     (double) (sband[SR_BAND7][curr_pix] * 0.5)) /
                    ((double) sband[SR_BAND5][curr_pix] +
                     (double) (sband[SR_BAND7][curr_pix] * 0.5));

            if (xndwi > ndwi_th1)
                xndwi = ndwi_th1;
            if (xndwi < ndwi_th2)
                xndwi = ndwi_th2;

            /* Initialize the band ratios */
            for (ib = 0; ib < NSR_BANDS; ib++)
            {
                erelc[ib] = -1.0;
                troatm[ib] = 0.0;
            }

            /* Compute the band ratio */
            erelc[DN_BAND1] = (xndwi * slprb1 + intrb1);
            erelc[DN_BAND2] = (xndwi * slprb2 + intrb2);
            erelc[DN_BAND4] = 1.0;
            erelc[DN_BAND7] = (xndwi * slprb7 + intrb7);

            /* Retrieve the TOA reflectance values for the current pixel */
            troatm[DN_BAND1] = aerob1[curr_pix] * SCALE_FACTOR;
            troatm[DN_BAND2] = aerob2[curr_pix] * SCALE_FACTOR;
            troatm[DN_BAND4] = aerob4[curr_pix] * SCALE_FACTOR;
            troatm[DN_BAND7] = aerob7[curr_pix] * SCALE_FACTOR;

            /* Determine the solar and view angles for the current pixel */
            xtv = vza[curr_pix] * 0.01;
            xmuv = cos(xtv * DEG2RAD);
            xts = sza[curr_pix] * 0.01;
            xmus = cos(xts * DEG2RAD);
            xfi = saa[curr_pix] * 0.01 - vaa[curr_pix] * 0.01 ;
            cosxfi = cos(xfi * DEG2RAD);

            /* Retrieve the aerosol information for eps 1.0 */
            iband1 = DN_BAND4;
            iband3 = DN_BAND1;
            pres = tp[curr_pix];
            uoz = tozi[curr_pix];
            uwv = twvi[curr_pix];
            eps = 1.0;
            iaots = 0;
            retval = subaeroret (iband1, iband3, xts, xtv, xmus, xmuv, xfi,
                cosxfi, pres, uoz, uwv, erelc, troatm, tpres, aot550nm,
                rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, tauray,
                ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                oztransa, &raot, &residual, &iaots, eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing aerosol retrieval.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Save the data */
            eps1 = eps;
            residual1 = residual;
            sraot1 = raot;

            /* Retrieve the aerosol information for eps 1.75 */
            eps = 1.75;
            retval = subaeroret (iband1, iband3, xts, xtv, xmus, xmuv, xfi,
                cosxfi, pres, uoz, uwv, erelc, troatm, tpres, aot550nm,
                rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, tauray,
                ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                oztransa, &raot, &residual, &iaots, eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing aerosol retrieval.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Save the data */
            eps2 = eps;
            residual2 = residual;
            sraot2 = raot;

            /* Retrieve the aerosol information for eps 2.5 */
            eps = 2.5;
            retval = subaeroret (iband1, iband3, xts, xtv, xmus, xmuv, xfi,
                cosxfi, pres, uoz, uwv, erelc, troatm, tpres, aot550nm,
                rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, tauray,
                ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                oztransa, &raot, &residual, &iaots, eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing aerosol retrieval.");
                error_handler (true, FUNC_NAME, errmsg);
                exit (ERROR);
            }

            /* Save the data */
            eps3 = eps;
            residual3 = residual;
            sraot3 = raot;

            /* Find the eps that minimizes the residual */
            xa = (eps1 * eps1) - (eps3 * eps3);
            xd = (eps2 * eps2) - (eps3 * eps3);
            xb = eps1 - eps3;
            xe = eps2 - eps3;
            xc = residual1 - residual3;
            xf = residual2 - residual3;
            coefa = (xc*xe - xb*xf) / (xa*xe - xb*xd);
            coefb = (xa*xf - xc*xd) / (xa*xe - xb*xd);
            epsmin = -coefb / (2.0 * coefa);
            eps = epsmin;

            if (epsmin >= 1.0 && epsmin <= 2.5)
            {
                retval = subaeroret (iband1, iband3, xts, xtv, xmus, xmuv,
                    xfi, cosxfi, pres, uoz, uwv, erelc, troatm, tpres,
                    aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                    tts, indts, ttv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, &raot,
                    &residual, &iaots, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing aerosol retrieval.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
            }
            else
            {
                if (epsmin <= 1.0)
                {
                    eps = eps1;
                    residual = residual1;
                    raot = sraot1;
                }
                else if (epsmin >= 2.5)
                {
                    eps = eps3;
                    residual = residual3;
                    raot = sraot3;
                }
            }

            teps[center_pix] = eps;
            taero[center_pix] = raot;
            corf = raot / xmus;

            /* Check the model residual.  Corf represents aerosol impact.
               Test the quality of the aerosol inversion. */
            if (residual < (0.015 + 0.005 * corf + 0.10 * troatm[DN_BAND7]))
            {
                /* Test if band 5 makes sense */
                iband = DN_BAND5;
                rotoa = aerob5[curr_pix] * SCALE_FACTOR;
                raot550nm = raot;
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband, pres, tpres, aot550nm, rolutt,
                    transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                    normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                    &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian "
                        "atmospheric correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                ros5 = roslamb;

                /* Test if band 4 makes sense */
                iband = DN_BAND4;
                rotoa = aerob4[curr_pix] * SCALE_FACTOR;
                raot550nm = raot;
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband, pres, tpres, aot550nm, rolutt,
                    transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                    normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                    &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian "
                        "atmospheric correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                ros4 = roslamb;

                /* Use the NDVI to validate the reflectance values */
                if ((ros5 > 0.1) && ((ros5 - ros4) / (ros5 + ros4) > 0))
                {
                    /* Clear pixel with valid aerosol retrieval */
                    taero[center_pix] = raot;
                    ipflag[center_pix] |= (1 << IPFLAG_CLEAR);
                }
                else
                {
                    /* This is water and the retrieval needs to be redone */
                    taero[center_pix] = raot;
                    ipflag[center_pix] |= (1 << IPFLAG_WATER);
                }
            }
            else
            {
                /* Retrieval needs to be redone */
                taero[center_pix] = raot;
                ipflag[center_pix] |= (1 << IPFLAG_WATER);
            }


            /* Redo the aerosol retrieval if flagged as water, which is
               basically anything that isn't a good aerosol retrieval and
               flagged as clear. */
            if (btest (ipflag[center_pix], IPFLAG_WATER))
            {
//printf ("DEBUG:    IPFLAG_WATER pixel at line, sample: %d, %d\n", i, j);
                /* Reset variables */
                erelc[DN_BAND1] = 1.0;
                erelc[DN_BAND2] = -1.0;
                erelc[DN_BAND3] = -1.0;
                erelc[DN_BAND4] = 1.0;
                erelc[DN_BAND5] = 1.0;
                erelc[DN_BAND6] = -1.0;
                erelc[DN_BAND7] = 1.0;
                erelc[DN_BAND8] = -1.0;
                troatm[DN_BAND1] = aerob1[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND4] = aerob4[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND5] = aerob5[curr_pix] * SCALE_FACTOR;
                troatm[DN_BAND7] = aerob7[curr_pix] * SCALE_FACTOR;
                iband1 = DN_BAND4;
                iband3 = DN_BAND1;
                pres = tp[curr_pix];
                uoz = tozi[curr_pix];
                uwv = twvi[curr_pix];
                eps = 1.5;
                iaots = 0;

                /* Aerosol retrieval over water */
                retval = subaeroretwat (iband1, iband3, xts, xtv, xmus, xmuv,
                    xfi, cosxfi, pres, uoz, uwv, erelc, troatm, tpres,
                    aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
                    xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                    tts, indts, ttv, tauray, ogtransa1, ogtransb0, ogtransb1,
                    wvtransa, wvtransb, oztransa, &raot, &residual, &iaots,
                    eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing aerosol retrieval over "
                        "water.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                teps[center_pix] = eps;
                taero[center_pix] = raot;
                corf = raot / xmus;

                /* Test band 1 reflectance to eliminate negative */
                iband = DN_BAND1;
                rotoa = aerob1[curr_pix] * SCALE_FACTOR;
                raot550nm = raot;
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, iband, pres, tpres, aot550nm, rolutt,
                    transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                    normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                    &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian "
                        "atmospheric correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }
                ros1 = roslamb;

                /* Test the quality of the aerosol inversion */
                if (residual > (0.01 + 0.005 * corf) || ros1 < 0)
                {
                    /* Set the retrieval failed bit, and unset the water bit */
                    ipflag[center_pix] &= ~(1 << IPFLAG_WATER);
                    ipflag[center_pix] |= (1 << IPFLAG_RETRIEVAL_FAIL);
//printf ("DEBUG:    IPFLAG_RETRIEVAL_FAIL pixel at line, sample: %d, %d\n", i, j);
                }
                else
                {
                    teps[center_pix] = eps;
                    taero[center_pix] = raot;
                    corf = raot / xmus;
                }
            }  /* end if water */

            /* Reset the looping variables to the center of the aerosol window
               versus the actual non-fill/non-cloud pixel that was processed
               so that we get the correct center for the next aerosol window */
            i = center_line;
            j = center_samp;
            curr_pix = center_pix;
        }  /* end for j */
    }  /* end for i */

#ifndef _OPENMP
    /* update status */
    printf ("100%%\n");
    fflush (stdout);
#endif

    /* Done with the aerob* arrays */
    free (aerob1);  aerob1 = NULL;
    free (aerob2);  aerob2 = NULL;
    free (aerob4);  aerob4 = NULL;
    free (aerob5);  aerob5 = NULL;
    free (aerob7);  aerob7 = NULL;

    /* Done with the ratiob* arrays */
    free (andwi);  andwi = NULL;
    free (sndwi);  sndwi = NULL;
    free (ratiob1);  ratiob1 = NULL;
    free (ratiob2);  ratiob2 = NULL;
    free (ratiob7);  ratiob7 = NULL;
    free (intratiob1);  intratiob1 = NULL;
    free (intratiob2);  intratiob2 = NULL;
    free (intratiob7);  intratiob7 = NULL;
    free (slpratiob1);  slpratiob1 = NULL;
    free (slpratiob2);  slpratiob2 = NULL;
    free (slpratiob7);  slpratiob7 = NULL;

    /* Done with the DEM, water vapor, and ozone arrays */
    free (dem);  dem = NULL;
    free (wv);  wv = NULL;
    free (oz);  oz = NULL;

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Interpolate the aerosol pixels which could not be inverted */
mytime = time(NULL);
printf ("DEBUG: Starting aerosol window interpolation: %s\n", ctime(&mytime));
    if (aerosol_window_interp (qaband, ipflag, taero, teps, nlines, nsamps) !=
        SUCCESS)
    {
        sprintf (errmsg, "Performing interpolation of the NxN window aerosol "
            "values.");
        error_handler (true, FUNC_NAME, errmsg);
        exit (ERROR);
    }

#ifdef WRITE_TAERO
    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag2.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols2.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);
#endif

    /* Use the center of the aerosol windows to interpolate the remaining
       pixels in the window */
mytime = time(NULL);
printf ("DEBUG: Starting aerosol interpolation: %s\n", ctime(&mytime));
    printf ("Starting aerosol interpolation ...\n");
    aerosol_interp (xml_metadata, qaband, taero, nlines, nsamps);

#ifdef GAIL
    /* Fill in the aerosol windows with the value computed for the center of
       the window */
mytime = time(NULL);
printf ("DEBUG: Start fill aerosol windows: %s\n", ctime(&mytime));
    printf ("Fill in each aerosol window ...\n");
    for (i = HALF_AERO_WINDOW; i < nlines; i += AERO_WINDOW)
    {
        curr_pix = i * nsamps + HALF_AERO_WINDOW;
        for (j = HALF_AERO_WINDOW; j < nsamps;
             j += AERO_WINDOW, curr_pix += AERO_WINDOW)
        {
            /* If this pixel is fill, then the whole window is fill. Nothing
               to do so move on. Otherwise populate the aerosol window using
               the aerosol value from the center of the window. */
            if (btest (ipflag[curr_pix], IPFLAG_FILL))
                continue;
            else
            {
                populate_aerosol_window (taero, teps, ipflag, nlines, nsamps,
                    i, j);
            }
        }  /* end for j */
    }  /* end for i */
#endif

#ifdef GAIL
    /* Before the aerosol interpolation increase the interpolation to the
       neighboring pixels (5x5) */
mytime = time(NULL);
printf ("DEBUG: Setting the neighboring interpolation pixels: %s\n", ctime(&mytime));
    printf ("Setting the neighboring interpolation pixels ...\n"); 
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* If this pixel needs interpolated then look at the pixels in the
               surrounding window */
            if (btest (ipflag[curr_pix], IPFLAG_RETRIEVAL_FAIL))
            {
                /* Check the 5x5 window around the current pixel */
                for (k = i-2; k <= i+2; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;

                    win_pix = k * nsamps + j-2;
                    for (l = j-2; l <= j+2; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;

                        /* If this is a clear ipflag then set it to be
                           interpolated. Unset the clear bit. */
                        if (btest (ipflag[win_pix], IPFLAG_CLEAR))
                        {
                            ipflag[win_pix] &= ~(1 << IPFLAG_CLEAR);
                            ipflag[win_pix] |= (1 << IPFLAG_TMP_NEIGHBOR);
                        }
                    }  /* for l */
                }  /* for k */
            }  /* if ipflag == retrieval failed */
        }  /* for j */
    }  /* for i */

    /* Reset any pixels flagged as neighbors to interpolation failed */
    for (i = 0; i < nlines*nsamps; i++)
    {
        if (btest (ipflag[i], IPFLAG_TMP_NEIGHBOR))
        {
            /* Set the retrieval failed, but unset the temp neighbor */
            ipflag[i] &= ~(1 << IPFLAG_TMP_NEIGHBOR);
            ipflag[i] |= (1 << IPFLAG_RETRIEVAL_FAIL);
        }
    }

    /* Compute the average of the 11x11 window of aersols for each pixel */
    nbpixnf = 0;
    nbpixtot = 0;
    for (i = 0; i < nlines; i++)
    {
        curr_pix = i * nsamps;
        for (j = 0; j < nsamps; j++, curr_pix++)
        {
            /* Process non-fill pixels */
            smflag[curr_pix] = false;
            if (!level1_qa_is_fill (qaband[curr_pix]))
            {
                /* Initialize the variables */
                nbpixtot++;
                nbaeroavg = 0;
                taeroavg = 0.0;
                tepsavg = 0.0;

                /* Check the 11x11 window around the current pixel */
                for (k = i-5; k <= i+5; k++)
                {
                    /* Make sure the line is valid */
                    if (k < 0 || k >= nlines)
                        continue;

                    win_pix = k * nsamps + j-5;
                    for (l = j-5; l <= j+5; l++, win_pix++)
                    {
                        /* Make sure the sample is valid */
                        if (l < 0 || l >= nsamps)
                            continue;

                        /* If the pixel has clear retrieval of aerosols or
                           valid water aerosols then add it to the total sum */
                        if (btest (ipflag[win_pix], IPFLAG_CLEAR) ||
                            btest (ipflag[win_pix], IPFLAG_WATER))
                        {
                            nbaeroavg++;
                            taeroavg += taero[win_pix];
                            tepsavg += teps[win_pix];
                        }
                    }  /* for l */
                }  /* for k */

                /* If the number of clear/aerosol pixels in the window is high
                   enough, then compute the average for the window */
                if (nbaeroavg > 20)
                {
                    taeroavg = taeroavg / nbaeroavg;
                    tepsavg = tepsavg / nbaeroavg;
                    taeros[curr_pix] = taeroavg;
                    tepss[curr_pix] = tepsavg;
                    smflag[curr_pix] = true;
                }
                else
                {
                    /* Keep track of the number of pixels where the count
                       was not high enough */
                    nbpixnf++;
                }
            }  /* if not fill */
        }  /* for j */
    }  /* for i */

    /* Handle the case where none of the pixels were able to be averaged by
       setting default values.  Otherwise take a second pass through the
       pixels. */
mytime = time(NULL);
printf ("DEBUG: Starting second pass for aerosol interpolation: %s\n", ctime(&mytime));
    printf ("Second pass for aerosol interpolation ...\n");
    if (nbpixnf == nbpixtot)
    {
        /* Set defaults */
        for (i = 0; i < nlines*nsamps; i++)
        {
            /* If this is not a fill pixel then set default values */
            if (!level1_qa_is_fill (qaband[i]))
            {
                taeros[i] = 0.05;
                tepss[i] = 1.5;
                smflag[i] = true;
            }
        }  /* for i */
    }
    else
    {
        /* Second pass */
        prev_nbpixnf = 0;
        while (nbpixnf != 0)
        {
            /* If nothing was gained in the last loop then just use default
               values for the remaining pixels */
            if (nbpixnf == prev_nbpixnf)
            {
                for (i = 0; i < nlines; i++)
                {
                    curr_pix = i * nsamps;
                    for (j = 0; j < nsamps; j++, curr_pix++)
                    {
                        /* If this is not a fill pixel and the aerosol average
                           was not computed for this pixel */
                        if (qaband[curr_pix] != 1 && !smflag[curr_pix])
                        {
                            taeros[i] = 0.05;
                            tepss[i] = 1.5;
                            smflag[i] = true;
                        }  /* if qaband and smflag */
                    }  /* for j */
                }  /* for i */

                /* Break out of the while loop */
                break;
            }  /* if nbpixnf == prev_nbpixnf */

            prev_nbpixnf = nbpixnf;
            nbpixnf = 0;
            for (i = 0; i < nlines; i++)
            {
                curr_pix = i * nsamps;
                for (j = 0; j < nsamps; j++, curr_pix++)
                {
                    /* If this is not a fill pixel and the aerosol average was
                       not computed for this pixel */
                    if (!level1_qa_is_fill (qaband[curr_pix]) &&
                        !smflag[curr_pix])
                    {
                        nbaeroavg = 0;
                        tepsavg = 0.0;
                        taeroavg = 0.0;

                        /* Check the 11x11 window around the current pixel */
                        for (k = i-5; k <= i+5; k++)
                        {
                            /* Make sure the line is valid */
                            if (k < 0 || k >= nlines)
                                continue;
    
                            win_pix = k * nsamps + j-5;
                            for (l = j-5; l <= j+5; l++, win_pix++)
                            {
                                /* Make sure the sample is valid */
                                if (l < 0 || l >= nsamps)
                                    continue;
    
                                /* If the pixel has valid averages, then use
                                   them in this second pass */
                                if (smflag[win_pix])
                                {
                                    nbaeroavg++;
                                    taeroavg += taeros[win_pix];
                                    tepsavg += tepss[win_pix];
                                }
                            }  /* for l */
                        }  /* for k */

                        /* If at least one pixel in the window had valid
                           averages, then compute the average for the window
                           using the previously-computed window averages */
                        if (nbaeroavg > 0)
                        {
                            taeroavg = taeroavg / nbaeroavg;
                            tepsavg = tepsavg / nbaeroavg;
                            taeros[curr_pix] = taeroavg;
                            tepss[curr_pix] = tepsavg;
                            smflag[curr_pix] = true;
                        }
                        else
                        {
                            /* Keep track of the number of pixels where the
                               averages were not computed */
                            nbpixnf++;
                        }
                    }  /* if qaband and smflag */
                }  /* for j */
            }  /* for i */
        }  /* while nbpixnf != 0 */
    }  /* else */

    /* Fill in the pixels where the aerosol retrieval failed */
    for (i = 0; i < nlines*nsamps; i++)
    {
        /* Note the FORTRAN code also checks the cloud QA as to whether
           this is a cirrus or cloud pixel, however the cloud QA array has not
           yet been set for cirrus clouds or regular clouds.  Furthermore that
           code has been removed altogether.  Therefore the cloud checks have
           been removed from this if statement. */
        if (btest (ipflag[i], IPFLAG_RETRIEVAL_FAIL))
        {
            taero[i] = taeros[i];
            teps[i] = tepss[i];

            /* Set the aerosol interpolated, but unset the retrieval failed */
            ipflag[i] |= (1 << IPFLAG_INTERP);
            ipflag[i] &= ~(1 << IPFLAG_RETRIEVAL_FAIL);
        }
    }
#endif

#ifdef WRITE_TAERO
    /* Write the aerosol values for comparison with other algorithms */
    aero_fptr = fopen ("aerosols3.img", "w");
    fwrite (taero, nlines*nsamps, sizeof (float), aero_fptr);
    fclose (aero_fptr);

    /* Write the ipflag values for comparison with other algorithms */
    aero_fptr = fopen ("ipflag3.img", "w");
    fwrite (ipflag, nlines*nsamps, sizeof (uint8), aero_fptr);
    fclose (aero_fptr);

#endif

    /* Free the smflag, taeros, and tepss as temporary arrays */
    free (smflag);  smflag = NULL;
    free (taeros);  taeros = NULL;
    free (tepss);  tepss = NULL;

    /* Perform the second level of atmospheric correction using the aerosols */
mytime = time(NULL);
printf ("DEBUG: Starting second second level of atmospheric correction using aerosols: %s\n", ctime(&mytime));
    printf ("Performing atmospheric correction ...\n");
    /* 0 .. DN_BAND7 is the same as 0 .. SR_BAND7 here, since the pan band
       isn't spanned */
    for (ib = 0; ib <= DN_BAND7; ib++)
    {
        printf ("  Band %d\n", ib+1);
#ifdef _OPENMP
        #pragma omp parallel for private (i, j, curr_pix, xtv, xts, xmus, xmuv, xfi, cosxfi, rsurf, rotoa, raot550nm, eps, pres, uwv, uoz, retval, tmpf, roslamb, tgo, roatm, ttatmg, satm, xrorayp, next)
#endif
        for (i = 0; i < nlines; i++)
        {
            curr_pix = i * nsamps;
            for (j = 0; j < nsamps; j++, curr_pix++)
            {
                /* If this pixel is fill, then don't process */
                if (level1_qa_is_fill (qaband[curr_pix]))
                    continue;

                /* Determine the solar and view angles for the current pixel */
                xtv = vza[curr_pix] * 0.01;
                xmuv = cos(xtv * DEG2RAD);
                xts = sza[curr_pix] * 0.01;
                xmus = cos(xts * DEG2RAD);
                xfi = saa[curr_pix] * 0.01 - vaa[curr_pix] * 0.01;
                cosxfi = cos(xfi * DEG2RAD);

                /* Correct all pixels */
                rsurf = sband[ib][curr_pix] * SCALE_FACTOR;
                rotoa = (rsurf * bttatmg[ib] / (1.0 - bsatm[ib] * rsurf) +
                    broatm[ib]) * btgo[ib];
                raot550nm = taero[curr_pix];
                eps = teps[curr_pix];
                pres = tp[curr_pix];
                uwv = twvi[curr_pix];
                uoz = tozi[curr_pix];
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt,
                    transt, xtsstep, xtsmin, xtvstep, xtvmin, sphalbt,
                    normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
                    ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                    ogtransb1, wvtransa, wvtransb, oztransa, rotoa,
                    &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
                    &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian "
                        "atmospheric correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    exit (ERROR);
                }

                /* If this is the coastal aerosol band then set the aerosol
                   bits in the QA band */
                if (ib == DN_BAND1)
                {
                    /* Set up aerosol QA bits */
                    tmpf = fabs (rsurf - roslamb);
                    if (tmpf <= 0.015)
                    {  /* Set the first aerosol bit (low aerosols) */
                        ipflag[curr_pix] |= (1 << AERO1_QA);
                    }
                    else
                    {
                        if (tmpf < 0.03)
                        {  /* Set the second aerosol bit (average aerosols) */
                            ipflag[curr_pix] |= (1 << AERO2_QA);
                        }
                        else
                        {  /* Set both aerosol bits (high aerosols) */
                            ipflag[curr_pix] |= (1 << AERO1_QA);
                            ipflag[curr_pix] |= (1 << AERO2_QA);
                        }
                    }
                }  /* end if this is the coastal aerosol band */

                /* Save the scaled surface reflectance value, but make sure it
                   falls within the defined valid range. */
                roslamb = roslamb * MULT_FACTOR;  /* scale the value */
                if (roslamb < MIN_VALID)
                    sband[ib][curr_pix] = MIN_VALID;
                else if (roslamb > MAX_VALID)
                    sband[ib][curr_pix] = MAX_VALID;
                else
                    sband[ib][curr_pix] = (int) (roundf (roslamb));
            }  /* end for j */
        }  /* end for i */
    }  /* end for ib */

    /* Free memory for arrays no longer needed */
    free (twvi);
    free (tozi);
    free (tp);
    free (taero);
    free (teps);
 
    /* Write the data to the output file */
mytime = time(NULL);
printf ("DEBUG: Writing SR results to output files: %s\n", ctime(&mytime));
    printf ("Writing surface reflectance corrected data to the output "
        "files ...\n");

    /* Open the output file */
    sr_output = open_output (xml_metadata, input, OUTPUT_SR);
    if (sr_output == NULL)
    {   /* error message already printed */
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Loop through the reflectance bands and write the data */
    for (ib = 0; ib <= DN_BAND7; ib++)
    {
        printf ("  Band %d: %s\n", ib+1,
            sr_output->metadata.band[ib].file_name);
        if (put_output_lines (sr_output, sband[ib], ib, 0, nlines,
            sizeof (int16)) != SUCCESS)
        {
            sprintf (errmsg, "Writing output data for band %d", ib);
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Create the ENVI header file this band */
        if (create_envi_struct (&sr_output->metadata.band[ib],
            &xml_metadata->global, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Creating ENVI header structure.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        /* Write the ENVI header */
        strcpy (envi_file, sr_output->metadata.band[ib].file_name);
        cptr = strchr (envi_file, '.');
        strcpy (cptr, ".hdr");
        if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        {
            sprintf (errmsg, "Writing ENVI header file.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }
    }

    /* Append the surface reflectance bands (1-7) to the XML file */
    if (append_metadata (7, sr_output->metadata.band, xml_infile) !=
        SUCCESS)
    {
        sprintf (errmsg, "Appending surface reflectance bands to the "
            "XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the aerosol QA band */
    printf ("  Band %d: %s\n", SR_AEROSOL+1,
            sr_output->metadata.band[SR_AEROSOL].file_name);
    if (put_output_lines (sr_output, ipflag, SR_AEROSOL, 0, nlines,
        sizeof (uint8)) != SUCCESS)
    {
        sprintf (errmsg, "Writing aerosol QA output data");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Free memory for ipflag data */
    free (ipflag);

    /* Create the ENVI header for the aerosol QA band */
    if (create_envi_struct (&sr_output->metadata.band[SR_AEROSOL],
        &xml_metadata->global, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Creating ENVI header structure.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Write the ENVI header */
    strcpy (envi_file, sr_output->metadata.band[SR_AEROSOL].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errmsg, "Writing ENVI header file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Append the aerosol QA band to the XML file */
    if (append_metadata (1, &sr_output->metadata.band[SR_AEROSOL],
        xml_infile) != SUCCESS)
    {
        sprintf (errmsg, "Appending aerosol QA band to XML file.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Close the output surface reflectance products */
    close_output (sr_output, OUTPUT_SR);
    free_output (sr_output);

    /* Free the spatial mapping pointer */
    free (space);

    /* Free the data arrays */
    free (rolutt);
    free (transt);
    free (sphalbt);
    free (normext);
    free (tsmax);
    free (tsmin);
    free (nbfic);
    free (nbfi);
    free (ttv);

    /* Successful completion */
mytime = time(NULL);
printf ("DEBUG: All DONE: %s\n", ctime(&mytime));
    return (SUCCESS);
}


/******************************************************************************
MODULE:  init_sr_refl

PURPOSE:  Initialization for the atmospheric corrections.  Initialization for
look up tables, auxiliary data, mapping, and geolocation information is used
for the surface reflectance correction.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error initializing the atmospheric parameters
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. The view angle is set to 0.0 and this never changes.
2. The DEM is used to calculate the surface pressure.
******************************************************************************/
int init_sr_refl
(
    int nlines,         /* I: number of lines in reflectance, thermal bands */
    int nsamps,         /* I: number of samps in reflectance, thermal bands */
    Input_t *input,     /* I: input structure for the Landsat product */
    Geoloc_t *space,    /* I: structure for geolocation information */
    char *anglehdf,     /* I: angle HDF filename */
    char *intrefnm,     /* I: intrinsic reflectance filename */
    char *transmnm,     /* I: transmission filename */
    char *spheranm,     /* I: spherical albedo filename */
    char *cmgdemnm,     /* I: climate modeling grid DEM filename */
    char *rationm,      /* I: ratio averages filename */
    char *auxnm,        /* I: auxiliary filename for ozone and water vapor */
    float *eps,         /* O: angstrom coefficient */
    int *iaots,         /* O: index for AOTs */ 
    float *xtv,         /* O: observation zenith angle (deg) */
    float *xmuv,        /* O: cosine of observation zenith angle */
    float *xfi,         /* O: azimuthal difference between sun and
                              observation (deg) */
    float *cosxfi,      /* O: cosine of azimuthal difference */
    float *raot550nm,   /* O: nearest value of AOT */
    float *pres,        /* O: surface pressure */
    float *uoz,         /* O: total column ozone */
    float *uwv,         /* O: total column water vapor (precipital water
                              vapor) */
    float *xtsstep,     /* O: solar zenith step value */
    float *xtsmin,      /* O: minimum solar zenith value */
    float *xtvstep,     /* O: observation step value */
    float *xtvmin,      /* O: minimum observation value */
    float *tsmax,       /* O: maximum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *tsmin,       /* O: minimum scattering angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float tts[22],      /* O: sun angle table */
    float *ttv,         /* O: view angle table
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int32 indts[22],    /* O: index for the sun angle table */
    float *rolutt,      /* O: intrinsic reflectance table
                          [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSOLAR_VALS] */
    float *transt,      /* O: transmission table 
                      [NSR_BANDS x NPRES_VALS x NAOT_VALS x NSUN_ANGLE_VALS] */
    float *sphalbt,     /* O: spherical albedo table 
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *normext,     /* O: aerosol extinction coefficient at the current
                              wavelength (normalized at 550nm)
                              [NSR_BANDS x NPRES_VALS x NAOT_VALS] */
    float *nbfic,       /* O: communitive number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    float *nbfi,        /* O: number of azimuth angles
                              [NVIEW_ZEN_VALS x NSOLAR_ZEN_VALS] */
    int16 *dem,         /* O: CMG DEM data array [DEM_NBLAT x DEM_NBLON] */
    int16 *andwi,       /* O: avg NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *sndwi,       /* O: standard NDWI [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob1,     /* O: mean band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob2,     /* O: mean band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *ratiob7,     /* O: mean band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob1,  /* O: integer band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob2,  /* O: integer band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *intratiob7,  /* O: integer band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob1,  /* O: slope band1 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob2,  /* O: slope band2 ratio [RATIO_NBLAT x RATIO_NBLON] */
    int16 *slpratiob7,  /* O: slope band7 ratio [RATIO_NBLAT x RATIO_NBLON] */
    uint16 *wv,         /* O: water vapor values [CMG_NBLAT x CMG_NBLON] */
    uint8 *oz           /* O: ozone values [CMG_NBLAT x CMG_NBLON] */
)
{
    char errmsg[STR_SIZE];                   /* error message */
    char FUNC_NAME[] = "init_sr_refl";       /* function name */
    int retval;          /* return status */
    int lcmg, scmg;      /* line/sample index for the CMG */
    int cmg_pix;         /* pixel location in the CMG array for [lcmg][scmg] */
    int dem_pix;         /* pixel location in the DEM array for [lcmg][scmg] */
    float xcmg, ycmg;    /* x/y location for CMG */

    /* Vars for forward/inverse mapping space */
    Img_coord_float_t img;        /* coordinate in line/sample space */
    Geo_coord_t geo;              /* coordinate in lat/long space */
    float center_lat, center_lon; /* lat/long for scene center */

    /* Initialize the look up tables */
    *eps = 1.0;
    *iaots = 0;
    *xtv = 0.0;
    *xmuv = cos (*xtv * DEG2RAD);
    *xfi = 0.0;
    *cosxfi = cos (*xfi * DEG2RAD);
    *xtsmin = 0;
    *xtsstep = 4.0;
    *xtvmin = 2.84090;
    *xtvstep = 6.52107 - *xtvmin;
    retval = readluts (tsmax, tsmin, ttv, tts, nbfic, nbfi, indts, rolutt,
        transt, sphalbt, normext, *xtsstep, *xtsmin, anglehdf,
        intrefnm, transmnm, spheranm);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the LUTs");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    printf ("The LUTs for urban clean case v2.0 have been read.  We can "
        "now perform atmospheric correction.\n");

    /* Read the auxiliary data files used as input to the reflectance
       calculations */
    retval = read_auxiliary_files (cmgdemnm, rationm, auxnm, dem, andwi, sndwi,
        ratiob1, ratiob2, ratiob7, intratiob1, intratiob2, intratiob7,
        slpratiob1, slpratiob2, slpratiob7, wv, oz);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Reading the auxiliary files");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Getting parameters for atmospheric correction */
    /* Update to get the parameter of the scene center */
    *raot550nm = 0.12;
    *pres = 1013.0;
    *uoz = 0.30;
    *uwv = 0.5;

    /* Use scene center (and center of the pixel) to compute atmospheric
       parameters */
    img.l = (int) (nlines * 0.5) - 0.5;
    img.s = (int) (nsamps * 0.5) + 0.5;
    img.is_fill = false;
    if (!from_space (space, &img, &geo))
    {
        sprintf (errmsg, "Mapping scene center to geolocation coords");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }
    center_lat = geo.lat * RAD2DEG;
    center_lon = geo.lon * RAD2DEG;
    printf ("Scene center line/sample: %f, %f\n", img.l, img.s);
    printf ("Scene center lat/long: %f, %f\n", center_lat, center_lon);

    /* Use the scene center lat/long to determine the line/sample in the
       CMG-related lookup tables, using the center of the UL pixel.
       Negative latitude values should be the largest line values in the CMG
       grid.  Negative longitude values should be the smallest sample values
       in the CMG grid. */
    ycmg = (89.975 - center_lat) * 20.0;    /* vs / 0.05 */
    xcmg = (179.975 + center_lon) * 20.0;   /* vs / 0.05 */
    lcmg = (int) roundf (ycmg);
    scmg = (int) roundf (xcmg);

    /* Handle the edges of the lat/long values in the CMG grid */
    if (lcmg < 0)
        lcmg = 0;
    else if (lcmg >= CMG_NBLAT)
        lcmg = CMG_NBLAT;

    if (scmg < 0)
        scmg = 0;
    else if (scmg >= CMG_NBLON)
        scmg = CMG_NBLON;

    cmg_pix = lcmg * CMG_NBLON + scmg;
    if (wv[cmg_pix] != 0)
        *uwv = wv[cmg_pix] / 200.0;
    else
        *uwv = 0.5;

    if (oz[cmg_pix] != 0)
        *uoz = oz[cmg_pix] / 400.0;
    else
        *uoz = 0.3;

    dem_pix = lcmg * DEM_NBLON + scmg;
    if (dem[dem_pix] != -9999)
        *pres = 1013.0 * exp (-dem[dem_pix] * ONE_DIV_8500);
    else
        *pres = 1013.0;
    *raot550nm = 0.05;

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  is_cloud

PURPOSE:  Determines if the pixel is a cloud (cloud, cloud shadow, or cirrus
cloud).  The Level-1 QA band is used.  A confidence of high for any of the
three QA types will result in the pixel being flagged as cloudy.

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           Pixel is not cloud
true            Pixel is cloud

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool is_cloud
(
    uint16_t l1_qa_pix      /* I: Level-1 QA value for current pixel */
)
{
    /* If the confidence level is high for cloud, cloud shadow, or cirrus, then
       flag this as a cloud */
    if (level1_qa_cloud_confidence (l1_qa_pix) == 3 ||
        level1_qa_cloud_shadow_confidence (l1_qa_pix) == 3 ||
        level1_qa_cirrus_confidence (l1_qa_pix) == 3)
        return (true);
    else
        return (false);
}


/******************************************************************************
MODULE:  is_water

PURPOSE:  Determines if the pixel is water.  The NDVI is used to determine if
this is a water pixel.

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           Pixel is not water
true            Pixel is water

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool is_water
(
    int16 band4_pix,     /* I: Band 4 reflectance for current pixel */
    int16 band5_pix      /* I: Band 5 reflectance for current pixel */
)
{
    double ndvi;             /* use NDVI for flagging water pixels */

    /* Calculate NDVI and flag water pixels */
    if (band5_pix < 100)
        ndvi = -0.01;
    else
        ndvi = ((double) band5_pix - (double) band4_pix) /
               ((double) band5_pix + (double) band4_pix);

    /* If the NDVI is low, then flag this as a water pixel */
    if (ndvi < 0.01)
        return (true);
    else
        return (false);
}


/******************************************************************************
MODULE:  find_closest_non_fill

PURPOSE:  Finds the closest non-fill pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool find_closest_non_fill
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-fill pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-fill pix in aerosol window */
)
{
    int i;                   /* looping variable for pixels */
    int line, samp;          /* looping variables for lines and samples */
    int aero_window;         /* looping variabel for the aerosol window */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        for (line = center_line - aero_window;
             line <= center_line + aero_window; line++)
        {
            /* Make sure the line is valid */
            if (line < 0 || line >= nlines)
                continue;

            i = line * nsamps;
            for (samp = center_samp - aero_window;
                 samp <= center_samp + aero_window; samp++)
            {
                /* Make sure the sample is valid */
                if (samp < 0 || samp >= nsamps)
                    continue;

                /* If this pixel is not fill, then mark it as the closest
                   non-fill pixel and return */
                if (!level1_qa_is_fill (qaband[i + samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return (true);
                }
            }
        }
    }

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  find_closest_non_cloud

PURPOSE:  Finds the closest non-cloud pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool find_closest_non_cloud
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
)
{
    int i;                   /* looping variable for pixels */
    int line, samp;          /* looping variables for lines and samples */
    int aero_window;         /* looping variabel for the aerosol window */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        for (line = center_line - aero_window;
             line <= center_line + aero_window; line++)
        {
            /* Make sure the line is valid */
            if (line < 0 || line >= nlines)
                continue;

            i = line * nsamps;
            for (samp = center_samp - aero_window;
                 samp <= center_samp + aero_window; samp++)
            {
                /* Make sure the sample is valid */
                if (samp < 0 || samp >= nsamps)
                    continue;

                /* If this pixel is not fill and is not cloud, then mark it as
                   the closest non-cloud pixel and return */
                if (!level1_qa_is_fill (qaband[i + samp]) &&
                    !is_cloud (qaband[i + samp]))
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return (true);
                }
            }
        }
    }

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  find_closest_non_water

PURPOSE:  Finds the closest non-water pixel in the aerosol window

RETURN VALUE:
Type = boolean
Value           Description
-----           -----------
false           No pixel found
true            Pixel found

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
bool find_closest_non_water
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int16 **sband,     /* I: input surface reflectance, nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    int *nearest_line, /* O: line for nearest non-cloud pix in aerosol window */
    int *nearest_samp  /* O: samp for nearest non-cloud pix in aerosol window */
)
{
    int i;                   /* looping variable for pixels */
    int line, samp;          /* looping variables for lines and samples */
    int aero_window;         /* looping variabel for the aerosol window */
    double ndvi;             /* use NDVI for flagging water pixels */

    /* Loop around the center pixel, moving outward with each loop, searching
       for a pixel that is not of the QA type specified and is not fill */
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        for (line = center_line - aero_window;
             line <= center_line + aero_window; line++)
        {
            /* Make sure the line is valid */
            if (line < 0 || line >= nlines)
                continue;

            i = line * nsamps;
            for (samp = center_samp - aero_window;
                 samp <= center_samp + aero_window; samp++)
            {
                /* Make sure the sample is valid */
                if (samp < 0 || samp >= nsamps)
                    continue;

                /* Calculate NDVI and flag water pixels */
                if (sband[SR_BAND5][i + samp] < 100)
                    ndvi = -0.01;
                else
                    ndvi = ((double) sband[SR_BAND5][i + samp] -
                            (double) sband[SR_BAND4][i + samp]) /
                           ((double) sband[SR_BAND5][i + samp] +
                            (double) sband[SR_BAND4][i + samp]);

                /* If this pixel is not fill, is not cloud, and is not water,
                   then mark it as the closest non-water pixel and return.
                   The NDVI will be used to flag water.  If ndvi < 0.01 then
                   the pixel is water. */
                if (!level1_qa_is_fill (qaband[i + samp]) &&
                    !is_cloud (qaband[i + samp]) && ndvi >= 0.01)
                {
                    *nearest_line = line;
                    *nearest_samp = samp;
                    return (true);
                }
            }
        }
    }

    /* No pixel was found that met the criteria */
    return (false);
}


/******************************************************************************
MODULE:  mask_aero_window

PURPOSE:  Masks the current pixel's quick use aerosol window for fill, cloud,
and water pixels.

RETURN VALUE: N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void mask_aero_window
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    int16 **sband,     /* I: input surface reflectance */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp,   /* I: sample for the center of the aerosol window */
    bool *quick_qa     /* O: quick QA for the current aerosol window,
                             AERO_WINDOW x AERO_WINDOW
                             (true=not clear, false=clear) */
)
{
    int curr_pix;            /* looping variable for current pixel in the
                                level-1 QA */
    int curr_qa_pix;         /* looping variable for current quick QA pixel */
    int line, samp;          /* looping variables for lines and samples */
    double ndvi;             /* use NDVI for flagging water pixels */

    /* Initialize the quick QA window to not clear, which includes pixels
       that go beyond the scene boundaries */
    for (curr_qa_pix = 0; curr_qa_pix < AERO_WINDOW * AERO_WINDOW;
         curr_qa_pix++)
        quick_qa[curr_qa_pix] = true;

    /* Loop around the current aerosol window flagging fill, cloudy, and water
       pixels */
    curr_qa_pix = 0;
    for (line = center_line - HALF_AERO_WINDOW;
         line <= center_line + HALF_AERO_WINDOW; line++)
    {
        /* Make sure the line is valid */
        if (line < 0 || line >= nlines)
            continue;

        curr_pix = line * nsamps + center_samp - HALF_AERO_WINDOW;
        for (samp = center_samp - HALF_AERO_WINDOW;
             samp <= center_samp + HALF_AERO_WINDOW;
             samp++, curr_pix++, curr_qa_pix++)
        {
            /* Make sure the sample is valid */
            if (samp < 0 || samp >= nsamps)
                continue;

            /* Calculate NDVI and flag water pixels */
            if (sband[SR_BAND5][curr_pix] < 100)
                ndvi = -0.01;
            else
                ndvi = ((double) sband[SR_BAND5][curr_pix] -
                        (double) sband[SR_BAND4][curr_pix]) /
                       ((double) sband[SR_BAND5][curr_pix] +
                        (double) sband[SR_BAND4][curr_pix]);

            /* If this pixel is not fill, is not cloud, and is not water, then
               mark it as clear.  The NDVI will be used to flag water.  If 
               ndvi < 0.01 then the pixel is water. */
            if (!level1_qa_is_fill (qaband[curr_pix]) &&
                !is_cloud (qaband[curr_pix]) && ndvi >= 0.01)
            { /* pixel is clear */
                quick_qa[curr_qa_pix] = false;
            }
        }  /* for samp */
    }  /* for line */
}


/******************************************************************************
MODULE:  populate_aerosol_window

PURPOSE:  Populates the entire aerosol window (for taero and teps) using the
value from the center of the aerosol window.

RETURN VALUE: N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void populate_aerosol_window
(
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps */
    int nlines,        /* I: number of lines in QA band */
    int nsamps,        /* I: number of samps in QA band */
    int center_line,   /* I: line for the center of the aerosol window */
    int center_samp    /* I: sample for the center of the aerosol window */
)
{
    int curr_pix;        /* current pixel in 1D arrays of nlines * nsamps */
    int center_pix;      /* current pixel in 1D arrays of nlines * nsamps for
                            the center of the aerosol window */
    int line, samp;      /* looping variables for lines and samples */
    int aero_window;     /* looping variabel for the aerosol window */

    /* Loop around the center pixel, moving outward with each loop, and populate
       the aerosol values for the entire aerosol window using the values from
       the center of the window */
    center_pix = center_line * nsamps + center_samp;
    for (aero_window = 1; aero_window <= HALF_AERO_WINDOW; aero_window++)
    {
        for (line = center_line - aero_window;
             line <= center_line + aero_window; line++)
        {
            /* Make sure the line is valid */
            if (line < 0 || line >= nlines)
                continue;

            curr_pix = line * nsamps;
            for (samp = center_samp - aero_window;
                 samp <= center_samp + aero_window; samp++)
            {
                /* Make sure the sample is valid */
                if (samp < 0 || samp >= nsamps)
                    continue;

                /* Populate the pixel from the center */
                teps[curr_pix + samp] = teps[center_pix];
                taero[curr_pix + samp] = taero[center_pix];
                ipflag[curr_pix + samp] = ipflag[center_pix];
            }
        }
    }
}
