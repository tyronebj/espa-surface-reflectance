/*****************************************************************************
FILE: subaeroret.c
  
PURPOSE: Contains functions for handling the atmosperic corrections.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "lut_subr.h"

/******************************************************************************
MODULE:  subaeroret

PURPOSE:  Main driver for the atmospheric correction.  This subroutine reads
the lookup table (LUT) and performs the atmospheric corrections.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int subaeroret
(
    int iband1,                      /* I: band 1 index (0-based) */
    int iband3,                      /* I: band 3 index (0-based) */
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xmus,                      /* I: cosine of solar zenith angle */
    float xmuv,                      /* I: cosine of observation zenith angle */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float cosxfi,                    /* I: cosine of azimuthal difference */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float erelc[NSR_BANDS],          /* I: band ratio variable */
    float troatm[NSR_BANDS],         /* I: atmospheric reflectance table */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: ????
                                           [NSR_BANDS][7][22] */
    float **tsmax,                   /* I: maximum scattering angle table
                                           [20][22] */
    float **tsmin,                   /* I: minimum scattering angle table
                                           [20][22] */
    float **nbfic,                   /* I: communitive number of azimuth angles
                                           [20][22] */
    float **nbfi,                    /* I: number of azimuth anglesi [20][22] */
    float tts[22],                   /* I: sun angle table */
    int32 indts[22],
    float **ttv,                     /* I: view angle table [20][22] */
    float tauray[NSR_BANDS],         /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS],     /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],      /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],      /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],      /* I: ozone transmission coeff */
    float *raot,                     /* O: AOT reflectance */
    float *residual,                 /* O: model residual */
    int *iaots,                      /* I/O: AOT index that is passed in and
                                             out for multiple calls (0-based) */
    float eps                        /* I: angstroem coefficient; spectral
                                           dependency of the AOT */
)
{
    char FUNC_NAME[] = "subaeroret";   /* function name */
    char errmsg[STR_SIZE];       /* error message */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int retval;             /* function return value */
    int ib;                 /* band index */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1, ros3;      /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    float raotsaved;        /* save the raot value */
    float next;             /* ???? */
    float tgo;              /* other gaseous transmittance */
    float roatm;            /* atmospheric intrinsic reflectance */
    float ttatmg;           /* total atmospheric transmission */
    float satm;             /* spherical albedo */
    float xrorayp;          /* reflectance of the atmosphere due to molecular
                               (Rayleigh) scattering */
    double residual1, residual2;  /* residuals for storing and comparing */
    double residualm;       /* local model residual */
    int nbval;              /* number of values meeting criteria */
    bool testth;            /* surface reflectance test variable */
    double xa, xb, xc, xd, xe, xf;  /* AOT ratio values */
    double coefa, coefb;    /* AOT ratio coefficients */
    double raotmin;         /* minimum AOT ratio */
    int iaot1, iaot2;       /* AOT indices (0-based) */
    float tth[NSR_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0, 1.0e-04,
                            0.0}; /* constant values for comparing against the
                                     surface reflectance */

    /* Correct band 3 and band 1 with increasing AOT (using pre till ratio is
       equal to erelc[2]) */
    iaot = *iaots;
    residual1 = 2000.0;
    residual2 = 1000.0;
    iaot2 = 0;
    iaot1 = 0;
    raot2 = 1.0e-06;
    raot1 = 0.0001;
    ros1 = 1.0;
    ros3 = 1.0;
    raot550nm = aot550nm[iaot];
    testth = false;

    /* Atmospheric correction for band 1 */
    ib = iband1;
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
        &next, eps);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (roslamb - tth[iband1] < 0.0)
        testth = true;
    ros1 = roslamb;

    /* Atmospheric correction for each band */
    nbval = 0;
    *residual = 0.0;
    for (ib = DN_BAND1; ib < DN_BAND8; ib++)
    {
        /* Don't reprocess iband1 */
        if ((erelc[ib] > 0.0) && (ib != iband1))
        {
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
                ib, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                ogtransb1, wvtransa, wvtransb, oztransa, troatm[ib], &roslamb,
                &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next, eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            if (roslamb - tth[ib] < 0.0)
                testth = true;
            *residual += (roslamb - erelc[ib] * ros1) *
                         (roslamb - erelc[ib] * ros1);
            nbval++;
        }
    }
    *residual = sqrt (*residual) / nbval;

    /* Loop until we converge on a solution */
    iaot++;
    while ((iaot < 22) && (*residual < residual1) && (!testth))
    {
        /* Reset variables for this loop */
        residual2 = residual1;
        iaot2 = iaot1;
        raot2 = raot1;
        residual1 = *residual;
        raot1 = raot550nm;
        iaot1 = iaot;
        raot550nm = aot550nm[iaot];

        /* Atmospheric correction for band 1 */
        ib = iband1;
        testth = false;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
            pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
            &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        *residual = 0.0;
        for (ib = DN_BAND1; ib < DN_BAND8; ib++)
        {
            /* Don't reprocess iband1 */
            if ((erelc[ib] > 0.0) && (ib != iband1))
            {
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
    
                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                *residual += (roslamb - erelc[ib] * ros1) *
                             (roslamb - erelc[ib] * ros1);
                nbval++;
            }
        }
        *residual = sqrt (*residual) / nbval;

        /* Move to the next AOT index */
        iaot++;
    }  /* while aot */

    /* If a minimum local was not reached for raot1, then just use the
       raot550nm value.  Otherwise continue to refine the raot. */
    if (iaot == 1)
    {
        *raot = raot550nm;
    }
    else
    {
        /* Refine the AOT ratio */
        *raot = raot550nm;
        raotsaved = *raot;
        xa = (raot1 * raot1) - (*raot * *raot);
        xd = (raot2 * raot2) - (*raot * *raot);
        xb = raot1 - *raot;
        xe = raot2 - *raot;
        xc = residual1 - *residual;
        xf = residual2 - *residual;
        coefa = (xc * xe - xb * xf) / (xa * xe - xb * xd);
        coefb = (xa * xf - xc * xd) / (xa * xe - xb * xd);
        raotmin = -coefb / (2.0 * coefa);

        /* Validate the min AOT ratio */
        if (raotmin < 0.01 || raotmin > 4.0)
            raotmin = *raot;

        /* Atmospheric correction for band 1 */
        raot550nm = raotmin;
        ib = iband1;
        testth = false;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
            pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
            &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        residualm = 0.0;
        for (ib = DN_BAND1; ib < DN_BAND8; ib++)
        {
            /* Don't reprocess iband1 */
            if ((erelc[ib] > 0.0) && (ib != iband1))
            {
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                residualm += (roslamb - erelc[ib] * ros1) *
                             (roslamb - erelc[ib] * ros1);
                nbval++;
            }
        }

        residualm = sqrt (residualm) / nbval;
        *raot = raot550nm;

        /* Check the residuals and reset the AOT ratio */
        if (residualm > *residual)
        {
            residualm = *residual;
            *raot = raotsaved;
        }
        if (residualm > residual1)
        {
            residualm = residual1;
            *raot = raot1;
        }
        if (residualm > residual2)
        {
            residualm = residual2;
            *raot = raot2;
        }

        *residual = residualm;
        *iaots = MAX ((iaot2 - 3), 0);
    }

    /* Successful completion */
    return (SUCCESS);
}


/******************************************************************************
MODULE:  subaeroretwat

PURPOSE:  Main driver for the atmospheric correction over water.  This
subroutine reads the lookup table (LUT) and performs the atmospheric
corrections. The residual computations are different in this case then they
are for non-water pixels.
TODO -- consider combining with subaeroret and just pass in a flag for handling
water pixels.

RETURN VALUE:
Type = int
Value          Description
-----          -----------
ERROR          Error occurred reading the LUT or doing the correction
SUCCESS        Successful completion

NOTES:
******************************************************************************/
int subaeroretwat
(
    int iband1,                      /* I: band 1 index (0-based) */
    int iband3,                      /* I: band 3 index (0-based) */
    float xts,                       /* I: solar zenith angle (deg) */
    float xtv,                       /* I: observation zenith angle (deg) */
    float xmus,                      /* I: cosine of solar zenith angle */
    float xmuv,                      /* I: cosine of observation zenith angle */
    float xfi,                       /* I: azimuthal difference between sun and
                                           observation (deg) */
    float cosxfi,                    /* I: cosine of azimuthal difference */
    float pres,                      /* I: surface pressure */
    float uoz,                       /* I: total column ozone */
    float uwv,                       /* I: total column water vapor (precipital
                                           water vapor) */
    float erelc[NSR_BANDS],          /* I: band ratio variable */
    float troatm[NSR_BANDS],         /* I: atmospheric reflectance table */
    float tpres[7],                  /* I: surface pressure table */
    float aot550nm[22],              /* I: AOT look-up table */
    float ****rolutt,                /* I: intrinsic reflectance table
                                           [NSR_BANDS][7][22][8000] */
    float ****transt,                /* I: transmission table
                                           [NSR_BANDS][7][22][22] */
    float xtsstep,                   /* I: solar zenith step value */
    float xtsmin,                    /* I: minimum solar zenith value */
    float xtvstep,                   /* I: observation step value */
    float xtvmin,                    /* I: minimum observation value */
    float ***sphalbt,                /* I: spherical albedo table
                                           [NSR_BANDS][7][22] */
    float ***normext,                /* I: ????
                                           [NSR_BANDS][7][22] */
    float **tsmax,                   /* I: maximum scattering angle table
                                           [20][22] */
    float **tsmin,                   /* I: minimum scattering angle table
                                           [20][22] */
    float **nbfic,                   /* I: communitive number of azimuth angles
                                           [20][22] */
    float **nbfi,                    /* I: number of azimuth anglesi [20][22] */
    float tts[22],                   /* I: sun angle table */
    int32 indts[22],
    float **ttv,                     /* I: view angle table [20][22] */
    float tauray[NSR_BANDS],         /* I: molecular optical thickness coeff */
    double ogtransa1[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb0[NSR_BANDS],     /* I: other gases transmission coeff */
    double ogtransb1[NSR_BANDS],     /* I: other gases transmission coeff */
    double wvtransa[NSR_BANDS],      /* I: water vapor transmission coeff */
    double wvtransb[NSR_BANDS],      /* I: water vapor transmission coeff */
    double oztransa[NSR_BANDS],      /* I: ozone transmission coeff */
    float *raot,                     /* O: AOT reflectance */
    float *residual,                 /* O: model residual */
    int *iaots,                      /* I/O: AOT index that is passed in and
                                             out for multiple calls (0-based) */
    float eps                        /* I: angstroem coefficient; spectral
                                           dependency of the AOT */
)
{
    char FUNC_NAME[] = "subaeroretwat";   /* function name */
    char errmsg[STR_SIZE];       /* error message */
    int iaot;               /* aerosol optical thickness (AOT) index */
    int retval;             /* function return value */
    int ib;                 /* band index */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1, ros3;      /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    float raotsaved;        /* save the raot value */
    float next;             /* ???? */
    float tgo;              /* other gaseous transmittance */
    float roatm;            /* atmospheric intrinsic reflectance */
    float ttatmg;           /* total atmospheric transmission */
    float satm;             /* spherical albedo */
    float xrorayp;          /* reflectance of the atmosphere due to molecular
                               (Rayleigh) scattering */
    double residual1, residual2;  /* residuals for storing and comparing */
    double residualm;       /* local model residual */
    int nbval;              /* number of values meeting criteria */
    bool testth;            /* surface reflectance test variable */
    double xa, xb, xc, xd, xe, xf;  /* AOT ratio values */
    double coefa, coefb;    /* AOT ratio coefficients */
    double raotmin;         /* minimum AOT ratio */
    int iaot1, iaot2;       /* AOT indices (0-based) */
    float tth[NSR_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 1.0e-03, 0.0,
                            1.0e-04, 0.0}; /* constant values for comparing
                                              against the surface reflectance */

    /* Correct band 3 and band 1 with increasing AOT (using pre till ratio is
       equal to erelc[2]) */
    iaot = *iaots;
    residual1 = 2000.0;
    residual2 = 1000.0;
    iaot2 = 0;
    iaot1 = 0;
    ros1 = 1.0;
    ros3 = 1.0;
    raot1 = 0.0001;
    raot2 = 1.0e-06;
    raot550nm = aot550nm[iaot];
    testth = false;

    /* Atmospheric correction for band 1 */
    ib = iband1;
    retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
        pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep, xtvmin,
        sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts, ttv, uoz,
        uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
        oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg, &satm, &xrorayp,
        &next, eps);
    if (retval != SUCCESS)
    {
        sprintf (errmsg, "Performing lambertian atmospheric correction "
            "type 2.");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    if (roslamb - tth[iband1] < 0.0)
        testth = true;
    ros1 = roslamb;

    /* Atmospheric correction for each band */
    nbval = 0;
    *residual = 0.0;
    for (ib = DN_BAND1; ib < DN_BAND8; ib++)
    {
        if (erelc[ib] > 0.0)
        {
            retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm,
                ib, pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin,
                xtvstep, xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi,
                tts, indts, ttv, uoz, uwv, tauray, ogtransa1, ogtransb0,
                ogtransb1, wvtransa, wvtransb, oztransa, troatm[ib], &roslamb,
                &tgo, &roatm, &ttatmg, &satm, &xrorayp, &next, eps);
            if (retval != SUCCESS)
            {
                sprintf (errmsg, "Performing lambertian atmospheric correction "
                    "type 2.");
                error_handler (true, FUNC_NAME, errmsg);
                return (ERROR);
            }

            if (roslamb - tth[ib] < 0.0)
                testth = true;
            *residual += roslamb * roslamb;
            nbval++;
        }
    }
    *residual = sqrt (*residual) / nbval;

    /* Loop until we converge on a solution */
    iaot++;
    while ((iaot < 22) && (*residual < residual1) && (!testth))
    {
        /* Reset variables for this loop */
        residual2 = residual1;
        iaot2 = iaot1;
        raot2 = raot1;
        residual1 = *residual;
        raot1 = raot550nm;
        iaot1 = iaot;
        raot550nm = aot550nm[iaot];

        /* Atmospheric correction for band 1 */
        ib = iband1;
        testth = false;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
            pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
            &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        *residual = 0.0;
        for (ib = DN_BAND1; ib < DN_BAND8; ib++)
        {
            if (erelc[ib] > 0.0)
            {
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }
    
                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                *residual += roslamb * roslamb;
                nbval++;
            }
        }
        *residual = sqrt (*residual) / nbval;

        /* Move to the next AOT index */
        iaot++;
    }

    /* If a minimum local was not reached for raot1, then just use the
       raot550nm value.  Otherwise continue to refine the raot. */
    if (iaot == 1)
    {
        *raot = raot550nm;
    }
    else
    {
        /* Refine the AOT ratio */
        *raot = raot550nm;
        raotsaved = *raot;
        xa = (raot1 * raot1) - (*raot * *raot);
        xd = (raot2 * raot2) - (*raot * *raot);
        xb = raot1 - *raot;
        xe = raot2 - *raot;
        xc = residual1 - *residual;
        xf = residual2 - *residual;
        coefa = (xc * xe - xb * xf) / (xa * xe - xb * xd);
        coefb = (xa * xf - xc * xd) / (xa * xe - xb * xd);
        raotmin = -coefb / (2.0 * coefa);

        /* Validate the min AOT ratio */
        if (raotmin < 0.01 || raotmin > 4.0)
            raotmin = *raot;

        /* Atmospheric correction for band 1 */
        raot550nm = raotmin;
        ib = iband1;
        testth = false;
        retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi, raot550nm, ib,
            pres, tpres, aot550nm, rolutt, transt, xtsstep, xtsmin, xtvstep,
            xtvmin, sphalbt, normext, tsmax, tsmin, nbfic, nbfi, tts, indts,
            ttv, uoz, uwv, tauray, ogtransa1, ogtransb0, ogtransb1, wvtransa,
            wvtransb, oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
            &satm, &xrorayp, &next, eps);
        if (retval != SUCCESS)
        {
            sprintf (errmsg, "Performing lambertian atmospheric correction "
                "type 2.");
            error_handler (true, FUNC_NAME, errmsg);
            return (ERROR);
        }

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        residualm = 0.0;
        for (ib = DN_BAND1; ib < DN_BAND8; ib++)
        {
            if (erelc[ib] > 0.0)
            {
                retval = atmcorlamb2 (xts, xtv, xmus, xmuv, xfi, cosxfi,
                    raot550nm, ib, pres, tpres, aot550nm, rolutt, transt,
                    xtsstep, xtsmin, xtvstep, xtvmin, sphalbt, normext, tsmax,
                    tsmin, nbfic, nbfi, tts, indts, ttv, uoz, uwv, tauray,
                    ogtransa1, ogtransb0, ogtransb1, wvtransa, wvtransb,
                    oztransa, troatm[ib], &roslamb, &tgo, &roatm, &ttatmg,
                    &satm, &xrorayp, &next, eps);
                if (retval != SUCCESS)
                {
                    sprintf (errmsg, "Performing lambertian atmospheric "
                        "correction type 2.");
                    error_handler (true, FUNC_NAME, errmsg);
                    return (ERROR);
                }

                if (roslamb - tth[ib] < 0.0)
                    testth = true;
                residualm += roslamb * roslamb;
                nbval++;
            }
        }

        residualm = sqrt (residualm) / nbval;
        *raot = raot550nm;

        /* Check the residuals and reset the AOT ratio */
        if (residualm > *residual)
        {
            residualm = *residual;
            *raot = raotsaved;
        }
        if (residualm > residual1)
        {
            residualm = residual1;
            *raot = raot1;
        }
        if (residualm > residual2)
        {
            residualm = residual2;
            *raot = raot2;
        }

        *residual = residualm;
        if (iaot == 1)
        {
            *iaots = 0;
/* TODO -- should this truly be an error and therefore return ERROR?? */
        }
        else
        {
            *iaots = MAX ((iaot2 - 3), 0);
        }
    }

/* TODO -- should this truly be an error and therefore return ERROR??  If so,
           should the same occur for subaeroret above? */
//    if (testth)
//        return ERROR;

    /* Successful completion */
    return (SUCCESS);
}
