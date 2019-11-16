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
MODULE:  subaeroret_new

PURPOSE:  Main driver for the atmospheric correction.  This subroutine uses
atmospheric coefficients to determine the atmospheric variables, then performs
the atmospheric corrections.


RETURN VALUE:
Type = N/A

NOTES:
******************************************************************************/
void subaeroret_new
(
    Sat_t sat,                             /* I: satellite */
    int iband1,                            /* I: band 1 index (0-based) */
    int iband3,                            /* I: band 3 index (0-based) */
    float erelc[NSR_BANDS],                /* I: band ratio variable */
    float troatm[NSR_BANDS],               /* I: toa reflectance */
    float tgo_arr[NREFL_BANDS],            /* I: per-band other gaseous
                                                 transmittance */
    float xrorayp_arr[NREFL_BANDS],        /* I: per-band reflectance of the
                                                 atmosphere due to molecular
                                                 (Rayleigh) scattering */
    int roatm_iaMax[NREFL_BANDS],          /* I: roatm_iaMax */
    float roatm_coef[NREFL_BANDS][NCOEF],  /* I: per band polynomial
                                                 coefficients for roatm */
    float ttatmg_coef[NREFL_BANDS][NCOEF], /* I: per band polynomial
                                                 coefficients for ttatmg */
    float satm_coef[NREFL_BANDS][NCOEF],   /* I: per band polynomial
                                                 coefficients for satm */
    float normext_p0a3_arr[NREFL_BANDS],   /* I: normext[iband][0][3] */
    float *raot,     /* O: AOT reflectance */
    float *residual, /* O: model residual */
    int *iaots,      /* I/O: AOT index that is passed in and out for multiple
                             calls (0-based) */
    float eps        /* I: angstroem coefficient; spectral dependency of AOT */
)
{
    int iaot;               /* aerosol optical thickness (AOT) index */
    int ib;                 /* band index */
    int start_band = 0;     /* starting band index for the loop */
    int end_band = 0;       /* ending band index for the loop */
    float raot550nm=0.0;    /* nearest input value of AOT */
    float roslamb;          /* lambertian surface reflectance */
    double ros1;            /* surface reflectance for bands */
    double raot1, raot2;    /* AOT ratios that bracket the predicted ratio */
    float raotsaved;        /* save the raot value */
    double residual1, residual2;  /* residuals for storing and comparing */
    double residualm;       /* local model residual */
    int nbval;              /* number of values meeting criteria */
    bool testth;            /* surface reflectance test variable */
    double xa, xb, xc, xd, xe, xf;  /* AOT ratio values */
    double coefa, coefb;    /* AOT ratio coefficients */
    double raotmin;         /* minimum AOT ratio */
    int iaot1, iaot2;       /* AOT indices (0-based) */
    float *tth = NULL;      /* pointer to the L8 or Sentinel tth array */
    float l8_tth[NSR_L8_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0,
                                  1.0e-04, 0.0};
                            /* constant values for comparing against the
                               L8 surface reflectance */
    float s2_tth[NSR_S2_BANDS] = {1.0e-03, 1.0e-03, 0.0, 1.0e-03, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 1.0e-04};
                            /* constant values for comparing against the
                               S2 surface reflectance (removed band 9&10) */
    float aot550nm[NAOT_VALS] = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6,
                                 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.6,
                                 3.0, 3.5, 4.0, 4.5, 5.0}; /* AOT values */

    /* Initialize variables based on the satellite type */
    if (sat == SAT_LANDSAT_8)
    {
        tth = l8_tth;
        start_band = DN_L8_BAND1;
        end_band = DN_L8_BAND7;
    }
    else if (sat == SAT_SENTINEL_2)
    {
        tth = s2_tth;
        start_band = DN_S2_BAND1;
        end_band = DN_S2_BAND12;
    }

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
    raot550nm = aot550nm[iaot];
    testth = false;

    /* Atmospheric correction for band 1 */
    ib = iband1;
    atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
        aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0], &ttatmg_coef[ib][0],
        &satm_coef[ib][0], raot550nm, ib, normext_p0a3_arr[ib], troatm[ib],
        &roslamb, eps);

    if (roslamb - tth[iband1] < 0.0)
        testth = true;
    ros1 = roslamb;

    /* Atmospheric correction for each band */
    nbval = 0;
    *residual = 0.0;
    for (ib = start_band; ib <= end_band; ib++)
    {
        /* Don't reprocess iband1 */
        if ((erelc[ib] > 0.0) && (ib != iband1))
        {
            atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
                aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);

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
    while ((iaot < NAOT_VALS) && (*residual < residual1) && (!testth))
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
        atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
            aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0], &ttatmg_coef[ib][0],
            &satm_coef[ib][0], raot550nm, ib, normext_p0a3_arr[ib], troatm[ib],
            &roslamb, eps);

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        *residual = 0.0;
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if ((erelc[ib] > 0.0) && (ib != iband1))
            {
                atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
                    aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                    &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                    normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);

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
        atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
            aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0], &ttatmg_coef[ib][0],
            &satm_coef[ib][0], raot550nm, ib, normext_p0a3_arr[ib], troatm[ib],
            &roslamb, eps);

        if (roslamb - tth[iband1] < 0.0)
            testth = true;
        ros1 = roslamb;

        /* Atmospheric correction for each band */
        nbval = 0;
        residualm = 0.0;
        for (ib = start_band; ib <= end_band; ib++)
        {
            /* Don't reprocess iband1 */
            if ((erelc[ib] > 0.0) && (ib != iband1))
            {
                atmcorlamb2_new (sat, tgo_arr[ib], xrorayp_arr[ib],
                    aot550nm[roatm_iaMax[ib]], &roatm_coef[ib][0],
                    &ttatmg_coef[ib][0], &satm_coef[ib][0], raot550nm, ib,
                    normext_p0a3_arr[ib], troatm[ib], &roslamb, eps);

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
}
