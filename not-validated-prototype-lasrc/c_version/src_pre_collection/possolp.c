/*****************************************************************************
FILE: possolp.c
  
PURPOSE: Contains functions for handling solar position calculations

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include "possolp.h"

/******************************************************************************
MODULE:  pos_fft

PURPOSE:  Determines the solar azimuth and zenith angles for the current day
and lat/long.

RETURN VALUE: N/A
******************************************************************************/
void pos_fft
(
    int doy,      /* I: day number in the current year */
    float tu,     /* I: acquisition hour */
    float xlon,   /* I: current longitude (deg) */
    float xlat,   /* I: current longitude (deg) */
    float *asol,  /* O: solar zenith angle (deg) */
    float *phi0   /* O: solar azimuth angle (deg) */
)
{
    float tsm;                         /* solar time */
    float xla;                         /* latitude in radians */
    float xj;                          /* DOY float value */
    float tet;                         /* mean solar time */
    float a1, a2, a3, a4, a5;          /* time coeffs */
    float b1, b2, b3, b4, b5, b6, b7;  /* solar declination coeffs */
    float et;
    float tsv;                  /* true solar time */
    float ah;                   /* hour angle */
    float delta;                /* solar declination */
    float amuzero;
    float elev;                 /* solar elevation */
    float az;
    float caz;
    float azim;                 /* solar azimuth (radians) */

    /* mean solar time (heure decimale) */
    tsm = tu + xlon / 15.0;
    xla = xlat * DEG2RAD;
    xj = (float) doy;
    tet = 2.0 * PI * xj / 365.0;
 
    /* time equation (in mn.dec) */
    a1 = 0.000075;
    a2 = 0.001868;
    a3 = 0.032077;
    a4 = 0.014615;
    a5 = 0.040849;
    et = a1 + a2 * cos(tet) - a3 * sin(tet) - a4 * cos(2.0 * tet) -
         a5 * sin(2.0 * tet);
    et = et * 12.0 * 60.0 / PI;
 
    /* true solar time */
    tsv = tsm + et / 60.0;
    tsv -= 12.0;
 
    /* hour angle */
    ah = tsv * 15.0 * DEG2RAD;
 
    /* solar declination (in radians) */
    b1 = 0.006918;
    b2 = 0.399912;
    b3 = 0.070257;
    b4 = 0.006758;
    b5 = 0.000907;
    b6 = 0.002697;
    b7 = 0.001480;
    delta = b1 - b2 * cos(tet) + b3 * sin(tet) - b4 * cos(2.0 * tet) +
            b5 * sin(2.0 * tet) - b6 * cos(3.0 * tet) + b7 * sin(3.0 * tet);
 
    /* Determine the elevation and azimuth */
    amuzero = sin(xla) * sin(delta) + cos(xla) * cos(delta) * cos(ah);
    elev = asin(amuzero);
    az = cos(delta) * sin(ah) / cos(elev);
    if (fabs(az) - 1.000 > 0.00000)
    {
        /* Set the azimuth to +1 or -1, depending on the sign of the angle */
        if (az >= 0.0)
            az = 1.0;
        else
            az = -1.0;
    }

    caz = (-cos(xla) * sin(delta) + sin(xla) * cos(delta) * cos(ah)) /
          cos(elev);
    azim = asin(az);

    if (caz <= 0.0)
        azim = PI - azim;
    if (caz > 0.0 && az <= 0.0)
        azim = 2.0 * PI + azim;

    azim += PI;
    if (azim > PIx2)
        azim -= PIx2;

    /* Conversion in degrees, and compute the azimuth from elevation */
    elev = elev * RAD2DEG;
    *asol = 90.0 - elev;
    *phi0 = azim * RAD2DEG;

    return;
}


/******************************************************************************
MODULE:  possolp

PURPOSE:  Determines the solar position (solar azimuth and solar zenith) for
the current day and lat/long.

RETURN VALUE: N/A
******************************************************************************/
void possolp
(
    int doy,      /* I: day number in the current year */
    float tu,     /* I: acquisition hour */
    float xlon,   /* I: current longitude (deg) */
    float xlat,   /* I: current longitude (deg) */
    float *asol,  /* O: solar zenith angle (deg) */
    float *phi0   /* O: solar azimuth angle (deg) */
)
{
    /* Basically just call pos_fft to determine the solar angles */
    pos_fft (doy, tu, xlon, xlat, asol, phi0);
  
    return;
}
