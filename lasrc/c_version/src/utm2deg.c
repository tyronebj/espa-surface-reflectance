#include "lasrc.h"
#include "time.h"
#include "math.h"
#include "espa_geoloc.h"


/******************************************************************************
MODULE:  utmtodeg

PURPOSE:  Convert UTM map coordinates in WGS84 to lat/long degrees.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error converting from UTM to degrees
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
int utmtodeg
(
    Space_def_t *space_def,  /* I: space definition structure */
    int line,                /* I: line */
    int sample,              /* I: sample */
    float *lat,              /* O: latitude */
    float *lon               /* O: longitude */
)
{
    char errmsg[STR_SIZE];           /* error message */
    char FUNC_NAME[] = "utmtodeg";   /* function name */
    int zone;                        /* abs(zone) */
    double x, y;                     /* projection coords for line, sample */
    double sa = 6378137.0;           /* WGS84 semi-major axis */
    double inv_flattening = 298.257223563;  /* inverse flattening */
    double sb, e, e2;                /* intermediate variables */
    double a1, a2, j2, j4, j6;       /* intermediate variables */
    double a, b, c, bm, s, v;        /* intermediate variables */
    double alfa, beta, delt, gama;   /* intermediate variables */
    double eps, epsi, senoheps;      /* intermediate variables */
    double nab, ta0;                 /* intermediate variables */
    double e2sq;                     /* e2 squared */
    double coslatsq;                 /* (cosine of latitude) squared */

    /* Make sure this is UTM projection or it's an error */
    if (space_def->proj_num != GCTP_UTM_PROJ)
    {
        sprintf (errmsg, "Projection must be UTM");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* Determine the projection coordinates of the line, sample */
    x = space_def->ul_corner.x + sample * space_def->pixel_size[0];
    y = space_def->ul_corner.y - line * space_def->pixel_size[1];

    sb = sa - (sa / inv_flattening);
    e = pow(sa, 2) - pow(sb, 2);
    e2 = sqrt(e) / sb;
    e2sq = pow(e2, 2);
    c = pow(sa, 2) / sb;

    /* Projection coordinate adjustments for the UTM zone.  Handle south UTM
       zones. */
    x -= 500000;
    zone = space_def->zone;
    if (space_def->zone < 0)
    {  /* south UTM zone */
        y -= 10000000.0;
        zone = -space_def->zone;
    }

    /* Perform corrections for UTM projection to lat/long */
    coslatsq = pow (cos(*lat), 2);
    s = (zone * 6.0) - 183.0;
    *lat = y / (6366197.724 * 0.9996);
    v = c / pow((1.0 + (e2sq * coslatsq)), 2) * 0.9996;
    a = x / v;
    a1 = sin(2.0 * *lat);
    a2 = a1 * coslatsq;
    j2 = *lat + (a1 * 0.5);
    j4 = ((3.0 * j2) + a2) * 0.25;
    j6 = ((5.0 * j4) + (a2 * coslatsq)) / 3.0;
    alfa = 0.75 * e2sq;
    beta = (5.0 / 3.0) * pow(alfa, 2);
    gama = (35.0 / 27.0) * pow(alfa, 3);
    bm = 0.9996 * c * (*lat - alfa*j2 + beta*j4 - gama*j6);
    b = (y - bm) / v;
    epsi = ((e2sq * pow(a, 2)) * 0.5) * coslatsq;
    eps = a * (1.0 - (epsi / 3.0));
    nab = *lat + (b * (1.0 - epsi));
    senoheps = (exp(eps) - exp(-eps)) * 0.5;
    delt = atan(senoheps / (cos(nab)));
    ta0 = atan(cos(delt) * tan(nab));
    *lon = (delt * DEG) + s;
    *lat = (*lat + (1 + e2sq * coslatsq -
                  (3.0 / 2.0) * e2sq * sin(*lat) * cos(*lat) * (ta0 - *lat)) *
                  (ta0 - *lat)) * DEG;

    /* Successful completion */
    return (SUCCESS);
}

