#ifndef POSSOLP_H
#define POSSOLP_H

#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "error.h"

void possolp
(
    int doy,      /* I: day number in the current year */
    float tu,     /* I: acquisition hour */
    float xlon,   /* I: current longitude (deg) */
    float xlat,   /* I: current longitude (deg) */
    float *asol,  /* O: solar zenith angle (deg) */
    float *phi0   /* O: solar azimuth angle (deg) */
);

#endif
