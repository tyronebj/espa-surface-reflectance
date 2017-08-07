#ifndef _AERO_INTERP_H_
#define _AERO_INTERP_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "lasrc.h"

void aerosol_interp
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will fill in the pixels for the remaining, non-center
                          pixels of the window. */
    int nlines,        /* I: number of lines in qaband, taero bands */
    int nsamps         /* I: number of samps in qaband, taero bands */
);

int aerosol_window_interp
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will interpolate/average the pixels of the windows
                          that failed the aerosol inversion (using ipflag) */
    float *teps,       /* I/O: eps (angstrom coefficient) for each pixel,
                               nlines x nsamps.  It is expected that the teps
                               values are computed for the center of the
                               aerosol windows. */
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
);

#endif
