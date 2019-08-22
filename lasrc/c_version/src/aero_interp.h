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
    Sat_t sat,         /* I: satellite */
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    int aero_window,   /* I: size of the aerosol window (S2 or L8) */
    int half_aero_window, /* I: size of half the aerosol window (S2 or L8) */
    int16 **sband,     /* I/O: input TOA reflectance */
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
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
);

float find_median_aerosol
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the center of the aerosol
                             windows. */
    float *taero,      /* I: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the center of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window (S2 or L8) */
    int half_aero_window, /* I: size of half the aerosol window (S2 or L8) */
    int nlines,        /* I: number of lines in taero band */
    int nsamps         /* I: number of samps in taero band */
);

void aerosol_fill_median
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the center of the
                               aerosol windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the center of the aerosol windows.  This routine
                          will interpolate/average the pixels of the windows
                          that failed the aerosol inversion (using ipflag) */
    int aero_window,   /* I: size of the aerosol window (S2 or L8) */
    int half_aero_window, /* I: size of half the aerosol window (S2 or L8) */
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
);

void ipflag_expand_urban
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
);

int aero_avg_urban
(
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps. It is expected that the ipflag
                               values are computed for all pixels. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                               updated with average values if needed */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps
                               updated with average values if needed */
    int nlines,        /* I: number of lines in ipflag band */
    int nsamps         /* I: number of samps in ipflag band */
);

#endif
