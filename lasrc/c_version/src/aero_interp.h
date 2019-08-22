#ifndef _AERO_INTERP_H_
#define _AERO_INTERP_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "lasrc.h"

void aerosol_interp_l8
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    int aero_window,   /* I: size of the aerosol window */
    int half_aero_window, /* I: size of half the aerosol window */
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

void aerosol_interp_s2
(
    int aero_window,   /* I: size of the aerosol window */
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the UL of the aerosol
                               windows to match taero. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the UL of the aerosol windows.  This routine
                          will fill in the pixels for the remaining pixels of
                          the window. */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
);

float find_median_aerosol_l8
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

float find_median_aerosol_s2
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the UL of the aerosol
                             windows. */
    float *taero,      /* I: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the UL of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window */
    int nlines,        /* I: number of lines in ipflag and taero band */
    int nsamps         /* I: number of samps in ipflag and taero band */
);

void aerosol_fill_median_l8
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

void aerosol_fill_median_s2
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the UL of the aerosol
                               windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the UL of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window (S2 or L8) */
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
);

#endif
