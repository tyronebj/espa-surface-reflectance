#include "aero_interp.h"
#include "quick_select.h"
#include "read_level1_qa.h"
#include "read_level2_qa.h"

typedef enum {
    FORWARD,
    REVERSE
} Fill_direction_t;

/******************************************************************************
MODULE:  aerosol_interp_l8

PURPOSE:  Interpolates the L8/L9 aerosol values throughout the image using the
aerosols that were calculated for each NxN window. Also cleans up the fill
pixels in the ipflag.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_interp_l8
(
    Espa_internal_meta_t *xml_metadata, /* I: XML metadata information */
    int aero_window,   /* I: size of the aerosol window */
    int half_aero_window, /* I: size of half the aerosol window */
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
    int nlines,        /* I: number of lines in qaband & taero bands */
    int nsamps         /* I: number of samps in qaband & taero bands */
)
{
    int i;                 /* looping variable for the bands */
    int line, samp;        /* looping variable for lines and samples */
    int curr_pix;          /* current pixel in 1D arrays of nlines * nsamps */
    int center_line;       /* line for the center of the aerosol window */
    int center_line1;      /* line+1 for the center of the aerosol window */
    int center_lindex;     /* center line array index */
    int center_lindex1;    /* center line 1 array index */
    int center_samp;       /* sample for the center of the aerosol window */
    int center_samp1;      /* sample+1 for the center of the aerosol window */
    int refl_indx = -99;   /* index of band 1 or first band */
    int aero_pix11;        /* pixel location for aerosol window values
                              [lcmg][scmg] */
    int aero_pix12;        /* pixel location for aerosol window values
                              [lcmg][scmg2] */
    int aero_pix21;        /* pixel location for aerosol window values
                              [lcmg2][scmg] */
    int aero_pix22;        /* pixel location for aerosol window values
                              [lcmg2][scmg2] */
    float xaero, yaero;    /* x/y location for aerosol pixel within the overall
                              larger aerosol window grid */
    float aero11;          /* aerosol value at window line, samp */
    float aero12;          /* aerosol value at window line, samp+1 */
    float aero21;          /* aerosol value at window line+1, samp */
    float aero22;          /* aerosol value at window line+1, samp+1 */
    float u, v;            /* line, sample fractional distance from current
                              pixel (weight applied to furthest line, sample) */
    float u_x_v;           /* u * v */
    int aero_window_index_step = aero_window*nsamps; /* aerosol window array
                                                        step size */
    float aero_step = 1.0/aero_window; /* fraction of window size representing
                                          1 pixel */

    /* Use band 1 band-related metadata for the reflectance information for
       Landsat (Level 1 products). If band 1 isn't available then just use the
       first band in the XML file. */
    for (i = 0; i < xml_metadata->nbands; i++)
    {
        if (!strcmp (xml_metadata->band[i].name, "b1") &&
            !strncmp (xml_metadata->band[i].product, "L1", 2))
        {
            /* this is the index we'll use for the band info */
            refl_indx = i;
        }
    }
    if (refl_indx == -99)
        refl_indx = 1;   /* default use second band in XML file */

    /* Interpolate the aerosol data for each pixel location */
    for (line = 0, curr_pix = 0; line < nlines; line++)
    {
        /* Determine the line of the representative center pixel in the
           aerosol NxN window array */
        center_line = (int) (line * aero_step) * aero_window + half_aero_window;
        center_lindex = center_line * nsamps;

        /* Determine fractional location of this line in the aerosol window.
           Negative values are at the top of the window. */
        yaero = (line - center_line) * aero_step;
        u = yaero - (int) yaero;

        /* Determine if this pixel is closest to the line below or the line
           above. If the fractional value is in the top part of the aerosol
           window, then use the line above.  Otherwise use the line below. */
        if (u < 0.0)
        {
            center_line1 = center_line - aero_window;
            center_lindex1 = center_lindex - aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 < 0)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }
        else
        {
            center_line1 = center_line + aero_window;
            center_lindex1 = center_lindex + aero_window_index_step;

            /* If the aerosol window line value is outside the bounds of the
               scene, then just use the same line in the aerosol window */
            if (center_line1 >= nlines-1)
            {
                center_line1 = center_line;
                center_lindex1 = center_lindex;
            }
        }

        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* If this pixel is fill, then don't process */
            if (level1_qa_is_fill (qaband[curr_pix]))
                continue;

            /* Determine the sample of the representative center pixel in the
               aerosol NxN window array */
            center_samp = (int) (samp * aero_step) * aero_window +
                half_aero_window;

            /* If the current line, sample are the same as the center line,
               sample, then skip to the next pixel.  We already have the
               aerosol value. */
            if (samp == center_samp && line == center_line)
                continue;

            /* Determine fractional location of this sample in the aerosol
               window.  Negative values are at the left of the window. */
            xaero = (samp - center_samp) * aero_step;
            v = xaero - (int) xaero;

            /* Determine if this pixel is closest to the sample to the left or
               the sample to the right.  If the fractional value is on the left
               side of the aerosol window, then use the sample to the left.
               Otherwise use the sample to the right. */
            if (v < 0.0)
            {
                center_samp1 = center_samp - aero_window;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 < 0)
                    center_samp1 = center_samp;
            }
            else
            {
                center_samp1 = center_samp + aero_window;
    
                /* If the aerosol window sample value is outside the bounds of
                   the scene, then just use the same sample in the aerosol
                   window */
                if (center_samp1 >= nsamps-1)
                    center_samp1 = center_samp;
            }

            /* Determine the four aerosol window pixels to be used for
               interpolating the current pixel */
            aero_pix11 = center_lindex + center_samp;
            aero_pix12 = center_lindex + center_samp1;
            aero_pix21 = center_lindex1 + center_samp;
            aero_pix22 = center_lindex1 + center_samp1;

            /* Get the aerosol values */
            aero11 = taero[aero_pix11];
            aero12 = taero[aero_pix12];
            aero21 = taero[aero_pix21];
            aero22 = taero[aero_pix22];

            /* From here make the fractional distance positive, regardless of
               where it is in the window. */
            u = fabsf(u);
            v = fabsf(v);
            u_x_v = u * v;

            /* Interpolate the aerosol */
            taero[curr_pix] = aero11  +
                              u * (aero21 - aero11) +
                              v * (aero12 - aero11) +
                              u_x_v * (aero11 - aero12 - aero21 + aero22);

            /* Set the aerosol to window interpolated. Clear anything else. */
            ipflag[curr_pix] = (1 << IPFLAG_INTERP_WINDOW);

            /* If any of the window pixels used in the interpolation were
               water pixels, then mask this pixel with water (in addition to
               the interpolation bit already set) */
            if (lasrc_qa_is_water (ipflag[aero_pix11]) ||
                lasrc_qa_is_water (ipflag[aero_pix12]) ||
                lasrc_qa_is_water (ipflag[aero_pix21]) ||
                lasrc_qa_is_water (ipflag[aero_pix22]))
                ipflag[curr_pix] |= (1 << IPFLAG_WATER);
        }  /* end for samp */
    }  /* end for line */

    /* Clean up the ipflag in the center of the NxN windows, for the fill
       pixels. If an NxN window is a mixture of fill and non-fill, the center
       of the window can be flagged as fill and some other QA based on the
       other pixels in that window. At the end, we want fill to be fill. */
    for (curr_pix = 0; curr_pix < nlines * nsamps; curr_pix++)
    {
        if (level1_qa_is_fill (qaband[curr_pix]))
            ipflag[curr_pix] = (1 << IPFLAG_FILL);
    }
}


/******************************************************************************
MODULE:  fill_with_local_average

PURPOSE:  Replaces the invalid aerosols with a local average of clear land or
water pixel values in a WxW window, for both the aerosol and eps values.
This function can process the image in a forward direction, from upper left
to lower right, or a reverse direction, from lower right to upper left.
This allows filled pixels to propagate throughout the entire image when
use_filled is true;

RETURN VALUE: None

NOTE: "fill" has two meanings here.
   The first meaning is for pixels, usually around the edge of the image,
   that do not have any valid image data; they are marked as fill in the
   Level 1 QA band. No aerosol retrieval is done for them.
   The second meaning is for pixels that have valid image data, but the aerosol
   retrieval failed. For these pixels, the aerosol value will be set ("filled")
   based on an average of nearby valid aerosol pixels.  This is the purpose
   of this function.
******************************************************************************/
static void fill_with_local_average
(
    Fill_direction_t direction, /* I: Direction to traverse the image */
    int required_clear,  /* I: Number of required clear pixels to use for the
                            local average */
    bool use_filled,     /* I: flag to use or not use previously-filled
                            pixels (as indicated by smflag) in the local
                            average */
    uint8 *ipflag,       /* I: QA flag to assist with aerosol interpolation,
                            nlines x nsamps.  */
    bool *smflag,        /* I/O: flag to indicate if the pixel was an invalid
                            aerosol and filled (true) or not filled (false) */
    float *taero,        /* I/O: aerosol values for each pixel, nlines x nsamps
                            It is expected that the aerosol values are computed
                            for the center of the aerosol windows */
    float *teps,         /* I/O: angstrom coeff for each pixel, nlines x nsamps
                            It is expected that the eps values are computed
                            for the center of the aerosol windows */
    int aero_window,     /* I: size of the aerosol window (NxN) */
    int half_aero_window,/* I: size of half the aerosol window */
    int nlines,          /* I: number of lines in taero band */
    int nsamps,          /* I: number of samps in taero band */
    int *nbpixnf,        /* O: Number of pixels left (not filled) */
    int *nbpixtot        /* O: Number of total non-fill pixels */
)
{
    int start_line, start_samp; /* starting line and sample, depending on
                                   direction the image is traversed */
    int step;             /* step to move thru the image forward or backward */
    int line, samp;       /* looping variable for lines and samples for
                             the center pixel in the NxN window */
    int iline, isamp;     /* looping variable for lines and samples for
                             the WxW window around the current center pixel,
                             in increments of NxN */
    int curr_pix;         /* current pixel in 1D array of center nlines
                             and nsamps */
    int ilpix;            /* start of current line in 1D array of expanded
                             WxW window */
    int ipix;             /* current pixel in 1D array of expanded WxW window */
    int nbclrpix;         /* number of clear aerosol pixels in this window */
    int window_offset = HALF_FIX_AERO_WINDOW - half_aero_window;
                          /* offset from the invalid pixel to the
                             first (UL) NxN aero cell in the WxW fix window */
    double sum_aero;      /* sum of the taero values in the window */
    double sum_eps;       /* sum of the teps values in the window */

    /* Set the line, sample of where we begin processing, and the step
       direction to move through the image */
    if (direction == FORWARD)
    {
        start_line = half_aero_window;
        start_samp = half_aero_window;
        step = aero_window;
    }
    else
    {
        /* When processing in reverse, find the lower right line & sample that
           are center pixels. This is where we will start processing. Set to
           the midpoint of the last NxN window whose center is in the image. */
        start_line = ((int)round((double)nlines / aero_window) - 1) *
            aero_window + half_aero_window;
        start_samp = ((int)round((double)nsamps / aero_window) - 1) *
            aero_window + half_aero_window;
        step = -aero_window;
    }

    /* Loop through the NxN center window values looking for invalid aerosol
       retrievals */
    *nbpixnf = 0;
    *nbpixtot = 0;
    for (line = start_line; line > 0 && line < nlines; line += step)
    {
        curr_pix = line * nsamps + start_samp;
        for (samp = start_samp; samp > 0 && samp < nsamps;
             samp += step, curr_pix += step)
        {
            /* Skip fill pixels */
            if (lasrc_qa_is_fill (ipflag[curr_pix]))
                continue;

            /* Increment the count of non-fill pixels */
            (*nbpixtot)++;

            /* Only process invalid, not-already-filled aerosols */
            if (lasrc_qa_is_valid_aerosol_retrieval (ipflag[curr_pix]) ||
                smflag[curr_pix])
                continue;

            /* Loop through the WxW pixels around the current center pixel,
               in increments of NxN, to look for valid aerosol retrievals
               amongst the surrounding center pixels */
            sum_aero = 0.0;
            sum_eps = 0.0;
            nbclrpix = 0;
            for (iline = line - window_offset;
                 iline <= line + window_offset;
                 iline += aero_window)
            {
                /* Make sure the window line is valid */
                if (iline < 0 || iline >= nlines)
                    continue;

                /* Determine the location of this line in our 2D array */
                ilpix = iline * nsamps;

                /* Loop through the samples of the window */
                for (isamp = samp - window_offset;
                     isamp <= samp + window_offset;
                     isamp += aero_window)
                {
                    /* Make sure the window samp is valid */
                    if (isamp < 0 || isamp >= nsamps)
                        continue;

                    /* Determine the location of this window pixel in our
                       2D array */
                    ipix = ilpix + isamp;

                    /* Check if this window pixel had valid aerosol retrieval
                       or has already been filled */
                    if (lasrc_qa_is_valid_aerosol_retrieval (ipflag[ipix]) ||
                        (use_filled && smflag[ipix]))
                    {
                        /* Add this aerosol and eps to the average */
                        nbclrpix++;
                        sum_aero += taero[ipix];
                        sum_eps += teps[ipix];
                    }
                }  /* for isamp */
            }  /* for iline */

            /* If there are enough clear pixels for computing this average
               then use the average. Otherwise increment the counter for
               number of pixels not filled within the scene. Leave the
               pixel flagged as failed aerosol thus the averaged values
               won't be used to compute other averages in this loop. */
            if (nbclrpix >= required_clear)
            {
                taero[curr_pix] = sum_aero / nbclrpix;
                teps[curr_pix] = sum_eps / nbclrpix;
                smflag[curr_pix] = true;
            }
            else
            {
                (*nbpixnf)++;
                taero[curr_pix] = DEFAULT_AERO;
                teps[curr_pix] = DEFAULT_EPS;
            }
        }  /* for samp */
    }  /* for line */
}


/******************************************************************************
MODULE:  aerosol_interp_s2

PURPOSE:  Interpolates the S2 aerosol values throughout the image using the
aerosols that were calculated for each NxN window.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_interp_s2
(
    int aero_window,   /* I: size of the aerosol window */
    uint16 *qaband,    /* I: QA band for the input image, nlines x nsamps */
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
)
{
    int line, samp;    /* looping variable for lines and samples */
    int iline, isamp;  /* looping variable for lines and samples in the
                          aerosol window */
    int curr_pix;      /* current pixel in 1D arrays of nlines * nsamps */
    int curr_win_pix;  /* current pixel in the nxn window for atm corr */
    int next_samp_pix; /* pixel location of the next sample */
    int next_line_pix; /* pixel location of the next line */
    int next_line_samp_pix; /* pixel location of next line and next sample */
    int awline;        /* line for the next aerosol window */
    int awsamp;        /* sample for the next aerosol window */
    int sq_aero_win;   /* square of the aerosol window */
    int awline_iline;  /* awline - iline */
    int iline_line;    /* iline - line */

    /* Use the UL corner of the aerosol windows to interpolate the remaining
       pixels in the window */
    sq_aero_win = aero_window * aero_window;
    for (line = 0; line < nlines; line++)
    {
        /* Determine the current pixel */
        curr_pix = line * nsamps;

        /* Determine the line for the next aerosol window */
        awline = line + aero_window;

        for (samp = 0; samp < nsamps; samp++, curr_pix++)
        {
            /* Determine the next line and next sample to be used for
               interpolating */
            next_samp_pix = line * nsamps + (samp + aero_window);
            next_line_pix = (line + aero_window) * nsamps + samp;
            next_line_samp_pix = (line + aero_window) * nsamps +
                (samp + aero_window);

            /* Determine the sample for the next aerosol window */
            awsamp = samp + aero_window;

            /* Loop through an NxN window with the current pixel being the UL
               corner of the window */
            for (iline = line; iline < awline; iline++)
            {
                /* Skip if this isn't a valid line */
                if (iline >= nlines) continue;

                awline_iline = awline - iline;
                iline_line = iline - line;
                for (isamp = samp; isamp < awsamp; isamp++)
                {
                    /* Skip if this isn't a valid sample */
                    if (isamp >= nsamps) continue;

                    curr_win_pix = iline * nsamps + isamp;
                    taero[curr_win_pix] = taero[curr_pix] *
                        (awline_iline) * (awsamp-isamp);

                    if ((awline < nlines) && (awsamp < nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[next_samp_pix] +
                            (awsamp-isamp)*(iline_line) * taero[next_line_pix] +
                            (isamp-samp)*(iline_line) *
                             taero[next_line_samp_pix];

                    else if ((awline >= nlines) && (awsamp < nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[next_samp_pix] +
                            (awsamp-isamp)*(iline_line) * taero[curr_pix] +
                            (isamp-samp)*(iline_line) * taero[next_samp_pix];

                    else if ((awline < nlines) && (awsamp >= nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[curr_pix] +
                            (awsamp-isamp)*(iline_line) * taero[next_line_pix] +
                            (isamp-samp)*(iline_line) * taero[next_line_pix];

                    else if ((awline >= nlines) && (awsamp >= nsamps))
                        taero[curr_win_pix] +=
                            (isamp-samp)*(awline_iline) * taero[curr_pix] +
                            (awsamp-isamp)*(iline_line) * taero[curr_pix] +
                            (isamp-samp)*(iline_line) * taero[curr_pix];

                    /* Compute the average */
                    taero[curr_win_pix] /= sq_aero_win;
                }
            }

            /* If not one of the representative window pixels, flag as
               interpolated */
            if ((line % aero_window != 0) || (samp % aero_window != 0))
                ipflag[curr_pix] = (1 << IPFLAG_INTERP_WINDOW);
        }  /* for samp */
    }  /* for line */

    /* Clean up the ipflag for the fill pixels */
    for (curr_pix = 0; curr_pix < nlines * nsamps; curr_pix++)
    {
        if (level1_qa_is_fill (qaband[curr_pix]))
            ipflag[curr_pix] = (1 << IPFLAG_FILL);
    }

}


/******************************************************************************
MODULE:  aerosol_fill_median_s2

PURPOSE:  Changes the S2 failed aerosol window values (i.e. water, cloud,
shadow, urban, barren, etc.) to be the median aerosol value of the clear
land pixels.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void aerosol_fill_median_s2
(
    uint8 *ipflag,     /* I/O: QA flag to assist with aerosol interpolation,
                               nlines x nsamps.  It is expected that the ipflag
                               values are computed for the UL of the aerosol
                               windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                          It is expected that the aerosol values are computed
                          for the UL of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window */
    float median_aero, /* I: median aerosol value of clear pixels */
    int nlines,        /* I: number of lines in ipflag & taero bands */
    int nsamps         /* I: number of samps in ipflag & taero bands */
)
{
    int line, samp;       /* looping variable for lines and samples */
    int curr_pix;         /* current pixel in 1D arrays of nlines * nsamps */

    /* Loop through the UL of the NxN window pixels */
    for (line = 0; line < nlines; line += aero_window)
    {
        curr_pix = line * nsamps;
        for (samp = 0; samp < nsamps;
             samp += aero_window, curr_pix += aero_window)
        {
            /* Find any pixel flagged as failed and reset the default aerosol
               value to that of the median aerosol value */
            if (btest (ipflag[curr_pix], IPFLAG_FAILED))
                taero[curr_pix] = median_aero;
        }
    }
}


/******************************************************************************
MODULE:  find_median_aerosol_s2

PURPOSE:  Finds the median Sentinel-2 aerosol value for the valid land aerosols.

RETURN VALUE:
Type = float
Value           Description
-----           -----------
zero            Error allocating memory for the aerosol array
non-zero        Median aerosol value of the array

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
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
)
{
    char errmsg[STR_SIZE];                         /* error message */
    char FUNC_NAME[] = "find_median_aerosol_s2";   /* function name */
    int line, samp;       /* looping variable for lines and samples */
    int curr_pix;         /* current pixel in 1D arrays of nlines * nsamps */
    int nbclrpix;         /* number of clear aerosol pixels in this array */
    int nwindows;         /* number of NxN windows in the image */
    float median;         /* median clear aerosol value */
    float *aero = NULL;   /* array of the clear aerosol values */

    /* Determine how many NxN windows there are in this array of data */
    nwindows = ceil ((float) nlines / aero_window) *
               ceil ((float) nsamps / aero_window);

    /* Allocate memory for the aerosols in each window */
    aero = calloc (nwindows, sizeof (float));
    if (aero == NULL)
    {
        sprintf (errmsg, "Error allocating memory for clear aerosol array");
        error_handler (true, FUNC_NAME, errmsg);
        return (0.0);
    }

    /* Loop through the NxN UL window values and write the clear aerosol
       values to the aerosol array for determining the median */
    nbclrpix = 0;
    for (line = 0; line < nlines; line += aero_window)
    {
        curr_pix = line * nsamps;
        for (samp = 0; samp < nsamps;
             samp += aero_window, curr_pix += aero_window)
        {
            /* Process clear aerosols */
            if (btest (ipflag[curr_pix], IPFLAG_CLEAR))
            {
                aero[nbclrpix] = taero[curr_pix];
                nbclrpix++;
            }  /* if pixel is clear */
        }  /* for samp */
    }  /* for line */

    /* If no clear aerosols were available, then just return a default value */
    if (nbclrpix == 0)
        median = DEFAULT_AERO;
    else
    {
        /* Get the median of the clear pixels */
        median = quick_select (aero, nbclrpix);
    }

    /* Free memory */
    free (aero);

    /* Successful completion */
    return (median);
}


/******************************************************************************
MODULE:  fix_invalid_aerosols_l8

PURPOSE:  Fixes the invalid aerosols using multiple passes through the image,
each one replacing invalid aerosols with some type of local average.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error fixing the invalid aerosols
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
int fix_invalid_aerosols_l8
(
    uint8 *ipflag,     /* I: QA flag to assist with aerosol interpolation,
                             nlines x nsamps.  It is expected that the ipflag
                             values are computed for the center of the aerosol
                             windows. */
    float *taero,      /* I/O: aerosol values for each pixel, nlines x nsamps
                             It is expected that the aerosol values are computed
                             for the center of the aerosol windows */
    float *teps,       /* I/O: angstrom coeff for each pixel, nlines x nsamps
                             It is expected that the eps values are computed
                             for the center of the aerosol windows */
    int aero_window,   /* I: size of the aerosol window (NxN) */
    int half_aero_window, /* I: size of half the aerosol window */
    int nlines,        /* I: number of lines in taero band */
    int nsamps         /* I: number of samps in taero band */
)
{
    char errmsg[STR_SIZE];                        /* error message */
    char FUNC_NAME[] = "fix_invalid_aerosols_l8"; /* function name */
    bool *smflag = NULL;  /* flag to indicate if the pixel was an invalid
                             aerosol and filled (true) or not filled (false) */
    int nbpixnf;          /* number of pixels not filled in this scene */
    int nbpixtot;         /* number of non-fill pixels in this scene */

    /* Allocate memory for filled aerosols flag and initialize to false */
    smflag = calloc (nlines * nsamps, sizeof (bool));
    if (smflag == NULL)
    {
        sprintf (errmsg, "Error allocating memory for the fill flag for "
            "invalid aerosols");
        error_handler (true, FUNC_NAME, errmsg);
        return (ERROR);
    }

    /* First pass, require at least MIN_CLEAR_PIX valid surrounding values */
    fill_with_local_average (FORWARD, MIN_CLEAR_PIX, false, ipflag, smflag,
        taero, teps, aero_window, half_aero_window, nlines, nsamps, &nbpixnf,
        &nbpixtot);
    printf ("First pass: %d pixels were not filled out of a total %d pixels\n",
        nbpixnf, nbpixtot);

    /* If the entire scene is invalid aerosols that aren't able to be filled,
       then all aerosol and eps have been set to default values and the
       ipflag is invalid.  Skip the remaining passes. */
    if (nbpixnf == nbpixtot)
    {
        printf ("Entire scene has invalid aerosols, using default taero and "
                "teps values\n");
        free (smflag);
        return (SUCCESS);
    }

    /** Second pass, require at least 1 valid surrounding value **/
    if (nbpixnf > 0)
    {
        fill_with_local_average (FORWARD, 1, true, ipflag, smflag, taero, teps,
            aero_window, half_aero_window, nlines, nsamps, &nbpixnf, &nbpixtot);
        printf ("Second pass: %d pixels were not filled out of a total %d "
            "pixels\n", nbpixnf, nbpixtot);
    }

    /** Final reverse pass for any remaining invalid retrievals (mostly the
        UL part of the image) **/
    if (nbpixnf > 0)
    {
        fill_with_local_average (REVERSE, 1, true, ipflag, smflag, taero, teps,
            aero_window, half_aero_window, nlines, nsamps, &nbpixnf, &nbpixtot);
        printf ("Final pass: %d pixels were not filled out of a total %d "
            "pixels\n", nbpixnf, nbpixtot);
    }

    /* Free the allocated memory */
    free (smflag);

    /* Successful completion */
    return (SUCCESS);
}
