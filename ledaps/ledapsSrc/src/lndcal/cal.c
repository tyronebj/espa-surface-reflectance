#include "cal.h"
#include "const.h"
#include "error.h"
#define nint(A)(A<0?(int)(A-0.5):(int)(A+0.5))

/* Functions */
/* !Revision:
 *
 * NOTES:
 * 1. TOA radiance and reflectance equations for Landsat 7 are available in
 *    http://landsathandbook.gsfc.nasa.gov/data_prod/prog_sect11_3.html
 * 2. The TOA reflectance gain/bias values from the MTL file (stored in the
 *    XML file after converting from LPGS to ESPA) do not account for the
 *    solar angle.  Thus the gain and bias need to be applied and then we
 *    still need to account for the solar angle.
 */

bool Cal(Param_t *param, Lut_t *lut, int iband, Input_t *input,
         unsigned char *line_in, int16 *line_in_sun_zen, int16 *line_out,
         unsigned char *line_out_qa, int iy) {
  int is,val;
  float rad_gain, rad_bias;           /* TOA radiance gain/bias */
  float refl_gain = 0.0,
        refl_bias = 0.0;              /* TOA reflectance gain/bias */
  float rad;                          /* TOA radiance value */
  float ref_conv = 0.0;               /* TOA reflectance conversion value */
  float ref;                          /* TOA reflectance value */
  float fval;                         /* temporary float value */
  float sun_zen;                      /* solar zenith angle for the current
                                         pixel (radians) */
  int nsamp= input->size.s;
  int ifill= (int)lut->in_fill;

  /* Get the TOA radiance gain/bias */
  rad_gain = lut->meta.rad_gain[iband];
  rad_bias = lut->meta.rad_bias[iband];

  /* Get the TOA reflectance gain/bias if they are available, otherwise use
     the TOA reflectance equation from the Landsat handbook. */
  if (input->meta.use_toa_refl_consts) {
    refl_gain = lut->meta.refl_gain[iband];
    refl_bias = lut->meta.refl_bias[iband];

    if ( iy==0 ) {
      printf("*** band=%1d refl gain=%f refl bias=%f "
            "cos_sun_zen(scene center)=%f\n", lut->meta.iband[iband],
            refl_gain, refl_bias, lut->cos_sun_zen);
      fflush(stdout);
    }
  }
  else {
    ref_conv = (PI * lut->dsun2) / (lut->esun[iband] * lut->cos_sun_zen);
  
    if ( iy==0 ) {
      printf("*** band=%1d rad gain=%f rad bias=%f dsun2=%f\n"
             "    ref_conv=%f=(PI*%f)/(%f*%f) ***\n", lut->meta.iband[iband],
             rad_gain, rad_bias, lut->dsun2, ref_conv, lut->dsun2,
             lut->esun[iband], lut->cos_sun_zen);
      fflush(stdout);
    }
  }

  /* Loop through the samples in the line */
  for (is = 0; is < nsamp; is++) {
    val = line_in[is];
    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      continue;
    }

    /* flag saturated pixels, added by Feng (3/23/09) */
    if (val == SATU_VAL[iband]) {
      line_out[is] = lut->out_satu;
      continue;
    }

    fval= (float)val;

    /* If the TOA reflectance gain/bias values are available, then use them.
       Otherwise compute the TOA radiance then reflectance, per the Landsat
       handbook equations. */
    rad = (rad_gain * fval) + rad_bias;
    if (input->meta.use_toa_refl_consts) {
      /* use per-pixel angles - convert the degree values to radians and then
         unscale */
      sun_zen = line_in_sun_zen[is] * 0.01 * RAD;
      ref = ((refl_gain * fval) + refl_bias) / cos (sun_zen);
    }
    else {
      ref = rad * ref_conv;
    }

    /* Apply the scaling (tied to the lut->scale_factor). Valid ranges are set
       up in lut.c as well. */
    line_out[is] = (int16)(ref * TOA_SCALE + 0.5);

    /* Cap the output using the min/max values.  Then reset the toa reflectance
       value so that it's correctly reported in the stats and the min/max
       range matches that of the image data. */
    if (line_out[is] < lut->valid_range_ref[0])
      line_out[is] = lut->valid_range_ref[0];
    else if (line_out[is] > lut->valid_range_ref[1])
      line_out[is] = lut->valid_range_ref[1];
  }  /* end for is */

  return true;
}

/* Brightness Temp Calibration for TM.  Band 6 is the only band available for
   computing brightness temp. */
bool Cal6(Lut_t *lut, Input_t *input, unsigned char *line_in, int16 *line_out,
          unsigned char *line_out_qa, int iy) {
  int is, val;
  float rad_gain, rad_bias, rad, temp;
  int nsamp= input->size_th.s;
  int ifill= (int)lut->in_fill;

  /* The gain and bias are only for band 6 in this case */
  rad_gain = lut->meta.rad_gain_th[0];
  rad_bias = lut->meta.rad_bias_th[0];
  
  if ( iy==0 ) {
    printf("*** band=%1d gain=%f bias=%f ***\n", 6, rad_gain, rad_bias);
  }

  /* Loop through the pixels in this line */
  for (is = 0; is < nsamp; is++) {
    val = line_in[is];

    /* for fill pixels */
    if (val == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      continue;
    }

    /* for saturated pixels */
    if (val >= SATU_VAL6) {
      line_out[is] = lut->out_satu;
      continue;
    }

    /* compute the TOA brightness temperature in Kelvin and apply scaling (tied
       to lut->scale_factor_th). valid ranges are set up in lut.c as well. */
    rad = (rad_gain * (float)val) + rad_bias;
    temp = lut->K2[0] / log(1.0 + (lut->K1[0]/rad));
    line_out[is] = (int16)(temp * BT_SCALE + 0.5);

    /* Cap the output using the min/max values */
    if (line_out[is] < lut->valid_range_th[0])
      line_out[is] = lut->valid_range_th[0];
    else if (line_out[is] > lut->valid_range_th[1])
      line_out[is] = lut->valid_range_th[1];
  }  /* end for is */

  return true;
}

/* Brightness Temp Calibration for ETM+.  Band 6H is the primary band used for
   brightness temp.  If Band 6H is saturated, then Band 6L is used.  If it is
   also saturated then the pixel is masked as saturated. */
#define BAND_6L 0
#define BAND_6H 1
bool Cal6_combined(Lut_t *lut, Input_t *input, unsigned char *line_b6l,
  unsigned char *line_b6h, int16 *line_out, unsigned char *line_out_qa,
  int iy) {
  int is;               /* looping variable for pixels */
  int cband;            /* index for 6L or 6H to be used in the combined band */
  int val[2];           /* input thermal value for band 6L, 6H */
  int nsamp = input->size_th.s;
  int ifill = (int)lut->in_fill;
  float rad_gain[2], rad_bias[2], rad, temp;
  float k1[2], k2[2];   /* thermal consts from the LUT */

  /* Set up the calibration coefficients */
  rad_gain[BAND_6L] = lut->meta.rad_gain_th[BAND_6L];
  rad_bias[BAND_6L] = lut->meta.rad_bias_th[BAND_6L];
  rad_gain[BAND_6H] = lut->meta.rad_gain_th[BAND_6H];
  rad_bias[BAND_6H] = lut->meta.rad_bias_th[BAND_6H];
  k1[BAND_6L] = lut->K1[BAND_6L];
  k2[BAND_6L] = lut->K2[BAND_6L];
  k1[BAND_6H] = lut->K1[BAND_6H];
  k2[BAND_6H] = lut->K2[BAND_6H];
  
  if ( iy==0 ) {
    printf("*** band=6L gain=%f bias=%f ***\n", rad_gain[BAND_6L],
      rad_bias[BAND_6L]);
    printf("*** band=6H gain=%f bias=%f ***\n", rad_gain[BAND_6H],
      rad_bias[BAND_6H]);
  }

  /* Loop through the pixels in this line */
  for (is = 0; is < nsamp; is++) {
    val[BAND_6L] = line_b6l[is];
    val[BAND_6H] = line_b6h[is];

    /* For saturated pixels. If one band is fill they are both fill. */
    if (val[BAND_6H] == ifill || line_out_qa[is]==lut->qa_fill ) {
      line_out[is] = lut->out_fill;
      continue;
    }

    /* If band 6H is not saturated, it will be used for the brightness temp.
       Otherwise if band 6L is not saturated, it will be used for the
       brightness temp. If both are saturated, then the output is saturated. */
    if ((val[BAND_6H] < SATU_VAL6) && (val[BAND_6H] > SATU_LOW_VAL6)) {
      /* use band 6H for the combined band */
      cband = BAND_6H;
    }
    else {
      if ((val[BAND_6L] < SATU_VAL6) && (val[BAND_6L] > SATU_LOW_VAL6)) {
        /* use band 6L for the combined band */
        cband = BAND_6L;
      }
      else {
        /* both bands are saturated. flag it as saturated (just in case we
           missed something by using the band 6L high saturation check).
           move on to the next pixel. */
        line_out[is] = lut->out_satu;
        line_out_qa[is] = ( 0x000001 << BAND6 );
        continue;
      }
    }

    /* Compute the TOA brightness temperature for the desired band, in Kelvin,
       and apply scaling (tied to lut->scale_factor_th). Valid ranges are set
       up in lut.c as well. */
    rad = (rad_gain[cband] * (float)val[cband]) + rad_bias[cband];
    temp = lut->K2[cband] / log(1.0 + (lut->K1[cband]/rad));
    line_out[is] = (int16)(temp * BT_SCALE + 0.5);

    /* Cap the output using the min/max values */
    if (line_out[is] < lut->valid_range_th[0])
      line_out[is] = lut->valid_range_th[0];
    else if (line_out[is] > lut->valid_range_th[1])
      line_out[is] = lut->valid_range_th[1];
  }  /* end for is */

  return true;
}
