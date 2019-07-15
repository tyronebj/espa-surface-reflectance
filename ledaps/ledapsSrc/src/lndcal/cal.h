#ifndef CAL_H
#define CAL_H

#include "lndcal.h"
#include "bool.h"
#include "lut.h"
#include "input.h"
static const int SATU_VAL[7]={255,255,255,255,255,255,255};
static const int SATU_VAL6=255;
static const int SATU_LOW_VAL6=1;

bool Cal(Param_t *param, Lut_t *lut, int iband, Input_t *input,
         unsigned char *line_in, int16 *line_in_sun_zen, int16 *line_out,
         unsigned char *line_out_qa, int iy);

bool Cal6(Lut_t *lut, Input_t *input, unsigned char *line_in, 
  int16 *line_out, unsigned char *line_out_qa, int iy);
bool Cal6_combined(Lut_t *lut, Input_t *input, unsigned char *line_b6l,
  unsigned char *line_b6h, int16 *line_out, unsigned char *line_out_qa,
  int iy);


#endif
