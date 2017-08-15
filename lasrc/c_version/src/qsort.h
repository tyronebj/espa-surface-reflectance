#ifndef _QSORT_H_
#define _QSORT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

void swap
(
    float *x,   /* I/O: first value to swap */
    float *y    /* I/O: second value to swap */
);

int partition
(
    float *list,    /* I/O: input floating point array to be sorted; output
                            sorted floating point array */
    int low,        /* I: starting element in array to be sorted */
    int high        /* I: ending element in array to be sorted */
);

void quicksort
(
    float *list,    /* I/O: input floating point array to be sorted; output
                            sorted floating point array */
    int low,        /* I: starting element in array to be sorted */
    int high        /* I: ending element in array to be sorted */
);

#endif
