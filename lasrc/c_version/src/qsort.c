/*****************************************************************************
FILE: qsort.c

PURPOSE: Contains functions for handling the quicksort algorithm for sorting
a floating point array.

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

LICENSE TYPE:  NASA Open Source Agreement Version 1.3

NOTES:
*****************************************************************************/
#include <qsort.h>

/******************************************************************************
MODULE:  swap

PURPOSE:  Swaps two items in the the floating point array.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
void swap
(
    float *x,   /* I/O: first value to swap */
    float *y    /* I/O: second value to swap */
)
{
    float temp;    /* temporary storage location */

    temp = *x;
    *x = *y;
    *y = temp;
}


/******************************************************************************
MODULE:  partition

PURPOSE:  This function takes last element as pivot, places the pivot element
at its correct position in the sorted array, and places all smaller (smaller
than the pivot) to left of pivot and all greater elements to right of pivot.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           Error computing the reflectance
SUCCESS         No errors encountered

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
******************************************************************************/
int partition
(
    float *list,    /* I/O: input floating point array to be sorted; output
                            sorted floating point array */
    int low,        /* I: starting element in array to be sorted */
    int high        /* I: ending element in array to be sorted */
)
{
    float pivot = list[high];    /* pivot value */
    int j;                       /* looping variable */
    int i = (low - 1);           /* index of smaller element */
 
    for (j = low; j <= high-1; j++)
    {
        /* If current element is smaller than or equal to pivot then swap it */
        if (list[j] <= pivot)
        {
            i++;    /* increment index of smaller element */
            swap (&list[i], &list[j]);
        }
    }

    /* Swap the next element of the lower list with the pivot */
    swap (&list[i+1], &list[high]);

    /* Return the index */
    return (i + 1);
}
 

/******************************************************************************
MODULE:  quicksort

PURPOSE:  Routine to sort the input array using the quicksort algorithm.

RETURN VALUE:
Type = N/A

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

NOTES:
1. Pick an element from the list, which is the pivot.
2. Reorder the list with a rule that all elements which are less than the
   pivot come before the pivot, and all elements that are higher than the
   pivot come after the pivot.  After partitioning the list, the pivot is in
   its position.
3. With the two sub-lists, apply the above steps recursively.
******************************************************************************/
void quicksort
(
    float *list,    /* I/O: input floating point array to be sorted; output
                            sorted floating point array */
    int low,        /* I: starting element in array to be sorted */
    int high        /* I: ending element in array to be sorted */
)
{
    int pi;         /* partitioning index */

    /* Sort as long as the starting element is lower than the ending element */
    if (low < high)
    {
        /* Partition the current list; upon return list[pi] will be in the
           right place */
        pi = partition (list, low, high);
 
        /* Sort the elements before the partition and after the partition */
        quicksort (list, low, pi - 1);
        quicksort (list, pi + 1, high);
    }
}
