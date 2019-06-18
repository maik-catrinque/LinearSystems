
#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

void debugArrayf(mprec* array, int n) {
    if (n > 0) {
        printf("[ %+.3f ", array[0]);
        int i;
        for (i = 1; i < n; i++) {
            printf("%+.3f ", array[i]);
        }
        printf("]\n");
    } else {
        printf("[ ]\n");
    }
}

void debugArrayi(int* array, int n) {
    if (n > 0) {
        printf("[ %+.3d ", array[0]);
        int i;
        for (i = 1; i < n; i++) {
            printf("%d ", array[i]);
        }
        printf("]\n");
    } else {
        printf("[ ]\n");
    }
}
