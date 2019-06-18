
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "crs_matrix.h"
#include "matrix.h"

struct matrix_crs_data {
    int column;
    mprec value;
};

struct matrix_crs_data_linked_list {
    struct matrix_crs_data_linked_list* next;

    struct matrix_crs_data data;
};

struct matrix_crs {
    int size;
    int nnz;

    int* rowPtr;
    int* pivotPtr;

    struct matrix_crs_data* data; // nnz elements
};

struct matrix_crs_builder {
    MatrixBuilder superInstance;

    struct matrix_crs_data_linked_list** rows;

    int nnz;
    int size;
};

MatrixCRS* crsClone(MatrixCRS* m) {
    MatrixCRS* clone = (MatrixCRS*) malloc(sizeof (MatrixCRS));

    int size = m->size;
    int nnz = clone->nnz;

    clone->data = (struct matrix_crs_data*) malloc(nnz * sizeof (struct matrix_crs_data));
    clone->rowPtr = (int*) malloc(size * sizeof (int));

    memcpy(clone->data, m->data, nnz * sizeof (struct matrix_crs_data));
    memcpy(clone->rowPtr, m->rowPtr, size * sizeof (int));

    return (MatrixCRS*) clone;
}

int crsGetSize(MatrixCRS* m) {
    return m->size;
}

void crsBuilderAdd(MatrixBuilder* instance,
        int row, int column, mprec value) {
    struct matrix_crs_builder* builder = (struct matrix_crs_builder*) instance;

    struct matrix_crs_data_linked_list* list = builder->rows[row];
    struct matrix_crs_data_linked_list* toInsert = (struct matrix_crs_data_linked_list*) malloc(sizeof (struct matrix_crs_data_linked_list));

    toInsert->data.column = column;
    toInsert->data.value = value;

    // insert the data in the column order
    if (list == NULL || list->data.column > column) {
        // it replaces the first element
        toInsert->next = list;

        builder->rows[row] = toInsert;
    } else {
        // search for the first data that has a greater column

        struct matrix_crs_data_linked_list* aux;
        do {
            aux = list;
            list = list->next;
        } while (list != NULL && list->data.column < column);

        toInsert->next = list;
        aux->next = toInsert;
    }
}

void debugCRS(MatrixCRS* inter) {
    MatrixCRS* m = (MatrixCRS*) inter;

    int size = m->size;

    int* rows = m->rowPtr;
    struct matrix_crs_data* data = m->data;

    int i, j;
    for (i = 0; i < size; i++) {
        int ptr = rows[i];
        int ptrEnd = rows[i + 1];
        printf("Row %i:", ptr);

        for (j = ptr; j < ptrEnd; j++) {
            printf("  [%d]=%.5f", data[j].column, data[j].value);
        }

        printf("\n");
    }
}

void crsGenerateB(MatrixCRS* dm, mprec* vb) {
    int size = dm->size;
    int i;
    for (i = 0; i < size; i++) {
        int j;
        int end = dm->rowPtr[i + 1];
        mprec b = 0;
        for (j = dm->rowPtr[i]; j < end; j++) {
            b += dm->data[j].value;
        }
        vb[i] = b;
    }
}

int crsBinarySearch(struct matrix_crs_data* data, int start, int end, int key) {
    int mid;
    while (start <= end) {
        mid = (start + end) / 2;
        if (key == data[mid].column)
            return mid;
        else if (key < data[mid].column) end = mid - 1;
        else start = mid + 1;
    }

    return -1;
}

struct matrix_crs_data_linked_list* crsAddToLink(struct matrix_crs_data_linked_list* last,
        int column, mprec value) {
    struct matrix_crs_data_linked_list* list = (struct matrix_crs_data_linked_list*) malloc(sizeof (struct matrix_crs_data_linked_list));

    last->next = list;

    list->data.column = column;
    list->data.value = value;

    return list;
}

mprec crsSORSum(mprec* xv, struct matrix_crs_data* row, int start, int end) {
    mprec sum = 0;
    while (start < end) {
        sum += xv[row[start].column] * row[start].value;
        start++;
    }
    return sum;
}

int crsSOR(MatrixCRS* m, int maxIterationCount,
        mprec tolerance, mprec omega,
        mprec* bv, mprec* xv) {
    int size = m->size;
    int nnz = m->nnz;

    int* rowPtr = m->rowPtr;
    struct matrix_crs_data* data = m->data;

    mprec* antxv = (mprec*) malloc(size * sizeof (mprec));
    mprec* pivotValueCache = (mprec*) malloc(size * sizeof (mprec));

    int it = 0;
    int i;

    mprec normaNum;
    mprec normaDen;

    mprec diff;
    mprec highestDifference;

    printf("Starting SOR method with %d iterations, %f as omega and %f as tolerance\n", maxIterationCount, omega, tolerance);

    for (i = 0; i < size; i++) { //making initial vector
        int pivotIndex = m->pivotPtr[i];
        mprec pivot;
        if (pivotIndex <= -1 || pivotIndex >= nnz) pivot = 0;
        else pivot = data[pivotIndex].value;

        antxv[i] = bv[i] / pivot;

        pivotValueCache[i] = pivot;
    }

    while (1) {
        for (i = 0; i < size; i++) { //computes xv
            int pivotIndex = m->pivotPtr[i];

            mprec xSum = crsSORSum(xv, data, rowPtr[i], pivotIndex) + crsSORSum(antxv, data, pivotIndex + 1, rowPtr[i + 1]);

            xv[i] = (1 - omega) * antxv[i] + omega / pivotValueCache[i] *(bv[i] - xSum);
        }
        normaDen = normaNum = 0;
        for (i = 0; i < size; i++) { //to calc the tolerance
            diff = fabsf(xv[i] - antxv[i]);
            if (diff > normaNum)
                normaNum = diff;
            if (fabsf(xv[i]) > normaDen)
                normaDen = fabsf(xv[i]);

            antxv[i] = xv[i]; // stores the previous vector xv in antxv
        }
        highestDifference = normaNum / normaDen;

        printf("Iteration %d with tolerance %.30f\n", it++, highestDifference);

        if (highestDifference <= tolerance) {
            printf("Loop ended because the tolerance was reached\n");
            break;
        } else if (it >= maxIterationCount) {
            printf("Loop ended because the number of iterations was reached\n");
            break;
        }
    }

    return 1;
}

void* crsCreateMatrix(MatrixBuilder* instance) {
    struct matrix_crs_builder* builder = (struct matrix_crs_builder*) instance;

    int size = builder->size;

    MatrixCRS* m = (MatrixCRS*) malloc(sizeof (MatrixCRS));

    int nnz = builder->nnz;

    m->size = size;

    m->data = (struct matrix_crs_data*) malloc(nnz * sizeof (struct matrix_crs_data));
    m->rowPtr = (int*) malloc((size + 1) * sizeof (int));
    m->pivotPtr = (int*) malloc(size * sizeof (int));

    int offset = 0;
    int i;

    for (i = 0; i < size; i++) {
        m->rowPtr[i] = offset;

        struct matrix_crs_data_linked_list* list = builder->rows[i];
        while (list != NULL) {
            m->data[offset++] = list->data;

            struct matrix_crs_data_linked_list* next = list->next;
            free(list);
            list = next;
        }

        int pivotIndex = crsBinarySearch(m->data, m->rowPtr[i], offset - 1, i);

        if (pivotIndex == -1 && i > m->data[offset - 1].column) {
            pivotIndex = nnz;
        }

        m->pivotPtr[i] = pivotIndex;
    }

    m->rowPtr[i] = offset;
    m->nnz = offset;

    return (MatrixCRS*) m;
}

MatrixBuilder* newMatrixCRSBuilder(int size, int nnz) {
    struct matrix_crs_builder* builder = (struct matrix_crs_builder*) malloc(sizeof (struct matrix_crs_builder));

    builder->rows = (struct matrix_crs_data_linked_list**) calloc(size, sizeof (struct matrix_crs_data_linked_list*));

    builder->nnz = nnz;
    builder->size = size;

    builder->superInstance.addValue = crsBuilderAdd;
    builder->superInstance.createMatrix = crsCreateMatrix;

    return (MatrixBuilder*) builder;
}
