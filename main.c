
/*
 * File:   main.c
 *
 * Created on 17 de Outubro de 2017, 14:29
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mmio.h"
#include "matrix.h"
#include "dense_matrix.h"
#include "crs_matrix.h"

#define MATRIX_STORAGE_DENSE    1
#define MATRIX_STORAGE_CRS      2

#define MATRIX_METHOD_LU  1
#define MATRIX_METHOD_SOR 2

#define MATRIX_NO_FLAG -1

// library:
// http://math.nist.gov/MatrixMarket/mmio-c.html

struct main_instance {
    char mStorageType;
    char mMethodType;

    char* mFilename;

    int sorIterations;
    mprec sorTolerance;
    mprec sorOmega;
};

/*
 * Return 1 if success. Success happens when it is valid to set the flag, which
 * is when it is not defined yet.
 */
int instanceSetFlag(char* msgType, char* flag, char value) {
    if (*flag != MATRIX_NO_FLAG) {
        printf("%s set more than once\n", msgType);
        return 1;
    }
    *flag = value;
    return 0;
}

int instanceSetMethod(struct main_instance* instance,
        char* argNameForErrMsg, char method) {
    if (instanceSetFlag("Method", &instance->mMethodType,
            method)) {
        printf("Method already set before argument %s\n", argNameForErrMsg);
        return 1;
    }

    return 0;
}

int instanceSetStorage(struct main_instance* instance,
        char* argNameForErrMsg, char storage) {
    if (instanceSetFlag("Storage", &instance->mStorageType,
            storage)) {
        printf("Storage already set before argument %s\n", argNameForErrMsg);
        return 1;
    }

    return 0;
}

int parseInt(char* src, int* destiny) {
    return sscanf(src, "%d", destiny) == 1;
}

int parseFloat(char* src, mprec* destiny) {
    return sscanf(src, "%f", destiny) == 1;
}

int parseArgumentTypes(int argc, char** argv, struct main_instance* instance) {
    // available arguments:
    // --dense: define the storage type as dense
    // --crs: define the storage type as compressed sparse row
    // --lu: use the LU factorization method
    // --sor: use successive over-relaxation method
    // --file [filename]: Specifies the filename to parse the matrix

    instance->mFilename = NULL;
    instance->mMethodType = MATRIX_NO_FLAG;
    instance->mStorageType = MATRIX_NO_FLAG;

    printf("DEBUG: Parsing arguments\n");
    int i;
    for (i = 1; i < argc; i++) {
        char* arg = argv[i];

        if (strcmp(arg, "--dense") == 0) {
            if (instanceSetStorage(instance, arg, MATRIX_STORAGE_DENSE))
                return 1;
        } else if (strcmp(arg, "--crs") == 0) {
            if (instanceSetStorage(instance, arg, MATRIX_STORAGE_CRS))
                return 1;
        } else if (strcmp(arg, "--lu") == 0) {
            if (instanceSetMethod(instance, arg, MATRIX_METHOD_LU))
                return 1;
        } else if (strcmp(arg, "--sor") == 0) {
            if (instanceSetMethod(instance, arg, MATRIX_METHOD_SOR)) {
                return 1;
            }

            if (i + 3 > argc ||
                    !parseFloat(argv[i + 1], &instance->sorOmega) ||
                    !parseFloat(argv[i + 2], &instance->sorTolerance) ||
                    !parseInt(argv[i + 3], &instance->sorIterations)) {
                printf("Use: --sor [omega] [tolerance] [iterations]\n");
                return 1;
            }

            i += 3;
        } else if (strcmp(arg, "--file") == 0) {
            if (++i >= argc) {
                printf("Missing file name: --file [filename]\n");
                return 1;
            }

            instance->mFilename = argv[i];
        } else {
            printf("WARN: Ignored argument: %s\n", arg);
        }
    }

    printf("DEBUG: Checking arguments that are not set\n");

    if (instance->mFilename == NULL) {
        printf("File not defined. Use --file [filename]\n");
        return 1;
    } else if (instance->mMethodType == MATRIX_NO_FLAG) {
        printf("Method type not defined. Use:\n");
        printf(" --lu                                       For LU factorization\n");
        printf(" --sor  [omega] [tolerance] [iterations]    For successive over-relaxation\n");
        return 1;
    } else if (instance->mStorageType == MATRIX_NO_FLAG) {
        printf("Storage type not defined. Use:\n");
        printf(" --dense    For dense storage\n");
        printf(" --crs      For compressed sparse rows\n");
        return 1;
    }

    return 0;
}

void* parseMatrix(MatrixBuilderFactory factory, FILE* f) {
    MM_typecode typecode;

    if (mm_read_banner(f, &typecode) != 0) {
        printf("ERROR: Error when parsing matrix file banner\n");
        return NULL;
    }

    if (mm_is_complex(typecode)) {
        printf("ERROR: Unsupported matrix type\n");
        return NULL;
    }

    int m, n, nnz;

    if (mm_read_mtx_crd_size(f, &m, &n, &nnz)) {
        printf("ERROR: Error while parsing matrix size\n");
        return NULL;
    } else if (m != n) {
        printf("ERROR: Only square matrix are supported\n");
        return NULL;
    }

    printf("DEBUG: Matrix with size %d and %d non zero elements\n", n, nnz);

    MatrixBuilder* builder = factory(n, nnz);

    int i, j;
    float value;

    int k;
    for (k = 0; k < nnz; k++) {
        if (fscanf(f, "%d%d%f\n", &i, &j, &value) != 3) {
            printf("ERROR: Non zero value at %d after size header\n", k);
            return NULL;
        }

        builder->addValue(builder, i - 1, j - 1, value);
    }

    return builder->createMatrix(builder);
}

int getCountOfExchanges(int* p, int n) {
    int count = 0;
    int i;
    for (i = 0; i < n; i++) {
        count += i != p[i];
    }
    return count;
}

int main(int argc, char** argv) {
    // struct to not use global variables
    struct main_instance instance;

    if (parseArgumentTypes(argc, argv, &instance)) {
        return 1;
    }

    FILE* file = fopen(instance.mFilename, "r");

    if (file == NULL) {
        printf("File \"%s\" could not be open\n", instance.mFilename);
        return 1;
    }

    int size;

    void* matrix;

    if (instance.mStorageType == MATRIX_STORAGE_DENSE) {
        DenseMatrix* denseMatrix = (DenseMatrix*) parseMatrix(newDenseBuilder, file);

        matrix = denseMatrix;

        //denseDebug(denseMatrix);
        size = denseGetSize(denseMatrix);
    } else {
        MatrixCRS* crs = (MatrixCRS*) parseMatrix(newMatrixCRSBuilder, file);

        matrix = crs;

        //debugCRS(crs);
        size = crsGetSize(crs);
    }

    if (fclose(file)) {
        printf("WARN: Error while closing file\n");
    }

    if (matrix == NULL) {
        return 1;
    }

    mprec* b = (mprec*) malloc(size * sizeof (mprec));
    mprec* x = (mprec*) malloc(size * sizeof (mprec));

    int invertible = 1;

    if (instance.mStorageType == MATRIX_STORAGE_DENSE) {
        DenseMatrix* denseMatrix = (DenseMatrix*) matrix;

        denseGenerateB(denseMatrix, b);

        if (instance.mMethodType == MATRIX_METHOD_LU) {
            // debug use
            //DenseMatrix* copy = denseClone(denseMatrix);

            int* interchangeRows = (int*) malloc(size * sizeof (int));
            invertible = denseLU(denseMatrix, interchangeRows);

            printf("Number of lines exchanged: %d\n", getCountOfExchanges(interchangeRows, size));

            mprec* y = (mprec*) malloc(size * sizeof (mprec));
            denseSolveLU(denseMatrix, x, y, interchangeRows, b);

            printf("Y:\n");
            debugArrayf(y, size);

            free(y);
            free(interchangeRows);
        } else {
            invertible = denseSOR(denseMatrix, instance.sorIterations, instance.sorTolerance, instance.sorOmega, b, x);
        }
    } else {
        // CRS format
        MatrixCRS* m = (MatrixCRS*) matrix;

        crsGenerateB(m, b);

        if (instance.mMethodType == MATRIX_METHOD_LU) {
            printf("LU factorization is not supported for CRS format\n");
            return 0;
        } else {
            invertible = crsSOR(m, instance.sorIterations, instance.sorTolerance, instance.sorOmega, b, x);
        }
    }

    if (!invertible) {
        printf("The matrix is not invertible\n");
        return 0;
    }

    printf("B:\n");
    debugArrayf(b, size);
    printf("X:\n");
    debugArrayf(x, size);

    mprec highestDifference = 0;
    int i;
    for (i = 0; i < size; i++) {
        mprec diff = fabsf(1.0 - x[i]);
        if (highestDifference < diff) {
            highestDifference = diff;
        }
    }

    printf("Highest difference in x: %.30f\n", highestDifference);

    free(b);
    free(x);

    return 0;
}
