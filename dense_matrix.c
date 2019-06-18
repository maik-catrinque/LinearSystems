
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "dense_matrix.h"

const mprec sModulusMask = ((-1 << 1) & ~1) >> 1;

struct dense_matrix {
    mprec* mMatrix;
    mprec** mRows;

    int size;
};

struct dense_matrix_builder {
    MatrixBuilder mSuper;

    mprec* mMatrix;
    mprec** mRows;
};

void denseSetValue(DenseMatrix* m, int row, int column, mprec value) {
    m->mRows[row][column] = value;
}

mprec denseSORSum(mprec* xv, mprec* row, int start, int end) {
    mprec sum = 0;
    while (start < end) {
        sum += xv[start] * row[start];

        start++;
    }
    return sum;
}

int denseSOR(DenseMatrix* m, int maxIterationCount,
        mprec tolerance, mprec omega,
        mprec* bv, mprec* xv) {
    int size = m->size;
    mprec** rows = m->mRows;
    mprec* antxv = (mprec*) malloc(size * sizeof (mprec));
    int it = 0;
    int i;
    mprec normaNum;
    mprec normaDen;
    mprec diff;
    mprec highestDifference;

    printf("Starting SOR method with %d iterations, %f as omega and %f as tolerance\n", maxIterationCount, omega,tolerance);

    for(i=0;i<size;i++){ //making initial vector
        antxv[i] = bv[i]/rows[i][i];
    }

    while (1) {

        for (i = 0; i < size; i++) { //computes xv
            mprec xSum = denseSORSum(xv, rows[i], 0, i) + denseSORSum(antxv, rows[i], i + 1, size);
            xv[i] = (1-omega)*antxv[i] + omega/rows[i][i] *(bv[i] - xSum);
        }
        normaDen = normaNum = 0;
        for(i=0;i < size;i++){ //to calc the tolerance
            diff = fabsf(xv[i] - antxv[i]);
            if(diff > normaNum)
                normaNum = diff;
            if(fabsf(xv[i]) > normaDen)
                normaDen = fabsf(xv[i]);

            antxv[i] = xv[i]; // stores the previous vector xv in antxv
        }
         highestDifference = normaNum/normaDen;

        printf("Iteration %d with tolerance %.30f\n", it++, highestDifference);

        if (highestDifference <= tolerance) {
            printf("Loop ended because the tolerance was reached\n");
            break;
        }
        else if(it >= maxIterationCount){
            printf("Loop ended because the number of iterations was reached\n");
            break;
        }
    }

    return 1;
}

int denseLU(DenseMatrix* dense, int* p) {
    int size = dense->size;

    mprec** rows = dense->mRows;

    int i, j, k;

    for (i = 0; i < size; i++) p[i] = i;

    for (i = 0; i < size; i++) {
        mprec pivot = rows[i][i];
        pivot = fabsf(pivot);
        int pivotRow = i;

        // search for the highest pivot in the column
        for (j = i + 1; j < size; j++) {
            mprec candidate = rows[j][i];
            candidate = fabsf(candidate);
            if (-candidate < pivot && pivot < candidate) {
                pivot = candidate;
                pivotRow = j;
            }
        }

        if (pivot == 0) {
            // the matrix is not invertible 
            return 0;
        }

        // exchange rows (if didn't found anything, then it is the same)
        mprec* row = rows[pivotRow];
        rows[pivotRow] = rows[i];
        rows[i] = row;

        // exchange matrix P
        j = p[i];
        p[i] = p[pivotRow];
        p[pivotRow] = j;

        // recover sign if lost when take the absolute value
        pivot = rows[i][i];
        // this might lose data
        //matrix_content pivotReciprocal = (matrix_content) 1.0 / pivot;

        for (j = i + 1; j < size; j++) {
            mprec* itRow = rows[j];

            // create U components
            mprec factor = itRow[i] / pivot;
            // compute line[j] -= factor * line[i]
            for (k = i + 1; k < size; k++) {
                itRow[k] -= factor * row[k];
            }

            // create L component
            itRow[i] = factor;
        }
    }

    return 1;
}

void denseFillRowArray0(mprec** rows, mprec* columnOffset, int size) {
    int row = 0;

    while (row < size) {
        rows[row] = columnOffset;

        columnOffset += size;
        row++;
    }
}

void denseMultiplyPv(mprec* vector, int* p, int n) {
    int i, index;
    for (i = 0; i < n; i++) {
        index = p[i];
        if (i < index) {
            mprec column = vector[i];
            vector[i] = vector[index];
            vector[index] = column;
        }
    }
}

void denseMultiplyPm(mprec** matrix, int* p, int n) {
    int i, index;
    for (i = 0; i < n; i++) {
        index = p[i];
        if (i < index) {
            mprec* row = matrix[i];
            matrix[i] = matrix[index];
            matrix[index] = row;
        }
    }
}

void denseFillRowArray(DenseMatrix* m) {
    denseFillRowArray0(m->mRows, m->mMatrix, m->size);
}

DenseMatrix* denseClone(DenseMatrix* m) {
    DenseMatrix* clone = (DenseMatrix*) malloc(sizeof (DenseMatrix));

    int n = m->size;

    clone->mMatrix = (mprec*) malloc(n * n * sizeof (mprec));
    clone->mRows = (mprec**) malloc(n * sizeof (mprec*));

    memcpy(clone->mMatrix, m->mMatrix, n * n * sizeof (mprec));

    denseFillRowArray(clone);

    return (DenseMatrix*) clone;
}

void denseDebugArray(mprec** rows, int size) {
    printf("Dense matrix with size = %d\n", size);
    int i, j;
    for (i = 0; i < size; i++) {
        mprec* row = rows[i];
        printf("[");
        for (j = 0; j < size; j++) {
            mprec elem = row[j];

            printf(" %+.3f", elem);
        }
        printf("]\n");
    }
    printf("\n");
}

void denseDebug(DenseMatrix* m) {
    denseDebugArray(m->mRows, m->size);
}

int denseCheckLU(DenseMatrix* ai, DenseMatrix* lu, int* p) {
    // the lower triangular matrix excluding the diagonal is the L matrix part
    // the high triangular matrix including the diagonal is the U matrix part

    // To check the matrixes, the product LU must be equals to PA

    mprec** aux = lu->mRows;

    int n = ai->size;

    mprec* lData = (mprec*) calloc(n * n, sizeof (mprec));
    mprec** lm = (mprec**) malloc(n * sizeof (mprec*));

    mprec* uData = (mprec*) calloc(n * n, sizeof (mprec));
    mprec** um = (mprec**) malloc(n * sizeof (mprec*));

    denseFillRowArray0(lm, lData, n);
    denseFillRowArray0(um, uData, n);

    int i;

    for (i = 0; i < n; i++) {
        lm[i][i] = 1;
        um[i][i] = aux[i][i];

        int j;
        for (j = i + 1; j < n; j++) {
            lm[j][i] = aux[j][i];
            um[i][j] = aux[i][j];
        }
    }

    mprec* multiplyData = (mprec*) malloc(n * n * sizeof (mprec));
    mprec** multiply = (mprec**) malloc(n * sizeof (mprec*));

    denseFillRowArray0(multiply, multiplyData, n);

    int k, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            mprec sum = 0;
            for (k = 0; k < n; k++) {
                sum += lm[i][k] * um[k][j];
            }

            multiply[i][j] = sum;
        }
    }

    mprec maxDifference = 0;

    mprec** original = ai->mRows;

    mprec** transformed = (mprec**) malloc(n * sizeof (mprec*));

    for (i = 0; i < n; i++) {
        transformed[i] = original[p[i]];
    }

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            mprec diff = multiply[i][j] - transformed[i][j];
            diff = fabsf(diff);
            if (diff > maxDifference) {
                maxDifference = diff;
            }
        }

    printf("Highest difference found: %f\n", maxDifference);

#define DEBUG_PRINT_MATRIX_STEPS

#ifdef DEBUG_PRINT_MATRIX_STEPS
    printf("Matrix A\n");
    denseDebug(ai);
    printf("Product PA\n");
    denseDebugArray(transformed, n);

    printf("Matrix L\n");
    denseDebugArray(lm, n);

    printf("Matrix U\n");
    denseDebugArray(um, n);

    printf("Matrix LU\n");
    denseDebug(lu);

    printf("Product LU:\n");
    denseDebugArray(multiply, n);
#endif

    printf("Vector P:\n");

    for (i = 0; i < n; i++) {
        printf("%d, ", p[i]);
    }
    printf("\n\n");

#ifdef DEBUG_PRINT_MATRIX_STEPS
    denseMultiplyPm(multiply, p, n);

    printf("Product PLU:\n");
    denseDebugArray(multiply, n);
#endif

    free(transformed);

    free(multiply);
    free(multiplyData);

    free(lData);
    free(lm);

    free(uData);
    free(um);

    return 1;
}

void denseSolveLU(DenseMatrix* m, mprec* xv, mprec* yv, int* p, mprec* b) {
    int n = m->size;
    mprec** matrix = m->mRows;

    int i, j;

    // solve y
    for (i = 0; i < n; i++) {
        mprec y = b[p[i]];
        mprec* row = matrix[i];
        for (j = 0; j < i; j++) {
            y -= row[j] * yv[j];
        }

        yv[i] = y;
    }

    // solve x
    for (i = n - 1; i >= 0; i--) {
        mprec x = yv[i];
        mprec* row = matrix[i];
        for (j = i + 1; j < n; j++) {
            x -= row[j] * xv[j];
        }
        x /= row[i];

        xv[i] = x;
    }
}

int denseGetSize(DenseMatrix* m) {
    return m->size;
}

void denseGenerateB(DenseMatrix* dm, mprec* vb) {
    int n = dm->size;
    mprec** matrix = dm->mRows;

    int i, j;
    for (i = 0; i < n; i++) {
        mprec* row = matrix[i];
        mprec b = row[0];
        for (j = 1; j < n; j++) {
            b += row[j];
        }

        vb[i] = b;
    }
}

void denseAddValue(MatrixBuilder* builder, int row, int column, mprec value) {
    ((struct dense_matrix_builder*) builder)->mRows[row][column] = value;
}

void* denseBuildMatrix(MatrixBuilder* aBuilder) {
    struct dense_matrix_builder* builder = (struct dense_matrix_builder*) aBuilder;
    DenseMatrix* instance = (DenseMatrix*) malloc(sizeof (DenseMatrix));

    instance->mRows = builder->mRows;
    instance->mMatrix = builder->mMatrix;

    instance->size = builder->mSuper.size;

    free(builder);

    return instance;
}

MatrixBuilder* newDenseBuilder(int size, int nnz) {
    struct dense_matrix_builder* builder = (struct dense_matrix_builder*) malloc(sizeof (struct dense_matrix_builder));

    builder->mSuper.size = size;

    builder->mSuper.addValue = denseAddValue;
    builder->mSuper.createMatrix = denseBuildMatrix;

    builder->mMatrix = (mprec*) calloc(size * size, sizeof (mprec));
    builder->mRows = (mprec**) malloc(size * sizeof (mprec*));

    denseFillRowArray0(builder->mRows, builder->mMatrix, size);

    return (MatrixBuilder*) builder;
}
