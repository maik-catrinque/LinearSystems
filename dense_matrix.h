
/* 
 * File:   dense_matrix.h
 * Author: guilh
 *
 * Created on 17 de Outubro de 2017, 17:29
 */

#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H

#include "matrix.h"


#ifdef __cplusplus
extern "C" {
#endif

    typedef struct dense_matrix DenseMatrix;

    int denseGetSize(DenseMatrix* m);

    void denseDebug(DenseMatrix* m);
    DenseMatrix* denseClone(DenseMatrix* m);
    int denseLU(DenseMatrix* dense, int* p);
    void denseGenerateB(DenseMatrix* dm, mprec* vb);
    void denseSolveLU(DenseMatrix* m, mprec* xv, mprec* yv, int* p, mprec* b);

    int denseSOR(DenseMatrix* m, int maxIterationCount,
            mprec tolerance, mprec omega,
            mprec* bv, mprec* xv);

    MatrixBuilder* newDenseBuilder(int size, int nnz);

#ifdef __cplusplus
}
#endif

#endif /* DENSE_MATRIX_H */

