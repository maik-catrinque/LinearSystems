
/* 
 * File:   crs_matrix.h
 * Author: guilh
 *
 * Created on 19 de Outubro de 2017, 16:08
 */

#ifndef CRS_MATRIX_H
#define CRS_MATRIX_H

#include "matrix.h"


#ifdef __cplusplus
extern "C" {
#endif

    typedef struct matrix_crs MatrixCRS;
    typedef struct matrix_crs_builder MatrixCRSBuilder;

    MatrixBuilder* newMatrixCRSBuilder(int size, int nnz);

    int crsGetSize(MatrixCRS* m);

    int crsSOR(MatrixCRS* m, int maxIterationCount,
            mprec tolerance, mprec omega,
            mprec* bv, mprec* xv);
    
    void crsGenerateB(MatrixCRS* dm, mprec* vb);
    
    void debugCRS(MatrixCRS* inter);

#ifdef __cplusplus
}
#endif

#endif /* CRS_MATRIX_H */

