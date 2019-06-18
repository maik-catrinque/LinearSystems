
/* 
 * File:   matrix.h
 * Author: guilh
 *
 * Created on 17 de Outubro de 2017, 17:06
 */

#ifndef MATRIX_H
#define MATRIX_H

#ifdef __cplusplus
extern "C" {
#endif

    typedef float mprec;

    typedef struct matrix_builder MatrixBuilder;

    struct matrix_builder {
        int size;
        
        void (*addValue)(MatrixBuilder* instance,
                int row, int column, mprec value);
        
        void* (*createMatrix)(MatrixBuilder* instance);
    };
    
    typedef MatrixBuilder* (*MatrixBuilderFactory)(int size,int nnz);

    void debugArrayf(mprec* array, int n);

    void debugArrayi(int* array, int n);
    
#ifdef __cplusplus
}
#endif

#endif /* MATRIX_H */

