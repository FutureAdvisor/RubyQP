#ifndef QP_GSL_H
#define QP_GSL_H

#include "ruby_qp.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#define QP_GSL_VECTOR_ALLOC(vector, n)  \
    vector = gsl_vector_alloc(n);       \
    QP_CLEANUP_ON_NULL(vector);

#define QP_GSL_VECTOR_CALLOC(vector, n) \
    vector = gsl_vector_calloc(n);      \
    QP_CLEANUP_ON_NULL(vector);

#define QP_GSL_MATRIX_ALLOC(matrix, nrow, ncol) \
    matrix = gsl_matrix_alloc(nrow, ncol);      \
    QP_CLEANUP_ON_NULL(matrix);

#define QP_GSL_MATRIX_CALLOC(matrix, nrow, ncol)    \
    matrix = gsl_matrix_calloc(nrow, ncol);         \
    QP_CLEANUP_ON_NULL(matrix);

#define QP_GSL_VECTOR_FREE(vector)  if (NULL != vector) gsl_vector_free(vector);
#define QP_GSL_MATRIX_FREE(matrix)  if (NULL != matrix) gsl_matrix_free(matrix);

#endif  // #ifndef QP_GSL_H
