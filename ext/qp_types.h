#ifndef QP_TYPES_H
#define QP_TYPES_H

#include "ruby_qp.h"
#include "qp_gsl.h"
#include "qp_ipopt.h"

//
// TYPES
//

// Encode the parameters to the minimization problem
//   min \\Ax - b||
// where \\.\\ is the weighted norm given by
//   \\x\\^2 = (w_1)^2(x_1)^2 + ... + (w_d)^2(x_d)^2
typedef struct {
    gsl_matrix *Amat;
    gsl_vector *bvec;
    gsl_matrix *Wmat;

    // linear constraints (constraining more than one coordinate)
    gsl_matrix *Gmat;

} qp_minimum_distance_problem;

#define QP_MININUM_DISTANCE_PROBLEM_ALLOC(prob) QP_MALLOC(prob, qp_minimum_distance_problem);

#define QP_MINIMUM_DISTANCE_PROBLEM_FREE(prob)  \
    QP_GSL_MATRIX_FREE(prob->Amat);             \
    QP_GSL_VECTOR_FREE(prob->bvec);             \
    QP_GSL_MATRIX_FREE(prob->Wmat);             \
    QP_GSL_MATRIX_FREE(prob->Gmat);             \
    QP_FREE(prob);

//
// TYPE CONVERSION
//

QP_STATUS
qp_gsl_vector_set_numbers(
    gsl_vector      *vector,
    const Number    *numbers,
    const size_t    n
    );

QP_STATUS
qp_rarray_vector_len(
    long        *pn,
    const VALUE rarray
    );

QP_STATUS
qp_rarray_matrix_size(
    long        *pnrow,
    long        *pncol,
    const VALUE rarray
    );

QP_STATUS
qp_gsl_vector_from_rarray(
    gsl_vector  *vector,
    const VALUE rarray,
    const long  n
    );

QP_STATUS
qp_numbers_from_rarray(
    Number      *numbers,
    const VALUE rarray,
    const long  n
    );

QP_STATUS
qp_gsl_matrix_from_rarray(
    gsl_matrix  *matrix,
    const VALUE rarray,
    const long  nrow,
    const long  ncol
    );

QP_STATUS
qp_gsl_matrix_set_diagonal(
    gsl_matrix  *matrix,
    const VALUE rarray,
    const long n
    );

QP_STATUS
qp_rarray_from_gsl_vector(
    VALUE               *prarray,
    const gsl_vector    *vector
    );

QP_STATUS
qp_rarray_from_gsl_matrix(
    VALUE               *prarray,
    const gsl_matrix    *matrix
    );

#endif  // #ifndef QP_TYPES_H
