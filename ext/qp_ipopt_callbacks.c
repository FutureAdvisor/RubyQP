#include "qp_ipopt.h"
#include "qp_types.h"

// In all callback functions that follow, user_data should point to an
// qp_minimum_distance_problem struct.

Bool
qp_eval_f(Index n, Number *x, Bool new_x,
          Number *obj_value, UserDataPtr user_data)
{
    qp_minimum_distance_problem *prob;
    gsl_vector *xvec, *yvec, *zvec;

    prob = (qp_minimum_distance_problem *) user_data;
    xvec = yvec = zvec = NULL;

    QP_GSL_VECTOR_ALLOC(xvec, n);
    QP_CLEANUP_ON_FAIL(qp_gsl_vector_set_numbers(xvec, x, n));

    QP_GSL_VECTOR_ALLOC(yvec, prob->bvec->size);
    QP_GSL_VECTOR_ALLOC(zvec, prob->bvec->size);

    gsl_vector_memcpy(yvec, prob->bvec);                                // y <- b
    gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec);    // y <- Ax - y
    gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec);     // z <- Wy
    gsl_blas_ddot(zvec, zvec, obj_value);                               // *obj_value <- (z^t)z

qp_cleanup:
    QP_GSL_VECTOR_FREE(xvec);
    QP_GSL_VECTOR_FREE(yvec);
    QP_GSL_VECTOR_FREE(zvec);

    return TRUE;
}

Bool
qp_eval_grad_f(Index n, Number *x, Bool new_x,
               Number *grad_f, UserDataPtr user_data)
{
    Index j;
    qp_minimum_distance_problem *prob;
    gsl_vector *xvec, *yvec, *zvec, *wvec;

    prob = (qp_minimum_distance_problem *) user_data;
    xvec = yvec = zvec = wvec = NULL;

    QP_GSL_VECTOR_ALLOC(xvec, n);
    QP_CLEANUP_ON_FAIL(qp_gsl_vector_set_numbers(xvec, x, n));

    QP_GSL_VECTOR_ALLOC(yvec, prob->bvec->size);
    QP_GSL_VECTOR_ALLOC(zvec, prob->bvec->size);

    // The jth entry of x should evaluate to 2<W^2 (Ax-b), A_j>,
    // where A_j denotes the jth column of A.

    gsl_vector_memcpy(yvec, prob->bvec);                                // y <- b
    gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec);    // y <- Ax - y
    gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec);     // z <- Wy
    gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, zvec, 0.0, yvec);     // y <- Wz

    // y is now equal to W^2(Ax - b)

    QP_GSL_VECTOR_ALLOC(wvec, prob->Amat->size1);
    for (j = 0; j < n; j++) {
        gsl_matrix_get_col(wvec, prob->Amat, j);    // z <- A_j
        gsl_blas_ddot(yvec, wvec, &grad_f[j]);      // grad_f[j] <- (y^t)z
        grad_f[j] = 2 * grad_f[j];                  // grad_f[j] <- 2 * grad_f[j]
    }

qp_cleanup:
    QP_GSL_VECTOR_FREE(xvec);
    QP_GSL_VECTOR_FREE(yvec);
    QP_GSL_VECTOR_FREE(zvec);
    QP_GSL_VECTOR_FREE(wvec);

    return TRUE;
}

Bool
qp_eval_g(Index n, Number *x, Bool new_x,
          Index m, Number *g, UserDataPtr user_data)
{
    Index i;
    qp_minimum_distance_problem *prob;
    gsl_vector *xvec, *yvec;

    i = 0;
    prob = (qp_minimum_distance_problem *) user_data;
    xvec = yvec = NULL;

    QP_GSL_VECTOR_ALLOC(xvec, n);
    QP_CLEANUP_ON_FAIL(qp_gsl_vector_set_numbers(xvec, x, n));

    QP_GSL_VECTOR_ALLOC(yvec, n);

    for (i = 0; i < m; i++) {
        gsl_matrix_get_row(yvec, prob->Gmat, i);    // y <- G_i (ith row of G)
        gsl_blas_ddot(xvec, yvec, &g[i]);           // g[i] <- (x^t)y
    }

qp_cleanup:
    QP_GSL_VECTOR_FREE(xvec);
    QP_GSL_VECTOR_FREE(yvec);

    return TRUE;
}

// TODO: the values computed here could be computed once and stored in the
// problem struct.
Bool
qp_eval_jac_g(Index n, Number *x, Bool new_x,
              Index m, Index nele_jac,
              Index *iRow, Index *jCol, Number *values,
              UserDataPtr user_data)
{
    // Our Jacobian is just prob->Gmat
    Index i, j, k;
    Number d;
    qp_minimum_distance_problem *prob;

    i = j = k = 0;
    d = 0;
    prob = (qp_minimum_distance_problem *) user_data;

    if (NULL == values) {
        // return the structure of the Jacobian:
        // record the indices of nonzero entries in prob->Gmat
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                if (gsl_matrix_get(prob->Gmat, i, j) != 0) {
                    iRow[k] = i;
                    jCol[k] = j;
                    k++;
                }
            }
        }
    } else {
        // return the values of the Jacobian of the constraints:
        // record the values of nonzero entries in prob->Gmat
        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                d = gsl_matrix_get(prob->Gmat, i, j);
                if (d != 0) {
                  values[k++] = d;
                }
            }
        }
    }

    return TRUE;
}

Bool
qp_eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
          Index m, Number *lambda, Bool new_lambda,
          Index nele_hess, Index *iRow, Index *jCol,
          Number *values, UserDataPtr user_data)
{
    Index i, j, k;
    qp_minimum_distance_problem *prob;
    gsl_vector *xvec, *yvec, *zvec;

    i = j = k = 0;
    prob = (qp_minimum_distance_problem *) user_data;
    xvec = yvec = zvec = NULL;

    // Since our constraints are linear, this is just the Hessian of the objective
    // function, which is 2<A_i, W^2 A_j>, where A_i and A_j refer to the ith and jth
    // column of A.

    if (NULL == values) {
        // This is a symmetric matrix, fill the lower left triangle only.
        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                iRow[k] = i;
                jCol[k] = j;
                k++;
            }
        }
    } else {
        QP_GSL_VECTOR_ALLOC(xvec, prob->Amat->size1);
        QP_GSL_VECTOR_ALLOC(yvec, prob->Amat->size1);
        QP_GSL_VECTOR_ALLOC(zvec, prob->Amat->size1);

        for (i = 0; i < n; i++) {
            for (j = 0; j <= i; j++) {
                gsl_matrix_get_col(xvec, prob->Amat, i);                        // x <- A_i
                gsl_matrix_get_col(yvec, prob->Amat, j);                        // y <- A_j
                gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec); // z <- Wy
                gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, zvec, 0.0, yvec); // y <- Wz
                gsl_blas_ddot(xvec, yvec, &values[k]);                          // values[k] <- (x^t)y
                values[k] = 2 * values[k];
                k++;
            }
        }
    }

qp_cleanup:
    QP_GSL_VECTOR_FREE(xvec);
    QP_GSL_VECTOR_FREE(yvec);
    QP_GSL_VECTOR_FREE(zvec);

    return TRUE;
}
