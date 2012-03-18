#include <ruby.h>
#include <IpStdCInterface.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#define QP_MAX_ITER 1000
#define QP_EPS_GAP 1e-10
#define QP_EPS_RESIDUALS 1e-10
#define QP_ERROR_MESSAGE_BUFFER_LEN 200

#define QP_MATRIX_ENTRY(VALUE,i,j) rb_ary_entry(rb_ary_entry(VALUE,i),j)
#define QP_MATRIX_NROW RARRAY_LEN
#define QP_MATRIX_NCOL(VALUE) RARRAY_LEN(rb_ary_entry(VALUE,0))

VALUE RubyQp = Qnil;

// ERROR HANDLING

char qp_error_message[QP_ERROR_MESSAGE_BUFFER_LEN];

// Copy error state to the error message.
//
void
qp_gsl_error_handler(const char *reason, const char *file, int line, int gsl_errno) {
  snprintf(qp_error_message, QP_ERROR_MESSAGE_BUFFER_LEN, 
           "[GSL error (%d) in %s:%d] %s", gsl_errno, file, line, reason);
}

// Return a ruby exception class appropriate to the GSL error code. 
//
VALUE
qp_error_class_for(int status) {
  VALUE error_class;

  switch(status) {
    case GSL_EBADLEN:
    case GSL_EDOM:
    case GSL_EINVAL:
      error_class = rb_eArgError;
      break;
    case GSL_ERANGE:
    case GSL_EUNDRFLW:
    case GSL_EOVRFLW:
      error_class = rb_eRangeError;
      break;
    case GSL_ENOMEM:
      error_class = rb_eNoMemError;
      break;
    case GSL_EZERODIV:
      error_class = rb_eZeroDivError;
      break;
    default:
      error_class = rb_eStandardError;
  }

  return error_class;
}

// Raise a ruby exception in case in case of a nonzero status.
//
void
qp_error_handler(int status) {
  if (status == GSL_CONTINUE) {
    rb_raise(rb_eException, "Maximum iterations reached, check problem constraints for consistency.");
  } else if (status < 0) {
    rb_raise(qp_error_class_for(status), gsl_strerror(status));
  } else if (status > 0) {
    // GSL error handler was invoked
    rb_raise(qp_error_class_for(status), qp_error_message);
  }
}

// Raise a ruby Exception if ary can't be interpreted as a GSL vector.
// 
void
qp_ensure_vector(const VALUE ary) {
  int i, len;

  // A vector argument must be a ruby Array
  Check_Type(ary, T_ARRAY);

  // A vector must have at least one element
  len = RARRAY_LEN(ary);
  if (len == 0) {
    rb_raise(rb_eArgError, "Vectors must have at least one element.");
  }

  // All entries in the Array must be Floats or Fixnums
  for (i = 0; i < len; i++) {
    if (rb_obj_is_kind_of(rb_ary_entry(ary, i), rb_cNumeric) == Qfalse) {
      rb_raise(rb_eArgError, "Vector entries must be Numeric.");
    }
  }
}

// TYPES

// Encode the parameters to the minimization problem
//   min \\Ax - b|| 
// where \\.\\ is the weighted norm given by
//   \\x\\^2 = (w_1)^2(x_1)^2 + ... + (w_d)^2(x_d)^2
typedef struct {
  int nrow;
  int ncol;
  gsl_matrix *Amat;
  gsl_vector *bvec;
  gsl_matrix *Wmat;

  // coordinate-wise constraints
  double *x_L;
  double *x_U;

  // linear constraints (constraining more than one coordinate)
  int nconstraint;
  double *g_L;
  double *g_U;
  gsl_matrix *Gmat;

} qp_minimum_distance_problem;

qp_minimum_distance_problem*
new_qp_minimum_distance_problem() {
  qp_minimum_distance_problem *prob;
  prob = malloc(sizeof(qp_minimum_distance_problem));

  prob->Amat = NULL;
  prob->bvec = NULL;
  prob->Wmat = NULL;
  prob->x_L  = NULL;
  prob->x_U  = NULL;
  prob->g_L  = NULL;
  prob->g_U  = NULL;
  prob->Gmat = NULL;

  return prob;
}

void
free_qp_minimum_distance_problem(qp_minimum_distance_problem* prob) {
  if (prob->Amat != NULL)
    gsl_matrix_free(prob->Amat);
  if (prob->bvec != NULL)
    gsl_vector_free(prob->bvec);
  if (prob->Wmat != NULL)
    gsl_matrix_free(prob->Wmat);
  if (prob->x_L != NULL)
    free(prob->x_L);
  if (prob->x_U != NULL)
    free(prob->x_U);
  if (prob->g_L != NULL)
    free(prob->g_L);
  if (prob->g_U != NULL)
    free(prob->g_U);
  if (prob->Gmat != NULL)
    gsl_matrix_free(prob->Gmat);

  free(prob);
}

// TYPE CONVERSION

// Copy memory into a (newly allocated) GSL vector
gsl_vector*
qp_new_vector(int n, double *x) {
  int i;
  gsl_vector *xvec;

  xvec = gsl_vector_alloc(n);
  for (i = 0; i < n; i++)
    gsl_vector_set(xvec, i, x[i]);

  return xvec;
}

// Raise a ruby Exception if ary can't be interpreted as a GSL matrix.
//
void
qp_ensure_matrix(const VALUE ary) {
  int i, j, nrow, ncol;

  // A matrix argument must be a ruby Array
  Check_Type(ary, T_ARRAY);

  // A matrix can't be empty
  nrow = RARRAY_LEN(ary);
  if (nrow == 0) {
    rb_raise(rb_eArgError, "Matrices must have at least one row.");
  }

  // Make sure the first value in ary is an Array and grab its length
  Check_Type(rb_ary_entry(ary, 0), T_ARRAY);
  ncol = QP_MATRIX_NCOL(ary);

  // Matrix columns can't be empty
  if (ncol == 0) {
    rb_raise(rb_eArgError, "Matrix columns must have at least one element.");
  }

  // Every entry in ary should be an Array. 
  // All sub-Arrays should have the same length
  // All entries in sub-Arrays should be Floats or Fixnums
  for (i = 1; i < nrow; i++) {
    Check_Type(rb_ary_entry(ary, i), T_ARRAY);
    if (RARRAY_LEN(rb_ary_entry(ary, i)) != ncol) {
      rb_raise(rb_eArgError, "Matrix array has rows of different lengths.");
    }
    for (j = 1; j < ncol; j++) {
      if (rb_obj_is_kind_of(QP_MATRIX_ENTRY(ary, i, j), rb_cNumeric) == Qfalse) {
        rb_raise(rb_eArgError, "Matrix entries must be Numeric.");
      }
    }
  }
}

// Construct a gsl_vector from a ruby Array. Guaranteed to work if 
// qp_ensure_vector(ary) doesn't raise an exception.
//
gsl_vector*
qp_ary_to_vector(const VALUE ary) {
  gsl_vector *vec;
  int i, len;

  len = RARRAY_LEN(ary);
  vec = gsl_vector_alloc(len);

  for (i = 0; i < len; i++) {
    gsl_vector_set(vec, i, NUM2DBL(rb_ary_entry(ary, i)));
  }

  return vec;
}

// Construct a gsl_matrix from a ruby Array. Guaranteed to work if
// qp_ensure_matrix(ary) doesn't raise an exception.
//
gsl_matrix*
qp_ary_to_matrix(const VALUE ary) {
  int i, j, nrow, ncol;
  gsl_matrix *mat;

  nrow = QP_MATRIX_NROW(ary);
  ncol = QP_MATRIX_NCOL(ary);
  mat = gsl_matrix_alloc(nrow, ncol);

  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++) {
      gsl_matrix_set(mat, i, j, NUM2DBL(QP_MATRIX_ENTRY(ary, i, j)));
    }
  }

  return mat;
}

// Construct a ruby Array from a gsl_vector.
//
VALUE
qp_vector_to_ary(const gsl_vector *vec) {
  int i;
  VALUE ary = rb_ary_new();
  for (i = 0; i < vec->size; i++) {
    rb_ary_push(ary, rb_float_new(gsl_vector_get(vec, i)));
  }
  return ary;
}

// Construct a ruby Array from a gsl_matrix. Each entry in the Array corresponds to a row
// of the matrix.
//
VALUE 
qp_matrix_to_ary(const gsl_matrix *mat) {
  int i, j;
  VALUE ary = rb_ary_new();
  for (i = 0; i < mat->size1; i++) {
    VALUE row_ary = rb_ary_new();
    for (j = 0; j < mat->size2; j++) {
      rb_ary_push(row_ary, rb_float_new(gsl_matrix_get(mat, i, j)));
    }
    rb_ary_push(ary, row_ary);
  }
  return ary;
}

// CALLBACKS

// in all callback functions that follow, user_data should point to an
// qp_minimum_distance_problem struct.

Bool
qp_eval_f(Index n, Number *x, Bool new_x, Number *obj_value, UserDataPtr user_data) {
  int i;
  qp_minimum_distance_problem *prob;
  gsl_vector *xvec, *yvec;

  prob = (qp_minimum_distance_problem*)user_data;
  xvec = qp_new_vector(n, x);
  yvec = gsl_vector_alloc(prob->nrow);
  gsl_vector_memcpy(yvec, prob->bvec);                             // y <- b
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec); // y <- Ax - y
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, yvec);  // y <- Wy
  gsl_blas_ddot(yvec, yvec, obj_value);                            // *obj_value <- (y^t)y

  gsl_vector_free(xvec);
  gsl_vector_free(yvec);
  return TRUE;
}

Bool
qp_eval_grad_f(Index n, Number *x, Bool new_x, Number *grad_f, UserDataPtr user_data) {
  int j;
  qp_minimum_distance_problem *prob;
  gsl_vector *xvec, *yvec, *zvec;

  // the jth entry of x should evaluate to 2<W^2 (Ax-b), A_j>,
  // where A_j denotes the jth column of A.

  prob = (qp_minimum_distance_problem*)user_data;
  xvec = qp_new_vector(n, x);
  yvec = gsl_vector_alloc(prob->nrow);
  gsl_vector_memcpy(yvec, prob->bvec);                             // y <- b
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec); // y <- Ax - y
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, yvec);  // y <- Wy
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, yvec);  // y <- Wy

  // y is now equal to W^2(Ax - b)

  zvec = gsl_vector_alloc(prob->nrow);
  for (j = 0; j < n; j++) {
    gsl_matrix_get_col(zvec, prob->Amat, j); // z <- A_j
    gsl_blas_ddot(yvec, zvec, &grad_f[j]);    // grad_f[j] <- (y^t)z
    grad_f[j] = 2 * grad_f[j];               // grad_f[j] <- 2 * grad_f[j]
  }

  gsl_vector_free(xvec);
  gsl_vector_free(yvec);
  gsl_vector_free(zvec);
  return TRUE;
}

Bool
qp_eval_g(Index n, Number *x, Bool new_x, Index m, Number *g, UserDataPtr user_data) {
  int i;
  qp_minimum_distance_problem *prob;
  gsl_vector *xvec, *yvec;

  prob = (qp_minimum_distance_problem*)user_data;
  xvec = qp_new_vector(n, x);
  yvec = gsl_vector_alloc(n);

  for (i = 0; i < m; i++) {
    gsl_matrix_get_row(yvec, prob->Gmat, i); // y <- G_i (ith row of G)
    gsl_blas_ddot(xvec, yvec, &g[i]);         // g[i] <- (x^t)y
  }

  gsl_vector_free(xvec);
  gsl_vector_free(yvec);
  return TRUE;
}

// TODO the values computed here could be computed once and stored in the
// problem struct.
Bool
qp_eval_jac_g(Index n, Number *x, Bool new_x, Index m, 
                Index nele_jac, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data) 
{
  // Our Jacobian is just prob->Gmat
  int i, j, k;
  double d;
  qp_minimum_distance_problem *prob;
  prob = (qp_minimum_distance_problem*)user_data;

  if (values == NULL) {
    // return the structure of the Jacobian:
    // record the indices of nonzero entries in prob->Gmat
    k = 0;
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
    k = 0;
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        if ((d = gsl_matrix_get(prob->Gmat, i, j)) != 0) {
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
  int i, j, k;

  // Since our constraints are linear, this is just the Hessian of the objective
  // function, which is 2<A_i, W^2 A_j>, where A_i and A_j refer to the ith and jth
  // column of A.

  if (values == NULL) {
    // This is a symmetric matrix, fill the lower left triangle only.
    k = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        iRow[k] = i;
        jCol[k] = j;
        k++;
      }
    }
  } else {
    qp_minimum_distance_problem *prob;
    gsl_vector *xvec, *yvec;

    prob = (qp_minimum_distance_problem*)user_data;
    xvec = gsl_vector_alloc(prob->nrow);
    yvec = gsl_vector_alloc(prob->nrow);

    k = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        gsl_matrix_get_col(xvec, prob->Amat, i);                         // x <- A_i
        gsl_matrix_get_col(yvec, prob->Amat, j);                         // y <- A_j
        gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, yvec);  // y <- Wy
        gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, yvec);  // y <- Wy
        gsl_blas_ddot(xvec, yvec, &values[k]);                           // values[k] <- (x^t)y
        values[k] = 2 * values[k];
        k++;
      }
    }

    gsl_vector_free(xvec);
    gsl_vector_free(yvec);
  }

  return TRUE;
}

// RUBY METHODS

// :call-seq:
//   RubyQp::solve_dist_full(m_mat, m_vec, a_mat, b_vec, c_mat, d_vec, w_vec = nil) => hash
//
// Find x which minimizes the expression
//   ||Mx - m||
// subject to the constraints
//   Ax = b 
//   Cx >= d,
// where +M+ is the matrix specified by _m_mat_, +m+ is the vector specified by _m_vec_, 
// and so forth for +A+, +b+, +C+, and +d+. _m_mat_ (and all matrix arguments) should be 
// an Array of Arrays. Sub-Arrays should all have the same length, and their entries must 
// be Numeric. _m_vec_ (and all vector arguments) should be an Array of Numerics. 
//
// The constraint +Cx = d+ means that each entry in the vector +Cx+ is greater than or
// equal to the corresponding entry in the vector +d+.
//
// The norm +||.||+ is a weighted vector norm with weights given by a vector +w+ (an 
// optional argument specified by _w_vec_). A larger number in the ith entry of _w_vec_
// means that the ith coordinate will have a greater effect on the distance computed by
// the norm. If _w_vec_ is not given, then +||.||+ is just the standard Euclidean norm.
// Formally,
//   ||v||^2 = (w_1)^2(v_1)^2 + ... + (w_n)^2(v_n)^2
// The weights are squared so that setting a weight of w_i in the ith coordinate means
// that distances in the ith coordinate are multiplied by w_i (instead of the square root
// of w_i).
//
// The underlying algorithm requires that there be at least one equality constraint and
// at least one inequality constraint. That is, +A+ and +C+ must have at least one row
// and +b+ and +d+ must have at least one entry.
//
// Returns a ruby Hash with the following keys and values set:
//   "solution"      => minimizing solution
//   "lagrange_eq"   => Lagrange multipliers corresponding to Ax = b
//   "lagrange_ineq" => Lagrange multipliers corresponding to Cx >= d
//   "iterations"    => number of iterations to find the solution
//
VALUE
qp_solve_dist_full(int argc, VALUE *argv, VALUE self) {
  VALUE Mary, mary, Aary, bary, Cary, dary, wary;
  gsl_matrix *Mmat, *Amat, *Cmat, *Wmat, *Qmat, *Tempmat;
  gsl_vector *mvec, *bvec, *dvec, *qvec, *tempvec;
  int i, j, len, status;
  double weight;
  VALUE qp_result;

  rb_scan_args(argc, argv, "61", &Mary, &mary, &Aary, &bary, &Cary, &dary, &wary);

  qp_ensure_matrix(Mary);
  qp_ensure_matrix(Aary);
  qp_ensure_matrix(Cary);
  qp_ensure_vector(mary);
  qp_ensure_vector(bary);
  qp_ensure_vector(dary);

  // Check argument dimensions
  if (QP_MATRIX_NROW(Mary) != RARRAY_LEN(mary))
  {
    rb_raise(rb_eArgError, "m_mat.length != m_vec.length");
  }

  if (argc > 6) {
    qp_ensure_vector(wary);
    if (QP_MATRIX_NROW(Mary) != RARRAY_LEN(wary)) {
      rb_raise(rb_eArgError, "m_mat.length != w_vec.length");
    }
  }

  Mmat = qp_ary_to_matrix(Mary);
  Amat = qp_ary_to_matrix(Aary);
  Cmat = qp_ary_to_matrix(Cary);
  mvec = qp_ary_to_vector(mary);
  bvec = qp_ary_to_vector(bary);
  dvec = qp_ary_to_vector(dary);

  if (argc > 6) {
    len = RARRAY_LEN(wary);
    Wmat = gsl_matrix_calloc(len, len);
    for (i = 0; i < len; i++) {
      weight = NUM2DBL(rb_ary_entry(wary, i));
      gsl_matrix_set(Wmat, i, i, weight*weight);
    }
  } else {
    len = Mmat->size1;
    Wmat = gsl_matrix_calloc(len, len);
    for (i = 0; i < len; i++) {
      gsl_matrix_set(Wmat, i, i, 1.0);
    }
  }

  // need to transform the problem to form (1/2)(x^t)Qx + (q^t)x
  // Q = (M^t)WM
  // q = -(M^t)Wm

  // Compute Q
  Tempmat = gsl_matrix_alloc(Wmat->size1, Mmat->size2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Wmat, Mmat, 0.0, Tempmat);
  Qmat = gsl_matrix_alloc(Mmat->size2, Mmat->size2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, Mmat, Tempmat, 0.0, Qmat);

  // Compute q
  tempvec = gsl_vector_alloc(Mmat->size1);
  gsl_blas_dgemv(CblasNoTrans, 1.0, Wmat, mvec, 0.0, tempvec);
  qvec = gsl_vector_alloc(Mmat->size2);
  gsl_blas_dgemv(CblasTrans, -1.0, Mmat, tempvec, 0.0, qvec);

  status = qp_call_cqp(&qp_result, Qmat, qvec, Amat, bvec, Cmat, dvec);

  gsl_matrix_free(Mmat);
  gsl_matrix_free(Amat);
  gsl_matrix_free(Cmat);
  gsl_matrix_free(Tempmat);
  gsl_matrix_free(Qmat);
  gsl_vector_free(mvec);
  gsl_vector_free(bvec);
  gsl_vector_free(dvec);
  gsl_vector_free(tempvec);
  gsl_vector_free(qvec);

  qp_error_handler(status);

  return qp_result;
}

// :call-seq:
//   RubyQp::solve_dist(m_mat, m_vec, a_mat, b_vec, c_mat, d_vec, w_vec = nil) => array
//
// Same as RubyQp::solve_dist_full, but returns only the solution vector as a ruby Array.
//
// VALUE
// qp_solve_dist(int argc, VALUE *argv, VALUE self) {
//   VALUE qp_result = qp_solve_dist_full(argc, argv, self);
//   return rb_hash_aref(qp_result, rb_str_new2("solution"));
// }

void Init_ruby_qp(void) {
  RubyQp = rb_define_module("RubyQp");
  rb_define_module_function(RubyQp, "solve_dist_full", qp_solve_dist_full, -1);

  gsl_set_error_handler(&qp_gsl_error_handler);
}
