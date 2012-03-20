#include <ruby.h>
#include <IpStdCInterface.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// If the following line is uncommented, then RubyQp may raise a ruby exception depending
// on the IpoptSolve return status. Otherwise, the return status is passed through in the
// result hash.
// #define QP_RAISE_ON_ERROR 1

#define QP_MATRIX_ENTRY(VALUE,i,j) rb_ary_entry(rb_ary_entry(VALUE,i),j)
#define QP_MATRIX_NROW RARRAY_LEN
#define QP_MATRIX_NCOL(VALUE) RARRAY_LEN(rb_ary_entry(VALUE,0))

VALUE RubyQp = Qnil;

// TYPES

// TODO do all these fields need to be here?
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
qp_new_vector(int n, const double *x) {
  int i;
  gsl_vector *xvec;

  xvec = gsl_vector_alloc(n);
  for (i = 0; i < n; i++)
    gsl_vector_set(xvec, i, x[i]);

  return xvec;
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
  gsl_vector *xvec, *yvec, *zvec;

  prob = (qp_minimum_distance_problem*)user_data;

  xvec = qp_new_vector(n, x);
  yvec = gsl_vector_alloc(prob->nrow);
  zvec = gsl_vector_alloc(prob->nrow);
  gsl_vector_memcpy(yvec, prob->bvec);                             // y <- b
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec); // y <- Ax - y
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec);  // z <- Wy
  gsl_blas_ddot(zvec, zvec, obj_value);                            // *obj_value <- (z^t)z

  gsl_vector_free(xvec);
  gsl_vector_free(yvec);
  gsl_vector_free(zvec);
  return TRUE;
}

Bool
qp_eval_grad_f(Index n, Number *x, Bool new_x, Number *grad_f, UserDataPtr user_data) {
  int j;
  qp_minimum_distance_problem *prob;
  gsl_vector *xvec, *yvec, *zvec, *wvec;

  // the jth entry of x should evaluate to 2<W^2 (Ax-b), A_j>,
  // where A_j denotes the jth column of A.

  prob = (qp_minimum_distance_problem*)user_data;
  xvec = qp_new_vector(n, x);
  yvec = gsl_vector_alloc(prob->nrow);
  zvec = gsl_vector_alloc(prob->nrow);
  gsl_vector_memcpy(yvec, prob->bvec);                             // y <- b
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Amat, xvec, -1.0, yvec); // y <- Ax - y
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec);  // z <- Wy
  gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, zvec, 0.0, yvec);  // y <- Wz

  // y is now equal to W^2(Ax - b)

  wvec = gsl_vector_alloc(prob->nrow);
  for (j = 0; j < n; j++) {
    gsl_matrix_get_col(wvec, prob->Amat, j);  // z <- A_j
    gsl_blas_ddot(yvec, wvec, &grad_f[j]);    // grad_f[j] <- (y^t)z
    grad_f[j] = 2 * grad_f[j];                // grad_f[j] <- 2 * grad_f[j]
  }

  gsl_vector_free(xvec);
  gsl_vector_free(yvec);
  gsl_vector_free(zvec);
  gsl_vector_free(wvec);
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
    gsl_matrix_get_row(yvec, prob->Gmat, i);  // y <- G_i (ith row of G)
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
    gsl_vector *xvec, *yvec, *zvec;

    prob = (qp_minimum_distance_problem*)user_data;
    xvec = gsl_vector_alloc(prob->nrow);
    yvec = gsl_vector_alloc(prob->nrow);
    zvec = gsl_vector_alloc(prob->nrow);

    k = 0;
    for (i = 0; i < n; i++) {
      for (j = 0; j <= i; j++) {
        gsl_matrix_get_col(xvec, prob->Amat, i);                         // x <- A_i
        gsl_matrix_get_col(yvec, prob->Amat, j);                         // y <- A_j
        gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, yvec, 0.0, zvec);  // z <- Wy
        gsl_blas_dgemv(CblasNoTrans, 1.0, prob->Wmat, zvec, 0.0, yvec);  // y <- Wz
        gsl_blas_ddot(xvec, yvec, &values[k]);                           // values[k] <- (x^t)y
        values[k] = 2 * values[k];
        k++;
      }
    }

    gsl_vector_free(xvec);
    gsl_vector_free(yvec);
    gsl_vector_free(zvec);
  }

  return TRUE;
}

VALUE
qp_hash_from_status(enum ApplicationReturnStatus status_code) {
  VALUE hash, status, error_class, message;

  error_class = Qnil;

  switch (status_code) {
    case Solve_Succeeded:
      status = rb_str_new2("Solve_Succeeded");
      message = rb_str_new2("EXIT: Optimal Solution Found.");
      break;
    case Solved_To_Acceptable_Level:
      status = rb_str_new2("Solved_To_Acceptable_Level");
      message = rb_str_new2("EXIT: Solved To Acceptable Level.");
      break;
    case Feasible_Point_Found:
      status = rb_str_new2("Feasible_Point_Found");
      message = rb_str_new2("EXIT: Feasible point for square problem found.");
      break;
    case Infeasible_Problem_Detected:
      status = rb_str_new2("Infeasible_Problem_Detected");
      message = rb_str_new2("EXIT: Feasible point for square problem found.");
      error_class = rb_const_get(RubyQp, rb_intern("InfeasibleProblemDetectedError"));
      break;
    case Search_Direction_Becomes_Too_Small:
      status = rb_str_new2("Search_Direction_Becomes_Too_Small");
      message = rb_str_new2("EXIT: Search Direction is becoming Too Small.");
      error_class = rb_const_get(RubyQp, rb_intern("SearchDirectionBecomesTooSmallError"));
      break;
    case Diverging_Iterates:
      status = rb_str_new2("Diverging_Iterates");
      message = rb_str_new2("EXIT: Iterates diverging; problem might be unbounded.");
      error_class = rb_const_get(RubyQp, rb_intern("DivergingIteratesError"));
      break;
    case User_Requested_Stop:
      status = rb_str_new2("User_Requested_Stop");
      message = rb_str_new2("EXIT: Sopping optimization at current point as requested by user.");
      break;
    case Maximum_Iterations_Exceeded:
      status = rb_str_new2("Maximum_Iterations_Exceeded");
      message = rb_str_new2("EXIT: Maximum Number of Iterations Exceeded.");
      error_class = rb_const_get(RubyQp, rb_intern("MaximumIterationsExceededError"));
      break;
    case Restoration_Failed:
      status = rb_str_new2("Restoration_Failed");
      message = rb_str_new2("EXIT: Restoration Failed!");
      error_class = rb_const_get(RubyQp, rb_intern("RestorationFailedError"));
      break;
    case Error_In_Step_Computation:
      status = rb_str_new2("Error_In_Step_Computation");
      message = rb_str_new2("EXIT: Error in step computation (regularization becomes to large?)!");
      error_class = rb_const_get(RubyQp, rb_intern("ErrorInStepComputationError"));
      break;
    case Invalid_Option:
      status = rb_str_new2("Invalid_Option");
      message = rb_str_new2("");
      error_class = rb_const_get(RubyQp, rb_intern("InvalidOptionError"));
      break;
    case Not_Enough_Degrees_Of_Freedom:
      status = rb_str_new2("Not_Enough_Degrees_Of_Freedom");
      message = rb_str_new2("EXIT: Problem has too few degrees of freedom.");
      error_class = rb_const_get(RubyQp, rb_intern("NotEnoughDegreesOfFreedomError"));
      break;
    case Invalid_Problem_Definition:
      status = rb_str_new2("Invalid_Problem_Definition");
      message = rb_str_new2("");
      error_class = rb_const_get(RubyQp, rb_intern("InvalidProblemDefinitionError"));
      break;
    case Unrecoverable_Exception:
      status = rb_str_new2("Unrecoverable_Exception");
      message = rb_str_new2("");
      error_class = rb_const_get(RubyQp, rb_intern("UnrecoverableExceptionError"));
      break;
    case NonIpopt_Exception_Thrown:
      status = rb_str_new2("NonIpopt_Exception_Thrown");
      message = rb_str_new2("Unknown Exception caught in Ipopt");
      error_class = rb_const_get(RubyQp, rb_intern("NonIpoptExceptionThrownError"));
      break;
    case Insufficient_Memory:
      status = rb_str_new2("Insufficient_Memory");
      message = rb_str_new2("EXIT: Not enough memory.");
      error_class = rb_eNoMemError;
      break;
    case Internal_Error:
      status = rb_str_new2("Internal_Error");
      message = rb_str_new2("EXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.");
      error_class = rb_const_get(RubyQp, rb_intern("InternalError"));
      break;
  }

  hash = rb_hash_new();
  rb_hash_aset(hash, rb_str_new2("status"), status);
  rb_hash_aset(hash, rb_str_new2("message"), message);

  if (error_class != Qnil) {
    rb_hash_aset(hash, rb_str_new2("error"), error_class);
  }

  return hash;
}

// Raise an exception corresponding to the Ipopt status code
void
qp_handle_error(VALUE status_hash) {
  VALUE error_class = rb_hash_aref(status_hash, rb_str_new2("error"));

  if (error_class != Qnil) {
    rb_raise(error_class, StrinvValueCStr(rb_hash_aref(status_hash, rb_str_new2("message"))));
  }
}

// RUBY METHODS

// :call-seq:
//   RubyQp::solve_dist_full(a_mat, b_vec, x_lower, x_upper, g_mat, g_lower, g_upper, x_init, w_vec = nil) => hash
//
//   ||Ax - b||
//   x_lower <= x <= x_upper
//   g_lower <= Gx <= g_upper
//
// where +A+ is the matrix specified by _a_mat_, +b+ is the vector specified by _b_vec_, 
// and so forth. _a_mat_ (and all matrix arguments) should be 
// an Array of Arrays. Sub-Arrays should all have the same length, and their entries must 
// be Numeric. _m_vec_ (and all vector arguments) should be an Array of Numerics. 
//
// The inequality constraints apply coordinate-wise. To add an equality constraint, set the
// corresponding values in x_lower and x_upper (or in g_lower and g_upper) to the same value.
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
// _xinit_ is the initial point for the optimization algorithm.
//
// Returns a ruby Hash with the following keys and values set:
//   "solution"      => minimizing solution
//   "status"        => Ipopt status code and message
//
VALUE
qp_solve_dist_full(int argc, VALUE *argv, VALUE self) {
  int i, j, nele_jac, nele_hess;
  enum ApplicationReturnStatus status;
  double *x;
  VALUE a_mat, b_vec, x_lower, x_upper, g_mat, g_lower, g_upper, x_init, w_vec;
  gsl_matrix *matrix;
  gsl_vector *vector;
  qp_minimum_distance_problem *prob;
  IpoptProblem nlp;
  VALUE qp_solution, qp_status, qp_result;

  rb_scan_args(argc, argv, "81", &a_mat, &b_vec, &x_lower, &x_upper, &g_mat, &g_lower, &g_upper, &x_init, &w_vec);

  qp_ensure_matrix(a_mat);
  qp_ensure_matrix(g_mat);
  qp_ensure_vector(b_vec);
  qp_ensure_vector(x_lower);
  qp_ensure_vector(x_upper);
  qp_ensure_vector(g_lower);
  qp_ensure_vector(g_upper);
  qp_ensure_vector(x_init);

  // TODO add more dimension checks; no longer handled by library

  // Check argument dimensions
  if (QP_MATRIX_NROW(a_mat) != RARRAY_LEN(b_vec))
    rb_raise(rb_eArgError, "a_mat.length != b_vec.length");

  if (QP_MATRIX_NCOL(a_mat) != RARRAY_LEN(x_lower))
    rb_raise(rb_eArgError, "a_mat[0].length != x_lower.length");

  if (QP_MATRIX_NCOL(a_mat) != RARRAY_LEN(x_upper))
    rb_raise(rb_eArgError, "a_mat[0].length != x_upper.length");

  if (QP_MATRIX_NCOL(a_mat) != QP_MATRIX_NCOL(g_mat))
    rb_raise(rb_eArgError, "a_mat[0].length != g_mat[0].length");

  if (QP_MATRIX_NROW(g_mat) != RARRAY_LEN(g_lower))
    rb_raise(rb_eArgError, "g_mat.length != g_lower.length");

  if (QP_MATRIX_NROW(g_mat) != RARRAY_LEN(g_upper))
    rb_raise(rb_eArgError, "g_mat.length != g_upper.length");

  if (QP_MATRIX_NCOL(a_mat) != RARRAY_LEN(x_init))
    rb_raise(rb_eArgError, "a_mat[0].length != x_init.length");

  if (argc > 8) {
    qp_ensure_vector(w_vec);
    if (QP_MATRIX_NROW(a_mat) != RARRAY_LEN(w_vec)) {
      rb_raise(rb_eArgError, "a_mat.length != w_vec.length");
    }
  }

  // Set up the problem parameters

  prob = new_qp_minimum_distance_problem();
  prob->nrow = QP_MATRIX_NROW(a_mat);
  prob->ncol = QP_MATRIX_NCOL(a_mat);
  prob->Amat = qp_ary_to_matrix(a_mat);
  prob->bvec = qp_ary_to_vector(b_vec);

  prob->x_L = malloc(prob->ncol * sizeof(double));
  for (i = 0; i < prob->ncol; i++)
    prob->x_L[i] = NUM2DBL(rb_ary_entry(x_lower, i));

  prob->x_U = malloc(prob->ncol * sizeof(double));
  for (i = 0; i < prob->ncol; i++)
    prob->x_U[i] = NUM2DBL(rb_ary_entry(x_upper, i));
                                  
  prob->nconstraint = QP_MATRIX_NROW(g_mat);
  prob->Gmat = qp_ary_to_matrix(g_mat);

  prob->g_L = malloc(prob->nconstraint * sizeof(double));
  for (i = 0; i < prob->nconstraint; i++)
    prob->g_L[i] = NUM2DBL(rb_ary_entry(g_lower, i));

  prob->g_U = malloc(prob->nconstraint * sizeof(double));
  for (i = 0; i < prob->nconstraint; i++)
    prob->g_U[i] = NUM2DBL(rb_ary_entry(g_upper, i));

  x = malloc(prob->ncol * sizeof(double));
  for (i = 0; i < prob->ncol; i++)
    x[i] = NUM2DBL(rb_ary_entry(x_init, i));

  if (argc > 8) {
    matrix = gsl_matrix_calloc(prob->nrow, prob->nrow);
    for (i = 0; i < prob->nrow; i++)
      gsl_matrix_set(matrix, i, i, NUM2DBL(rb_ary_entry(w_vec, i)));
    prob->Wmat = matrix;
  } else {
    matrix = gsl_matrix_calloc(prob->nrow, prob->nrow);
    for (i = 0; i < prob->nrow; i++)
      gsl_matrix_set(matrix, i, i, 1.0);
    prob->Wmat = matrix;
  }

  // count the number of nonzero elements in the constraint Jacobian
  // (i.e., the number of nonzero elements in Gmat)
  nele_jac = 0;
  for (i = 0; i < prob->nconstraint; i++)
    for (j = 0; j < prob->ncol; j++)
      if (gsl_matrix_get(prob->Gmat, i, j) != 0)
        nele_jac++;

  // our Hessian is dense, so nele_hess is just the number of entries below the diagonal
  // (inclusive) in the Hessian.
  nele_hess = 0;
  for (i = 1; i <= prob->ncol; i++)
    nele_hess += i;

  nlp = CreateIpoptProblem(prob->ncol,
                           prob->x_L,
                           prob->x_U,
                           prob->nconstraint,
                           prob->g_L,
                           prob->g_U,
                           nele_jac,
                           nele_hess,
                           0,
                           &qp_eval_f,
                           &qp_eval_g,
                           &qp_eval_grad_f,
                           &qp_eval_jac_g,
                           &qp_eval_h);

  // TODO these options should be set globally in a ruby-accessible manner
  AddIpoptIntOption(nlp, "print_level", 0);
  AddIpoptIntOption(nlp, "max_iter", 200);
  AddIpoptStrOption(nlp, "mehrotra_algorithm", "yes");

  // TODO capture more of the output
  status = IpoptSolve(nlp, x, NULL, NULL, NULL, NULL, NULL, prob);

  qp_solution = rb_ary_new();
  for (i = 0; i < prob->ncol; i++)
    rb_ary_push(qp_solution, rb_float_new(x[i]));

  qp_status = qp_hash_from_status(status);
  qp_result = rb_hash_new();
  rb_hash_aset(qp_result, rb_str_new2("solution"), qp_solution);
  rb_hash_aset(qp_result, rb_str_new2("status"), qp_status);

  FreeIpoptProblem(nlp);
  free(x);
  free_qp_minimum_distance_problem(prob);

  // memory referenced by matrix and vector will be freed by free_qp_minimum_distance_problem

#ifdef QP_RAISE_ON_ERROR
  qp_handle_error(qp_status);
#endif

  return qp_result;
}

// :call-seq:
//   RubyQp::solve_dist(m_mat, m_vec, a_mat, b_vec, c_mat, d_vec, w_vec = nil) => array
//
// Same as RubyQp::solve_dist_full, but returns only the solution vector as a ruby Array.
//
VALUE
qp_solve_dist(int argc, VALUE *argv, VALUE self) {
  VALUE qp_result = qp_solve_dist_full(argc, argv, self);
  return rb_hash_aref(qp_result, rb_str_new2("solution"));
}

void 
Init_ruby_qp(void) {
  RubyQp = rb_define_module("RubyQp");
  rb_require("ruby_qp/errors");
  rb_define_module_function(RubyQp, "solve_dist_full", qp_solve_dist_full, -1);
  rb_define_module_function(RubyQp, "solve_dist", qp_solve_dist, -1);
}
