#include <ruby.h>
#include <IpStdCInterface.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

// If the following line is uncommented, then RubyQp may raise a ruby exception depending
// on the IpoptSolve return status. Otherwise, the return status is passed through in the
// result hash.
// #define QP_RAISE_ON_ERROR 1

#define QP_ERROR_BUFFER_LENGTH 256
char qp_error_buffer[QP_ERROR_BUFFER_LENGTH+1];

#define QP_ASSERT(TEST,ERROR_MESSAGE) if (!(TEST)) { \
                                        strncpy(qp_error_buffer, ERROR_MESSAGE, QP_ERROR_BUFFER_LENGTH); \
                                        return -1; \
                                      }

#define QP_MATRIX_ENTRY(VALUE,i,j) rb_ary_entry(rb_ary_entry(VALUE,i),j)
#define QP_MATRIX_NROW RARRAY_LEN
#define QP_MATRIX_NCOL(VALUE) RARRAY_LEN(rb_ary_entry(VALUE,0))

VALUE RubyQp = Qnil;

// TYPES

// TODO do all these fields need to be here? Callbacs don't need x_L, x_U, g_L, g_U
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

// Allocate a new gsl_vector and initialize it with the values in array.
// Returns 0 upon success, -1 otherwise. In either case, the caller needs to ensure that
// the memory allocated for vector is freed.
//
int
qp_vector_from_array(gsl_vector **pvector, const VALUE array) {
  int i, len;
  VALUE entry;

  QP_ASSERT(rb_obj_is_kind_of(array, rb_cArray) == Qtrue, "Vectors must be Arrays.");

  len = RARRAY_LEN(array);
  QP_ASSERT(len > 0, "Vectors must have at least one element.");

  *pvector = gsl_vector_alloc(len);
  for (i = 0; i < len; i++) {
    entry = rb_ary_entry(array, i);
    QP_ASSERT(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entries must be Numeric.");
    gsl_vector_set(*pvector, i, NUM2DBL(entry));
  }

  return 0;
}

// Allocate memory to a double* and initialize it with the values in array.
// Returns number of bytes written upon success, -1 otherwise. In either case, the caller 
// needs to ensure that the memory allocated for vector is freed.
//
int
qp_ptr_from_array(double **ptr, const VALUE array) {
  int i, len;
  VALUE entry;

  QP_ASSERT(rb_obj_is_kind_of(array, rb_cArray) == Qtrue, "Vectors must be Arrays.");

  len = RARRAY_LEN(array);
  QP_ASSERT(len > 0, "Vectors must have at least one element.");

  *ptr = malloc(len * sizeof(double));
  for (i = 0; i < len; i++) {
    entry = rb_ary_entry(array, i);
    QP_ASSERT(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entries must be Numeric.");
    (*ptr)[i] = NUM2DBL(entry);
  }

  return len;
}

// Allocate a new gsl_matrix and initialize it with the values in array.
// Returns 0 upon success, -1 otherwise. In either case, the caller needs to ensure that
// the memory allocated for vector is freed.
//
int
qp_matrix_from_array(gsl_matrix **pmatrix, const VALUE array) {
  int i, j, nrow, ncol;
  VALUE row, entry;

  QP_ASSERT(rb_obj_is_kind_of(array, rb_cArray) == Qtrue, "Matrices must be Arrays.");
  nrow = RARRAY_LEN(array);
  QP_ASSERT(nrow > 0, "Matrices must have at least one row.");
  row = rb_ary_entry(array, 0);
  QP_ASSERT(rb_obj_is_kind_of(row, rb_cArray) == Qtrue, "Matrix rows must be Arrays.");
  ncol = RARRAY_LEN(row);
  QP_ASSERT(ncol > 0, "Matrices must have at least one column.");

  *pmatrix = gsl_matrix_alloc(nrow, ncol);

  // redundant check that row 0 is an Array
  for (i = 0; i < nrow; i++) {
    row = rb_ary_entry(array, i);
    QP_ASSERT(rb_obj_is_kind_of(row, rb_cArray) == Qtrue, "Matrix rows must be Arrays.");
    QP_ASSERT(RARRAY_LEN(row) == ncol, "Matrix array has rows of different lengths.");

    for (j = 0; j < ncol; j++) {
      entry = rb_ary_entry(row, j);
      QP_ASSERT(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Matrix entries must be Numeric.");
      gsl_matrix_set(*pmatrix, i, j, NUM2DBL(entry));
    }
  }

  return 0;
}

// Allocate a new gsl_matrix and initialize its diagonal entries with values from array.
// Returns 0 upon success, -1 otherwise. In either case, the caller needs to ensure that
// the memory allocated for vector is freed.
//
int
qp_diagonal_matrix_from_array(gsl_matrix **pmatrix, const VALUE array) {
  int i, len;
  VALUE entry;

  QP_ASSERT(rb_obj_is_kind_of(array, rb_cArray) == Qtrue, "Vectors must be Arrays.");

  len = RARRAY_LEN(array);
  QP_ASSERT(len > 0, "Vectors must have at least one element.");

  *pmatrix = gsl_matrix_calloc(len, len);
  for (i = 0; i < len; i++) {
    entry = rb_ary_entry(array, i);
    QP_ASSERT(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entries must be Numeric.");
    gsl_matrix_set(*pmatrix, i, i, NUM2DBL(entry));
  }

  return 0;
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
    VALUE error_message = rb_hash_aref(status_hash, rb_str_new2("message"));
    rb_raise(error_class, StringValueCStr(error_message));
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
  int i, j, nele_jac, nele_hess, arg_status;
  enum ApplicationReturnStatus ipopt_status;
  double *x_var;
  VALUE a_mat, b_vec, x_lower, x_upper, g_mat, g_lower, g_upper, x_init, w_vec;
  qp_minimum_distance_problem *prob;
  IpoptProblem nlp;
  VALUE qp_solution, qp_status, qp_result;

  rb_scan_args(argc, argv, "81", &a_mat, &b_vec, &x_lower, &x_upper, &g_mat, &g_lower, &g_upper, &x_init, &w_vec);

  // Set up the problem parameters and check dimensions of arguments for consistency

  prob = new_qp_minimum_distance_problem();
  x_var = NULL;
  nlp = NULL;

  // printf("trying a_mat... ");
  if (arg_status = qp_matrix_from_array(&prob->Amat, a_mat)) goto finalize_qp_solve_dist_full;
  // printf("ok\n");
  prob->nrow = prob->Amat->size1;
  prob->ncol = prob->Amat->size2;
  // printf("trying b_vec... ");
  if (arg_status = qp_vector_from_array(&prob->bvec, b_vec)) goto finalize_qp_solve_dist_full;
  // printf("ok\n");

  if (prob->Amat->size1 != prob->bvec->size) {
    strncpy(qp_error_buffer, "a_mat.length != b_vec.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  // printf("trying x_lower... ");
  if ((arg_status = qp_ptr_from_array(&prob->x_L, x_lower)) < 0) goto finalize_qp_solve_dist_full;
  // printf("ok\n");

  // arg_status contains x_lower.length
  if (prob->Amat->size2 != arg_status) {
    strncpy(qp_error_buffer, "a_mat[0].length != x_lower.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  if ((arg_status = qp_ptr_from_array(&prob->x_U, x_upper)) < 0) goto finalize_qp_solve_dist_full;

  // arg_status contains x_upper.length
  if (prob->Amat->size2 != arg_status) {
    strncpy(qp_error_buffer, "a_mat[0].length != x_upper.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }
                                  
  if (arg_status = qp_matrix_from_array(&prob->Gmat, g_mat)) goto finalize_qp_solve_dist_full;
  prob->nconstraint = prob->Gmat->size1;

  if (prob->Amat->size2 != prob->Gmat->size2) {
    strncpy(qp_error_buffer, "a_mat[0].length != g_mat[0].length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  if ((arg_status = qp_ptr_from_array(&prob->g_L, g_lower)) < 0) goto finalize_qp_solve_dist_full;

  // arg_status contains g_lower.length
  if (prob->nconstraint != arg_status) {
    strncpy(qp_error_buffer, "g_mat.length != g_lower.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  if ((arg_status = qp_ptr_from_array(&prob->g_U, g_upper)) < 0) goto finalize_qp_solve_dist_full;

  // arg_status contains g_upper.length
  if (prob->nconstraint != arg_status) {
    strncpy(qp_error_buffer, "g_mat.length != g_upper.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  if ((arg_status = qp_ptr_from_array(&x_var, x_init)) < 0) goto finalize_qp_solve_dist_full;

  // arg_status contains x_init.length
  if (prob->Amat->size2 != arg_status) {
    strncpy(qp_error_buffer, "a_mat[0].length != x_init.length", QP_ERROR_BUFFER_LENGTH);
    arg_status = -1;
    goto finalize_qp_solve_dist_full;
  }

  if (argc > 8) {
    if (qp_diagonal_matrix_from_array(&prob->Wmat, w_vec)) goto finalize_qp_solve_dist_full;

    if (prob->Amat->size1 != prob->Wmat->size1) {
      strncpy(qp_error_buffer, "a_mat.length != w_vec.length", QP_ERROR_BUFFER_LENGTH);
      arg_status = -1;
      goto finalize_qp_solve_dist_full;
    }
  } else {
    prob->Wmat = gsl_matrix_calloc(prob->nrow, prob->nrow);
    for (i = 0; i < prob->nrow; i++)
      gsl_matrix_set(prob->Wmat, i, i, 1.0);
  }

  printf("OKOKOKOK\n");

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
  ipopt_status = IpoptSolve(nlp, x_var, NULL, NULL, NULL, NULL, NULL, prob);

  qp_solution = rb_ary_new();
  for (i = 0; i < prob->ncol; i++)
    rb_ary_push(qp_solution, rb_float_new(x_var[i]));

  qp_status = qp_hash_from_status(ipopt_status);
  qp_result = rb_hash_new();
  rb_hash_aset(qp_result, rb_str_new2("solution"), qp_solution);
  rb_hash_aset(qp_result, rb_str_new2("status"), qp_status);

finalize_qp_solve_dist_full:

  if (nlp != NULL) FreeIpoptProblem(nlp);
  if (x_var != NULL) free(x_var);

  free_qp_minimum_distance_problem(prob);

  printf("arg_status: %d\n", arg_status);

  if (arg_status) {
    rb_raise(rb_eArgError, qp_error_buffer);
  }

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
