#include <ruby.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include "gsl_cqp.h"

#define QP_MAX_ITER 1000
#define QP_EPS_GAP 1e-10
#define QP_EPS_RESIDUALS 1e-10

#define QP_MATRIX_ENTRY(VALUE,i,j) rb_ary_entry(rb_ary_entry(VALUE,i),j)
#define QP_MATRIX_NROW RARRAY_LEN
#define QP_MATRIX_NCOL(VALUE) RARRAY_LEN(rb_ary_entry(VALUE,0))

VALUE RubyQp = Qnil;

// Raise a ruby Exception if ary can't be interpreted as a GSL vector.
// 
void
qp_ensure_vector(const VALUE ary) {
  int i;

  // A vector argument must be a ruby Array
  Check_Type(ary, T_ARRAY);

  // All entries in the Array must be Floats or Fixnums
  for (i = 0; i < RARRAY_LEN(ary); i++) {
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

  // If it's an empty Array, there's nothing to check
  nrow = RARRAY_LEN(ary);
  if (nrow == 0) return;

  // Make sure the first value in ary is an Array and grab its length
  Check_Type(rb_ary_entry(ary, 0), T_ARRAY);
  ncol = QP_MATRIX_NCOL(ary);

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
      error_class = rb_eRangeError;
      break;
    case GSL_ENOMEM:
      error_class = rb_eNoMemError;
      break;
    default:
      error_class = rb_eException;
  }

  return error_class;
}

// Call the quadratic programming routine on the arguments. The calling function should
// check that the parameters are consistent. Failure to do so may result in an error and
// leave memory unfreed. Returns a ruby Hash with the following keys and values set:
//   "solution"      => minimizing solution
//   "lagrange_eq"   => Lagrange multipliers corresponding to Ax = b
//   "lagrange_ineq" => Lagrange multipliers corresponding to Cx >= d
//   "iterations"    => number of iterations to find the solution
//   "status"        => final status (GSL_SUCCESS or GSL_CONTINUE if maximum number of 
//                      iterations was reached)
//
int
qp_call_cqp(VALUE *qp_result, gsl_matrix *Qmat, gsl_vector *qvec, 
            gsl_matrix *Amat, gsl_vector *bvec, 
            gsl_matrix *Cmat, gsl_vector *dvec) 
{
  gsl_cqp_data *cqp_data;
  size_t iter;
  int status;
  const gsl_cqpminimizer_type *T;
  gsl_cqpminimizer *s;

  cqp_data = malloc(sizeof(gsl_cqp_data));
  cqp_data->Q = Qmat;
  cqp_data->q = qvec;
  cqp_data->A = Amat;
  cqp_data->b = bvec;
  cqp_data->C = Cmat;
  cqp_data->d = dvec;

  T = gsl_cqpminimizer_mg_pdip;
  s = gsl_cqpminimizer_alloc(T, Qmat->size1, Amat->size1, Cmat->size1);
  status = gsl_cqpminimizer_set(s, cqp_data);

  iter = 0;
  do {
    iter++;
    status = gsl_cqpminimizer_iterate(s);
    status = gsl_cqpminimizer_test_convergence(s, QP_EPS_GAP, QP_EPS_RESIDUALS);
  } while (status == GSL_CONTINUE && iter <= QP_MAX_ITER);

  *qp_result = rb_hash_new();
  rb_hash_aset(*qp_result, rb_str_new2("solution"), qp_vector_to_ary(gsl_cqpminimizer_x(s)));
  rb_hash_aset(*qp_result, rb_str_new2("lagrange_eq"), qp_vector_to_ary(gsl_cqpminimizer_lm_eq(s)));
  rb_hash_aset(*qp_result, rb_str_new2("lagrange_ineq"), qp_vector_to_ary(gsl_cqpminimizer_lm_ineq(s)));
  rb_hash_aset(*qp_result, rb_str_new2("iterations"), ULONG2NUM(iter));
  rb_hash_aset(*qp_result, rb_str_new2("status"), INT2FIX(status));

  gsl_cqpminimizer_free(s);
  free(cqp_data);

  return status;
}

// :call-seq:
//   RubyQp::solve_full(q_mat, q_vec, a_mat, b_vec, c_mat, d_vec) => hash
//
// Find +x+ which minimizes the expression
//   (1/2)(x^t)Qx + (q^t)x
// subject to the constraints
//   Ax = b 
//   Cx >= d,
// where +Q+ is the matrix specified by _q_mat_, +q+ is the vector specified by _q_vec_, 
// and so forth for +A+, +b+, +C+, and +d+. _q_mat_ (and all matrix arguments) should be 
// an Array of Arrays. Sub-Arrays should all have the same length, and their entries must 
// be Numeric. _q_vec_ (and all vector arguments) should be an Array of Numerics. 
//
// The constraint +Cx = d+ means that each entry in the vector +Cx+ is greater than or
// equal to the corresponding entry in the vector +d+.
//
// Returns a ruby Hash with the following keys and values set:
//   "solution"      => minimizing solution
//   "lagrange_eq"   => Lagrange multipliers corresponding to Ax = b
//   "lagrange_ineq" => Lagrange multipliers corresponding to Cx >= d
//   "iterations"    => number of iterations to find the solution
//   "status"        => final status (GSL_SUCCESS or GSL_CONTINUE if maximum number of 
//                      iterations was reached)
//
VALUE
qp_solve_full(VALUE self, VALUE Qary, VALUE qary, VALUE Aary, VALUE bary, 
              VALUE Cary, VALUE dary) 
{
  // copy inputs into GSL matrix and vector types, initialize CQP data
  int ncol, status;
  gsl_matrix *Qmat, *Amat, *Cmat;
  gsl_vector *qvec, *bvec, *dvec;
  VALUE qp_result;

  qp_ensure_matrix(Qary);
  qp_ensure_matrix(Aary);
  qp_ensure_matrix(Cary);
  qp_ensure_vector(qary);
  qp_ensure_vector(bary);
  qp_ensure_vector(dary);

  // Check argument dimensions
  ncol = QP_MATRIX_NCOL(Qary);
  if (RARRAY_LEN(qary) != ncol ||
      QP_MATRIX_NCOL(Aary) != ncol ||
      RARRAY_LEN(bary) != QP_MATRIX_NROW(Aary) ||
      QP_MATRIX_NCOL(Cary) != ncol ||
      RARRAY_LEN(dary) != QP_MATRIX_NROW(Cary)) 
  {
    rb_raise(rb_eArgError, "Arguments have inconsistent dimensions");
  }

  Qmat = qp_ary_to_matrix(Qary);
  Amat = qp_ary_to_matrix(Aary);
  Cmat = qp_ary_to_matrix(Cary);
  qvec = qp_ary_to_vector(qary);
  bvec = qp_ary_to_vector(bary);
  dvec = qp_ary_to_vector(dary);

  status = qp_call_cqp(&qp_result, Qmat, qvec, Amat, bvec, Cmat, dvec);

  gsl_matrix_free(Qmat);
  gsl_matrix_free(Amat);
  gsl_matrix_free(Cmat);
  gsl_vector_free(qvec);
  gsl_vector_free(bvec);
  gsl_vector_free(dvec);

  if (status) {
    rb_raise(qp_error_class_for(status), gsl_strerror(status));
  }

  return qp_result;
}

// :call-seq: 
//   RubyQp::solve(q_mat, q_vec, a_mat, b_vec, c_mat, d_vec) => array
//
// Same as RubyQp::solve_full, but returns only the solution vector as a ruby Array.
//
VALUE
qp_solve(VALUE self, VALUE Qary, VALUE qary, VALUE Aary, VALUE bary, 
              VALUE Cary, VALUE dary) 
{
  VALUE qp_result = qp_solve_full(self, Qary, qary, Aary, bary, Cary, dary);
  return rb_hash_aref(qp_result, rb_str_new2("solution"));
}

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
//
// Returns a ruby Hash with the following keys and values set:
//   "solution"      => minimizing solution
//   "lagrange_eq"   => Lagrange multipliers corresponding to Ax = b
//   "lagrange_ineq" => Lagrange multipliers corresponding to Cx >= d
//   "iterations"    => number of iterations to find the solution
//   "status"        => final status (GSL_SUCCESS or GSL_CONTINUE if maximum number of 
//                      iterations was reached)
//
VALUE
qp_solve_dist_full(int argc, VALUE *argv, VALUE self) {
  VALUE Mary, mary, Aary, bary, Cary, dary, wary;
  gsl_matrix *Mmat, *Amat, *Cmat, *Wmat, *Qmat, *Tempmat;
  gsl_vector *mvec, *bvec, *dvec, *qvec, *tempvec;
  int i, j, len, status;
  VALUE qp_result;

  rb_scan_args(argc, argv, "61", &Mary, &mary, &Aary, &bary, &Cary, &dary, &wary);

  qp_ensure_matrix(Mary);
  qp_ensure_matrix(Aary);
  qp_ensure_matrix(Cary);
  qp_ensure_vector(mary);
  qp_ensure_vector(bary);
  qp_ensure_vector(dary);

  // Check argument dimensions
  if (QP_MATRIX_NROW(Mary) != RARRAY_LEN(mary) ||
      QP_MATRIX_NCOL(Mary) != QP_MATRIX_NCOL(Aary) ||
      QP_MATRIX_NROW(Aary) != RARRAY_LEN(bary) ||
      QP_MATRIX_NCOL(Mary) != QP_MATRIX_NCOL(Cary) ||
      QP_MATRIX_NROW(Cary) != RARRAY_LEN(dary))
  {
    rb_raise(rb_eArgError, "Arguments have inconsistent dimensions");
  }

  if (argc > 6) {
    qp_ensure_vector(wary);
    if (QP_MATRIX_NROW(Mary) != RARRAY_LEN(wary)) {
      rb_raise(rb_eArgError, "Weight vector length inconsistent with dimension of M");
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
      gsl_matrix_set(Wmat, i, i, NUM2DBL(rb_ary_entry(wary, i)));
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

  if (status) {
    rb_raise(qp_error_class_for(status), gsl_strerror(status));
  }

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

void Init_ruby_qp(void) {
  RubyQp = rb_define_module("RubyQp");
  rb_define_module_function(RubyQp, "solve_full", qp_solve_full, 6);
  rb_define_module_function(RubyQp, "solve", qp_solve, 6);
  rb_define_module_function(RubyQp, "solve_dist_full", qp_solve_dist_full, -1);
  rb_define_module_function(RubyQp, "solve_dist", qp_solve_dist, -1);

  gsl_set_error_handler_off();
}
