#include "ruby_qp.h"
#include "qp_types.h"

VALUE rb_mRubyQp = Qnil;
char qp_error_buffer[QP_ERROR_BUFFER_LENGTH + 1];

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
// _x_init_ is the initial point for the optimization algorithm.
//
// Returns a Ruby Hash with the following keys and values set:
//   :solution  => minimizing solution
//   :status    => Ipopt status code and message
//
VALUE
qp_solve_dist_full(int argc, VALUE *argv, VALUE self) {
    QP_STATUS_INIT;
    VALUE a_mat, b_vec, x_lower, x_upper, g_mat, g_lower, g_upper, x_init, w_vec, ipopt_status_hash, solution, result;
    qp_minimum_distance_problem *prob;
    long nrowa, ncola, nb, nx_lower, nx_upper, nrowg, ncolg, ng_lower, ng_upper, nx_init, nw, i, j;
    Index nvar, nconstraint, nele_jac, nele_hess;
    Number *x_L, *x_U, *g_L, *g_U, *x_var, minimum_distance;
    IpoptProblem nlp;
    enum ApplicationReturnStatus ipopt_status;

    rb_scan_args(argc, argv, "81", &a_mat, &b_vec, &x_lower, &x_upper, &g_mat, &g_lower, &g_upper, &x_init, &w_vec);
    ipopt_status_hash = solution = result = Qnil;
    prob = NULL;
    nrowa = ncola = nb = nx_lower = nx_upper = nrowg = ncolg = ng_lower = ng_upper = nx_init = nw = i = j = 0;
    nvar = nconstraint = nele_jac = nele_hess = 0;
    x_L = x_U = g_L = g_U = x_var = NULL;
    minimum_distance = 0;
    nlp = NULL;
    ipopt_status = Solve_Succeeded;

    // Set up the problem parameters and check dimensions of arguments for consistency.
    QP_MININUM_DISTANCE_PROBLEM_ALLOC(prob);

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_matrix_size(&nrowa, &ncola, a_mat));
    QP_CHECK(ncola <= QP_IPOPT_INDEX_MAX, "Too many optimization variables.");
    nvar = (Index) ncola;
    QP_GSL_MATRIX_ALLOC(prob->Amat, nrowa, ncola);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_gsl_matrix_from_rarray(prob->Amat, a_mat, nrowa, ncola));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&nb, b_vec));
    QP_CHECK(nrowa == nb, "a_mat.length != b_vec.length")
    QP_GSL_VECTOR_ALLOC(prob->bvec, nb);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_gsl_vector_from_rarray(prob->bvec, b_vec, nb));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&nx_lower, x_lower));
    QP_CHECK(ncola == nx_lower, "a_mat[0].length != x_lower.length");
    QP_MALLOC_N(x_L, Number, nx_lower);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_numbers_from_rarray(x_L, x_lower, nx_lower));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&nx_upper, x_upper));
    QP_CHECK(ncola == nx_upper, "a_mat[0].length != x_upper.length");
    QP_MALLOC_N(x_U, Number, nx_upper);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_numbers_from_rarray(x_U, x_upper, nx_upper));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_matrix_size(&nrowg, &ncolg, g_mat));
    QP_CHECK(ncola == ncolg, "a_mat[0].length != g_mat[0].length");
    QP_CHECK(nrowg <= QP_IPOPT_INDEX_MAX, "Too many constraints.");
    QP_CHECK(ncolg <= QP_IPOPT_INDEX_MAX / nrowg, "Too many elements in g_mat.");  // nele_jac can potentially overflow later if this condition is not satisfied.
    nconstraint = (Index) nrowg;
    QP_GSL_MATRIX_ALLOC(prob->Gmat, nrowg, ncolg);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_gsl_matrix_from_rarray(prob->Gmat, g_mat, nrowg, ncolg));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&ng_lower, g_lower));
    QP_CHECK(nrowg == ng_lower, "g_mat.length != g_lower.length");
    QP_MALLOC_N(g_L, Number, ng_lower);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_numbers_from_rarray(g_L, g_lower, ng_lower));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&ng_upper, g_upper));
    QP_CHECK(nrowg == ng_upper, "g_mat.length != g_upper.length");
    QP_MALLOC_N(g_U, Number, ng_upper);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_numbers_from_rarray(g_U, g_upper, ng_upper));

    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&nx_init, x_init));
    QP_CHECK(ncola == nx_init, "a_mat[0].length != x_init.length");
    QP_MALLOC_N(x_var, Number, nx_init);
    QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_numbers_from_rarray(x_var, x_init, nx_init));

    QP_GSL_MATRIX_CALLOC(prob->Wmat, nrowa, nrowa);
    if (NIL_P(w_vec)) {
        for (i = 0; i < nrowa; i++) {
            gsl_matrix_set(prob->Wmat, i, i, 1.0);
        }
    } else {
        QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_rarray_vector_len(&nw, w_vec));
        QP_CHECK(nrowa == nw, "a_mat.length != w_vec.length");
        QP_CLEANUP_WITH_STATUS_ON_FAIL(qp_gsl_matrix_set_diagonal(prob->Wmat, w_vec, nw));
    }

    // Count the number of nonzero elements in the constraint Jacobian
    // (i.e., the number of nonzero elements in Gmat).
    for (i = 0; i < nrowg; i++) {
        for (j = 0; j < ncolg; j++) {
            if (gsl_matrix_get(prob->Gmat, i, j) != 0) {
                nele_jac++;
            }
        }
    }

    // Our Hessian is dense, so nele_hess is just the number of entries below the diagonal
    // (inclusive) in the Hessian.
    for (j = 0; j < ncolg; j++) {
        nele_hess += j;
    }

    QP_IPOPT_PROBLEM_CREATE(nlp, nvar, x_L, x_U, nconstraint, g_L, g_U, nele_jac, nele_hess, 0);

    // TODO: These options should be set globally in a ruby-accessible manner.
    AddIpoptIntOption(nlp, "print_level", 0);
    AddIpoptIntOption(nlp, "max_iter", 200);
    AddIpoptStrOption(nlp, "mehrotra_algorithm", "yes");

    // TODO: Capture more of the output.
    ipopt_status = IpoptSolve(nlp, x_var, NULL, &minimum_distance, NULL, NULL, NULL, prob);
    ipopt_status_hash = qp_ipopt_status_hash(ipopt_status);

    solution = rb_ary_new2(nx_init);
    for (i = 0; i < nx_init; i++) {
        rb_ary_store(solution, i, rb_float_new(x_var[i]));
    }

    result = rb_hash_new();
    rb_hash_aset(result, QP_RESULT_SOLUTION, solution);
    rb_hash_aset(result, QP_RESULT_MINIMUM_DISTANCE, rb_float_new(minimum_distance));
    rb_hash_aset(result, QP_RESULT_STATUS, ipopt_status_hash);

qp_cleanup:
    QP_MINIMUM_DISTANCE_PROBLEM_FREE(prob);
    QP_FREE(x_L);
    QP_FREE(x_U);
    QP_FREE(g_L);
    QP_FREE(g_U);
    QP_FREE(x_var);
    QP_IPOPT_PROBLEM_FREE(nlp);

    if (QP_STATUS_FAILED) {
        rb_raise(rb_eArgError, qp_error_buffer);
    }

// TODO: Errors from Ipopt are hidden since the QP_RAISE_ON_ERROR definition is normally commented out; figure out a better way to handle errors.
#ifdef QP_RAISE_ON_ERROR
    qp_handle_ipopt_status_hash(ipopt_status_hash);
#endif

    return result;
}

// :call-seq:
//   RubyQp::solve_dist(m_mat, m_vec, a_mat, b_vec, c_mat, d_vec, w_vec = nil) => array
//
// Same as RubyQp::solve_dist_full, but returns only the solution vector as a ruby Array.
//
VALUE
qp_solve_dist(int argc, VALUE *argv, VALUE self) {
    VALUE result = qp_solve_dist_full(argc, argv, self);
    return rb_hash_aref(result, QP_RESULT_SOLUTION);
}

void 
Init_ruby_qp(void) {
    rb_mRubyQp = rb_define_module("RubyQp");
    memset(qp_error_buffer, 0, sizeof(qp_error_buffer));

    rb_require("ruby_qp/errors");
    rb_define_module_function(rb_mRubyQp, "solve_dist_full", qp_solve_dist_full, -1);
    rb_define_module_function(rb_mRubyQp, "solve_dist", qp_solve_dist, -1);
}
