#ifndef QP_IPOPT_H
#define QP_IPOPT_H

#include "ruby_qp.h"

#include <IpStdCInterface.h>

//
// CALLBACKS
//

// Eval_F_CB
Bool
qp_eval_f(Index n, Number *x, Bool new_x,
          Number *obj_value, UserDataPtr user_data);

// Eval_Grad_F_CB
Bool
qp_eval_grad_f(Index n, Number *x, Bool new_x,
               Number *grad_f, UserDataPtr user_data);

// Eval_G_CB
Bool
qp_eval_g(Index n, Number *x, Bool new_x,
          Index m, Number *g, UserDataPtr user_data);

// Eval_Jac_G_CB
Bool
qp_eval_jac_g(Index n, Number *x, Bool new_x,
              Index m, Index nele_jac,
              Index *iRow, Index *jCol, Number *values,
              UserDataPtr user_data);

// Eval_H_CB
Bool
qp_eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
          Index m, Number *lambda, Bool new_lambda,
          Index nele_hess, Index *iRow, Index *jCol,
          Number *values, UserDataPtr user_data);

//
// PROBLEM CREATION
//

#define QP_IPOPT_PROBLEM_CREATE(nlp, n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style)                                                                \
    nlp = CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style, &qp_eval_f, &qp_eval_g, &qp_eval_grad_f, &qp_eval_jac_g, &qp_eval_h);  \
    QP_CLEANUP_ON_NULL(nlp);

#define QP_IPOPT_PROBLEM_FREE(nlp) if (NULL != nlp) FreeIpoptProblem(nlp);

//
// ERROR HANDLING
//

// Maximum value for an Ipopt Index variable.
#define QP_IPOPT_INDEX_MAX INT_MAX

// If the following line is uncommented, then a Ruby exception may be raised depending
// on the IpoptSolve return status. Otherwise, the return status is passed through in the
// result hash.
// #define QP_RAISE_ON_ERROR 1

// Ruby symbols used to access values in the Ruby hash representing an Ipopt status.
#define QP_IPOPT_STATUS_HASH_STATUS         ID2SYM(rb_intern("status"))
#define QP_IPOPT_STATUS_HASH_MESSAGE        ID2SYM(rb_intern("message"))
#define QP_IPOPT_STATUS_HASH_ERROR_CLASS    ID2SYM(rb_intern("error_class"))

VALUE
qp_ipopt_status_hash(enum ApplicationReturnStatus ipopt_status);

void
qp_handle_ipopt_status_hash(VALUE ipopt_status_hash);

#endif  // #ifndef QP_IPOPT_H
