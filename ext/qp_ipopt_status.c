#include "qp_ipopt.h"

// Translates an Ipopt status into a Ruby hash.
//
VALUE
qp_ipopt_status_hash(enum ApplicationReturnStatus ipopt_status) {
    VALUE ipopt_status_hash, status, message, error_class;

    ipopt_status_hash = status = message = error_class = Qnil;

    switch (ipopt_status) {
        case Solve_Succeeded:
            status      = rb_str_new2("Solve_Succeeded");
            message     = rb_str_new2("EXIT: Optimal Solution Found.");
            break;
        case Solved_To_Acceptable_Level:
            status      = rb_str_new2("Solved_To_Acceptable_Level");
            message     = rb_str_new2("EXIT: Solved To Acceptable Level.");
            break;
        case Feasible_Point_Found:
            status      = rb_str_new2("Feasible_Point_Found");
            message     = rb_str_new2("EXIT: Feasible point for square problem found.");
            break;
        case Infeasible_Problem_Detected:
            status      = rb_str_new2("Infeasible_Problem_Detected");
            message     = rb_str_new2("EXIT: Feasible point for square problem found.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("InfeasibleProblemDetectedError"));
            break;
        case Search_Direction_Becomes_Too_Small:
            status      = rb_str_new2("Search_Direction_Becomes_Too_Small");
            message     = rb_str_new2("EXIT: Search Direction is becoming Too Small.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("SearchDirectionBecomesTooSmallError"));
            break;
        case Diverging_Iterates:
            status      = rb_str_new2("Diverging_Iterates");
            message     = rb_str_new2("EXIT: Iterates diverging; problem might be unbounded.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("DivergingIteratesError"));
            break;
        case User_Requested_Stop:
            status      = rb_str_new2("User_Requested_Stop");
            message     = rb_str_new2("EXIT: Sopping optimization at current point as requested by user.");
            break;
        case Maximum_Iterations_Exceeded:
            status      = rb_str_new2("Maximum_Iterations_Exceeded");
            message     = rb_str_new2("EXIT: Maximum Number of Iterations Exceeded.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("MaximumIterationsExceededError"));
            break;
        case Restoration_Failed:
            status      = rb_str_new2("Restoration_Failed");
            message     = rb_str_new2("EXIT: Restoration Failed!");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("RestorationFailedError"));
            break;
        case Error_In_Step_Computation:
            status      = rb_str_new2("Error_In_Step_Computation");
            message     = rb_str_new2("EXIT: Error in step computation (regularization becomes to large?)!");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("ErrorInStepComputationError"));
            break;
        case Invalid_Option:
            status      = rb_str_new2("Invalid_Option");
            message     = rb_str_new2("");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("InvalidOptionError"));
            break;
        case Not_Enough_Degrees_Of_Freedom:
            status      = rb_str_new2("Not_Enough_Degrees_Of_Freedom");
            message     = rb_str_new2("EXIT: Problem has too few degrees of freedom.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("NotEnoughDegreesOfFreedomError"));
            break;
        case Invalid_Problem_Definition:
            status      = rb_str_new2("Invalid_Problem_Definition");
            message     = rb_str_new2("");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("InvalidProblemDefinitionError"));
            break;
        case Unrecoverable_Exception:
            status      = rb_str_new2("Unrecoverable_Exception");
            message     = rb_str_new2("");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("UnrecoverableExceptionError"));
            break;
        case NonIpopt_Exception_Thrown:
            status      = rb_str_new2("NonIpopt_Exception_Thrown");
            message     = rb_str_new2("Unknown Exception caught in Ipopt");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("NonIpoptExceptionThrownError"));
            break;
        case Insufficient_Memory:
            status      = rb_str_new2("Insufficient_Memory");
            message     = rb_str_new2("EXIT: Not enough memory.");
            error_class = rb_eNoMemError;
            break;
        case Internal_Error:
            status      = rb_str_new2("Internal_Error");
            message     = rb_str_new2("EXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.");
            error_class = rb_const_get(rb_mRubyQp, rb_intern("InternalError"));
            break;
    }

    ipopt_status_hash = rb_hash_new();
    rb_hash_aset(ipopt_status_hash, QP_IPOPT_STATUS_HASH_STATUS, status);
    rb_hash_aset(ipopt_status_hash, QP_IPOPT_STATUS_HASH_MESSAGE, message);

    if (!NIL_P(error_class)) {
        rb_hash_aset(ipopt_status_hash, QP_IPOPT_STATUS_HASH_ERROR_CLASS, error_class);
    }

    return ipopt_status_hash;
}

// Raises an exception corresponding to the Ipopt status hash if there was an error.
//
void
qp_handle_ipopt_status_hash(VALUE ipopt_status_hash) {
    VALUE error_class, message;

    error_class = rb_hash_aref(ipopt_status_hash, QP_IPOPT_STATUS_HASH_ERROR_CLASS);
    if (!NIL_P(error_class)) {
        message = rb_hash_aref(ipopt_status_hash, QP_IPOPT_STATUS_HASH_MESSAGE);
        rb_raise(error_class, StringValueCStr(message));
    }
}
