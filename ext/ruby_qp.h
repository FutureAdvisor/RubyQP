#ifndef RUBY_QP_H
#define RUBY_QP_H

#include <ruby.h>

extern VALUE rb_mRubyQp;
extern char qp_error_buffer[];
#define QP_ERROR_BUFFER_LENGTH 256

typedef int QP_STATUS;
#define QP_SUCCESS  (0)
#define QP_FAIL     (-1)

#define QP_FAILED(STATUS)   ((STATUS) < 0)
#define QP_SUCCEEDED        !QP_FAILED

#define QP_STATUS_INIT          QP_STATUS qp_status = QP_SUCCESS;
#define QP_STATUS_SET(STATUS)   (qp_status = STATUS)
#define QP_STATUS_FAILED        QP_FAILED(qp_status)
#define QP_STATUS_SUCCEEDED     QP_SUCCEEDED(qp_status)
#define QP_STATUS_RETURN        return qp_status;

#define QP_CLEANUP                              goto qp_cleanup;
#define QP_CLEANUP_ON_NULL(x)                   if (NULL == x) QP_CLEANUP;
#define QP_CLEANUP_ON_FAIL(STATUS)              if (QP_FAILED(STATUS)) QP_CLEANUP;
#define QP_CLEANUP_WITH_STATUS_ON_FAIL(STATUS)  QP_CLEANUP_ON_FAIL(QP_STATUS_SET(STATUS))
    // Similar to QP_CLEANUP_ON_FAIL, but saves the status code in qp_status.

#define QP_CHECK(TEST, ERROR_MESSAGE)                                       \
    if (!(TEST)) {                                                          \
        strncpy(qp_error_buffer, ERROR_MESSAGE, QP_ERROR_BUFFER_LENGTH);    \
        QP_STATUS_SET(QP_FAIL);                                             \
        QP_CLEANUP;                                                         \
    }

#define QP_MALLOC(x, type)              \
    x = (type *) malloc(sizeof(type));  \
    QP_CLEANUP_ON_NULL(x);              \
    memset(x, 0, sizeof(type));

#define QP_MALLOC_N(x, type, n)             \
    x = (type *) malloc(sizeof(type) * n);  \
    QP_CLEANUP_ON_NULL(x);                  \
    memset(x, 0, sizeof(type) * n);

#define QP_FREE(x) if (NULL != x) free(x);

// Ruby symbols used to access values in the Ruby hash representing the solver's result.
#define QP_RESULT_SOLUTION          ID2SYM(rb_intern("solution"))
#define QP_RESULT_MINIMUM_DISTANCE  ID2SYM(rb_intern("minimum_distance"))
#define QP_RESULT_STATUS            ID2SYM(rb_intern("status"))

#endif  // #ifndef RUBY_QP_H
