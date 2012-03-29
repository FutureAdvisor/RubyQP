#include "qp_types.h"

//
// TYPE CONVERSION
//

// Converts an array of Numbers from Ipopt into a GSL vector.
//
QP_STATUS
qp_gsl_vector_set_numbers(
    gsl_vector      *vector,
    const Number    *numbers,
    const size_t    n
    )
{
    QP_STATUS_INIT;
    size_t i;

    for (i = 0; i < n; i++) {
        gsl_vector_set(vector, i, numbers[i]);
    }

qp_cleanup:
    QP_STATUS_RETURN;
}

// Retrieves the length of a non-empty Ruby array.
//
QP_STATUS
qp_rarray_vector_len(
    long        *pn,
    const VALUE rarray
    )
{
    QP_STATUS_INIT;
    long n;

    QP_CHECK(rb_obj_is_kind_of(rarray, rb_cArray) == Qtrue, "Vector must be an Array.");

    n = RARRAY_LEN(rarray);
    QP_CHECK(n > 0, "Vector must have at least one element.");

    *pn = n;

qp_cleanup:
    QP_STATUS_RETURN;
}

// Retrieves the size of a non-empty Ruby array of arrays.
//
QP_STATUS
qp_rarray_matrix_size(
    long        *pnrow,
    long        *pncol,
    const VALUE rarray
    )
{
    QP_STATUS_INIT;
    long nrow, ncol, i;
    VALUE row;

    QP_CHECK(rb_obj_is_kind_of(rarray, rb_cArray) == Qtrue, "Matrix must be an Array.");

    nrow = RARRAY_LEN(rarray);
    QP_CHECK(nrow > 0, "Matrix must have at least one row.");

    row = rb_ary_entry(rarray, 0);
    QP_CHECK(rb_obj_is_kind_of(row, rb_cArray) == Qtrue, "Matrix row must be an Array.");

    ncol = RARRAY_LEN(row);
    QP_CHECK(ncol > 0, "Matrix row must have at least one column.");

    *pnrow = nrow;
    *pncol = ncol;

qp_cleanup:
    QP_STATUS_RETURN;
}

// Converts a Ruby array into a GSL vector.
//
QP_STATUS
qp_gsl_vector_from_rarray(
    gsl_vector  *vector,
    const VALUE rarray,
    const long  n
    )
{
    QP_STATUS_INIT;
    long i;
    VALUE entry;

    for (i = 0; i < n; i++) {
        entry = rb_ary_entry(rarray, i);
        QP_CHECK(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entry must be a Numeric.");
        gsl_vector_set(vector, i, NUM2DBL(entry));
    }

qp_cleanup:
    QP_STATUS_RETURN;
}

// Converts a Ruby array into an array of Numbers for Ipopt.
//
QP_STATUS
qp_numbers_from_rarray(
    Number      *numbers,
    const VALUE rarray,
    const long  n
    )
{
    QP_STATUS_INIT;
    long i;
    VALUE entry;

    for (i = 0; i < n; i++) {
        entry = rb_ary_entry(rarray, i);
        QP_CHECK(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entry must be a Numeric.");
        numbers[i] = NUM2DBL(entry);
    }

qp_cleanup:
    QP_STATUS_RETURN;
}

// Converts a Ruby array of arrays into a GSL matrix.
//
QP_STATUS
qp_gsl_matrix_from_rarray(
    gsl_matrix  *matrix,
    const VALUE rarray,
    const long  nrow,
    const long  ncol
    )
{
    QP_STATUS_INIT;
    long i, j;
    VALUE row, entry;

    for (i = 0; i < nrow; i++) {
        row = rb_ary_entry(rarray, i);
        QP_CHECK(rb_obj_is_kind_of(row, rb_cArray) == Qtrue, "Matrix row must be an Array.");
        QP_CHECK(RARRAY_LEN(row) == ncol, "Matrix has rows of different lengths.");

        for (j = 0; j < nrow; j++) {
            entry = rb_ary_entry(row, j);
            QP_CHECK(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Matrix entry must be a Numeric.");
            gsl_matrix_set(matrix, i, j, NUM2DBL(entry));
        }
    }

qp_cleanup:
    QP_STATUS_RETURN;
}

// Initializes a GSL matrix with values from the specified Ruby array as the diagonal values.
//
QP_STATUS
qp_gsl_matrix_set_diagonal(
    gsl_matrix  *matrix,
    const VALUE rarray,
    const long n
    )
{
    QP_STATUS_INIT;
    long i;
    VALUE entry;

    for (i = 0; i < n; i++) {
        entry = rb_ary_entry(rarray, i);
        QP_CHECK(rb_obj_is_kind_of(entry, rb_cNumeric) == Qtrue, "Vector entry must be a Numeric.");
        gsl_matrix_set(matrix, i, i, NUM2DBL(entry));
    }

qp_cleanup:
    QP_STATUS_RETURN;
}

// Converts a GSL vector into a Ruby array.
//
QP_STATUS
qp_rarray_from_gsl_vector(
    VALUE               *prarray,
    const gsl_vector    *vector
    )
{
    QP_STATUS_INIT;
    long n, i;
    VALUE rarray;

    QP_CHECK(vector->size <= LONG_MAX, "Vector has too many elements.");
    n = (long) vector->size;

    rarray = rb_ary_new2(n);
    for (i = 0; i < n; i++) {
        rb_ary_store(rarray, i, rb_float_new(gsl_vector_get(vector, i)));
    }

    *prarray = rarray;

qp_cleanup:
    QP_STATUS_RETURN;
}

// Converts a GSL matrix into a Ruby array of arrays.
//
QP_STATUS
qp_rarray_from_gsl_matrix(
    VALUE               *prarray,
    const gsl_matrix    *matrix
    )
{
    QP_STATUS_INIT;
    long nrow, ncol, i, j;
    VALUE rarray, row;

    QP_CHECK(matrix->size1 <= LONG_MAX, "Matrix has too many rows.");
    nrow = (long) matrix->size1;

    QP_CHECK(matrix->size2 <= LONG_MAX, "Matrix has too many columns.");
    ncol = (long) matrix->size2;

    rarray = rb_ary_new2(nrow);
    for (i = 0; i < nrow; i++) {
        row = rb_ary_new2(ncol);
        for (j = 0; j < ncol; j++) {
            rb_ary_store(row, j, rb_float_new(gsl_matrix_get(matrix, i, j)));
        }
        rb_ary_store(rarray, i, row);
    }

    *prarray = rarray;

qp_cleanup:
    QP_STATUS_RETURN;
}
