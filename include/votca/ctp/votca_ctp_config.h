#ifndef _VOTCA_CTP_CONFIG_H
#define _VOTCA_CTP_CONFIG_H

/* Linear algebra packages */
/* #undef EIGEN */
#define MKL
#define GSL

#define NDEBUG

#if defined(MKL)
    #include <votca/ctp/mkl_boost_ublas_matrix_prod.h>
#elif defined(GSL)
    #include <votca/ctp/gsl_boost_ublas_matrix_prod.h>
#endif

#endif // _VOTCA_CTP_CONFIG_H