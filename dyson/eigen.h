/* Eigenvector and eigenvalue functions.
 *
 * Copyright Mark Burnett, April 2009
 */
#ifndef _NUMERIC_LINALG_EIGEN_H_
#define _NUMERIC_LINALG_EIGEN_H_

#include <exception>
#include <types.h>
#include <boost/foreach.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/bindings/lapack/lapack.h>
//#include <lapack.h>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>

#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>

#include <utility>
#include <vector>
#include <algorithm>
#include <iostream>

namespace ublas = boost::numeric::ublas;

namespace numeric { namespace linalg {

// Real versions
cvector_t eigenvalues( const matrix_t &mat );
std::pair< cvector_t, matrix_t > eig( const matrix_t &mat );

// Complex versions
cvector_t eigenvalues( const cmatrix_t &mat );
std::pair< cvector_t, cmatrix_t > eig( const cmatrix_t &mat );

// Sorted eigenvalue versions
std::vector< double >
sorted_eigenvalues( const matrix_t &m );

void mult_mtx(const cmatrix_t &A, const cmatrix_t &B, cmatrix_t &C);

// Complex Matrix Inversion
void 
inverse_mtx( cmatrix_t &m );

// Real Matrix Inversion
void 
inverse_mtx( matrix_t &m );

} // end namespace linalg
} // end namespace numeric

#endif // _NUMERIC_LINALG_EIGEN_H_
