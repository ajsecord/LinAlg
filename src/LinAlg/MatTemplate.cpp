/*
 * LinAlg: A fixed-length vector and matrix library.
 * Mat$D.cpp
 * Copyright 2005-2008 Adrian Secord.
 * Contact: <a href="http://www.google.com/search?q=%22Adrian%20Secord%22">Adrian Secord</a>
 */

#include <iostream>
#include <limits>
#include <sstream>

#if $D == 3
#   include <vector>
#endif

#include <LinAlg/Mat$D.h>
#include <LinAlg/Vec$D.h>

// This helper class that includes 2D and 3D specializations 
// of some matrix functions and float/double specialiazations of Lapack
// calls.
#include "TypeHelpers.h"

#ifndef M_PI
#   define M_PI 3.1415926535897932384626433832795
#endif

namespace LinAlg {
    
    // Private namespace for this file
namespace {
    
    // This helper class will be specialized for floats or doubles
    // to farm off calls to Lapack, which requires differently-named
    // functions for the different types.
    // This could be extended to complex numbers if necessary.
    template <typename T>
    struct LapackHelper$D {
        static void eigs(const Mat$D<T>& A, 
                         Vec$D<T>& values, Mat$D<T>& vectors) {
            assert(!"Cannot calculate eigenvalues for this type");
        }
    };
    
#if !defined(LINALG_NO_LAPACK)
    
#   include "Lapack.h"
    
    /// Get the optimum block size for this machine.  Note that we should only 
    /// call this routine once, and after that it should just return the 
    /// precalculated block size.
    int getBlockSize(const int order) {
        int param = 1;
        int negOne = -1;
        char name[] = "DGETRI";
        char opts[] = "";
        int o = order;
        
        static int blockSize = ilaenv_(&param, name, opts, &o, &o, 
                                       &negOne, &negOne);
        return blockSize;
    }
    
    // Specialize for float types
    template <>
    struct LapackHelper$D<float> {
        static void eigs(const Mat$D<float>& A, 
                         Vec$D<float>& values, Mat$D<float>& vectors) {
            char jobz = 'V';          // Compute both eigenvalues and eigenvectors
            char uplo = 'L';          // Lower triangle of A (at least) is stored
            int n     = $D;            // Order of the matrix
            int lda   = $D;            // Leading dimension
            int lwork = (getBlockSize($D)+2)*$D; // Size of additional work space
            int info;
            
            float* work = new float[lwork];
            
            vectors = A;
            
            ssyev_(&jobz, &uplo, &n, vectors.data(), &lda, 
                   values.data(), work, &lwork, &info);
            
            assert(info == 0);
            
            delete [] work;
        }
    };
    
    // Specialize for double types
    template <>
    struct LapackHelper$D<double> {
        static void eigs(const Mat$D<double>& A, 
                         Vec$D<double>& values, Mat$D<double>& vectors) {
            char jobz = 'V';          // Compute both eigenvalues and eigenvectors
            char uplo = 'L';          // Lower triangle of A (at least) is stored
            int n     = $D;            // Order of the matrix
            int lda   = $D;            // Leading dimension
            int lwork = (getBlockSize($D)+2)*$D; // Size of additional work space
            int info;
            
            double* work = new double[lwork];
            
            vectors = A;
            
            dsyev_(&jobz, &uplo, &n, vectors.data(), &lda, 
                   values.data(), work, &lwork, &info);
            
            assert(info == 0);
            
            delete [] work;
        }
    };
#endif
    
    // Dimension-specific eigenvalue routines
#if $D == 2

    /// Matrix operations that only make sense for floating-point types.
    /// This will get specialized to only handle floating-point types.
    template <typename T, bool IntegerType = std::numeric_limits<T>::is_integer>
    struct Mat2FloatingPointHelper {
        static void eigs(const Mat2<T>& a, Vec2<T>& values, Mat2<T>& vectors) {
            assert(!"Taking the eigenvalues of an integer array doesn't make sense.  Convert to floating point if required.");
        }
    };
    
    template <typename T>
    struct Mat2FloatingPointHelper<T, false> {
        static void eigs(const Mat2<T>& a, Vec2<T>& values, Mat2<T>& vectors) {
            // Symmetricity check
            const T eps = T(1e3) * std::numeric_limits<T>::epsilon();
            assert(std::abs(a(0,1) - a(1,0)) < eps);
            
            // Check for already-diagonal matrices.  The following code does not
            // deal well with them.
            if (std::abs(a(0,0)) >= eps && std::abs(a(1,0)) <  eps &&
                std::abs(a(0,1)) <  eps && std::abs(a(1,1)) >= eps) {
                if (std::abs(a(0,0)) > std::abs(a(1,1))) {
                    values(0) = a(0,0);
                    values(1) = a(1,1);
                    vectors(0,0) = 1; vectors(1,0) = 0;
                    vectors(0,1) = 0; vectors(1,1) = 1;
                } else {
                    values(0) = a(1,1);
                    values(1) = a(0,0);
                    vectors(0,0) = 0; vectors(1,0) = 1;
                    vectors(0,1) = 1; vectors(1,1) = 0;
                }
                
                return;
            }
            
            // Solve the equation det(a - gamma * I) = 0
            const T b = -a(0,0) - a(1,1);
            const T c = a(0,0) * a(1,1) - a(0,1) * a(1,0);
            
            // Solve the quadratic equation gamma^2 + b * gamma + c = 0
            // to get the evals.
            // From Numerical Recipes in C, 2nd. ed.
            const T sgn_b = (T)(b < 0 ? -1 : 1);
            const T sqrt_arg = b * b - 4 * c;
            
            // In the case of two identical evals, sqrt_arg should be zero, but
            // can be numerically very small and negative, which is bad.
            T q;
            if (std::abs(sqrt_arg) >= eps)
                q = (T)(-0.5 * (b + sgn_b * sqrt(sqrt_arg)));
            else 
                q = (T)(-0.5 * b);
            
            values(0) = q;
            values(1) = c / q;
            
            // Make the first eval the dominant one.
            if (values(0) < values(1)) {
                const T tmp = values(0);
                values(0) = values(1);
                values(1) = tmp;
            }
            
            // Arbitrary choice for the first evect -- choose so 
            // the length of the evect as large as possible...
            if (std::abs(values(0) - a(0,0)) > std::abs(values(0) - a(1,1))) {
                vectors(0,0) = a(1,0);
                vectors(1,0) = values(0) - a(0,0); 
            } else {
                vectors(0,0) = values(0) - a(1,1);
                vectors(1,0) = a(1,0);
            }
            
            // ... because we are normalising the vectors.
            const T length = sqrt(vectors(0,0) * vectors(0,0) + 
                                  vectors(1,0) * vectors(1,0));
            //assert(std::abs(length) > eps);
            assert(length != 0.0);
            
            vectors(0,0) /= length;
            vectors(1,0) /= length;
            
            // Choose 2nd evect to be perpendicular to the first;
            vectors(0,1) = -vectors(1,0);
            vectors(1,1) =  vectors(0,0);
        }
    };
    
#elif $D == 3
    /// Find an eigenvector corresponding to an eigenvalue of multiplicity 1;
    template <typename T> 
    void findEigVect1(const Mat3<T>& A, const T value, Vec3<T>& vec1) {
        Mat3<T> M;
        
        M = A - value * Mat3<T>::ident();
        
        static Mat3<T> U;
        U(0,0) = M(1,1) * M(2,2) - M(1,2) * M(1,2);
        U(1,0) = M(0,2) * M(1,2) - M(0,1) * M(2,2);
        U(2,0) = M(0,1) * M(1,2) - M(0,2) * M(1,1);
        
        U(0,1) = U(1,0);
        U(1,1) = M(0,0) * M(2,2) - M(0,2) * M(0,2);
        U(2,1) = M(0,2) * M(0,1) - M(1,2) * M(0,0);
        
        U(0,2) = U(2,0);
        U(1,2) = U(2,1);
        U(2,2) = M(0,0) * M(1,1) - M(0,1) * M(0,1);
        
        // Find max element,  and its column.
        int maxCol = 0;
        T maxElem = -std::numeric_limits<T>::max();
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                if (U(i,j) > maxElem) {
                    maxElem = U(i,j);
                    maxCol = j;
                }
            }
        }
        
        vec1 = U.getCol(maxCol);
        vec1 /= sqrt(inner(vec1, vec1));
    }
    
    /// Find an eigenvector corresponding to an eigenvalue of multiplicity 2.
    template <typename T> 
    void findEigVect2(const Mat3<T>& A, const T value, 
                      Vec3<T>& vec1, Vec3<T>& vec2) {
        static Mat3<T> M;
        
        M = A - value * Mat3<T>::ident();
        
        // Find max element, its row and column
        int maxCol = 0, maxRow = 0;
        T maxElem = -std::numeric_limits<T>::max();
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                if (M(i,j) > maxElem) {
                    maxElem = M(i,j);
                    maxRow = i;
                    maxCol = j;
                }
            }
        }
        
        if ((maxRow == 0 && maxCol == 0) ||
            (maxRow == 0 && maxCol == 1) ||
            (maxRow == 1 && maxCol == 0)) {
            vec1(0) = -M(0,1);
            vec1(1) =  M(0,0);
            vec1(2) =  0.0;
            
            vec2(0) = -M(0,2) * M(0,0);
            vec2(1) = -M(0,2) * M(0,1);
            vec2(2) =  M(0,0) * M(0,0) + M(0,1) * M(0,1);
            
        } else if ((maxRow == 0 && maxCol == 2) ||
                   (maxRow == 2 && maxCol == 0)) {
            vec1(0) =  M(0,2);
            vec1(1) =  0.0;
            vec1(2) = -M(0,0);
            
            vec2(0) = -M(0,1) * M(0,0);
            vec2(1) =  M(0,0) * M(0,0) + M(0,2) * M(0,2);
            vec2(2) = -M(0,1) * M(0,2);
            
        } else if ((maxRow == 1 && maxCol == 1) ||
                   (maxRow == 1 && maxCol == 2) ||
                   (maxRow == 2 && maxCol == 1)) {
            vec1(0) =  0.0;
            vec1(1) = -M(1,2);
            vec1(2) =  M(1,1);
            
            vec2(0) =  M(1,1) * M(1,1) + M(1,2) * M(1,2);
            vec2(1) = -M(0,1) * M(1,1);
            vec2(2) = -M(0,1) * M(1,2);
            
        } else if (maxRow == 2 && maxCol == 2) {
            vec1(0) =  0.0;
            vec1(1) = -M(2,2);
            vec1(2) =  M(1,2);
            
            vec2(0) =  M(1,2) * M(1,2) + M(2,2) * M(2,2);
            vec2(1) = -M(0,2) * M(1,2);
            vec2(2) = -M(0,2) * M(2,2);
            
        } else {
            assert(false);
        }
        
        vec1 /= sqrt(inner(vec1, vec1));
        vec2 /= sqrt(inner(vec2, vec2));
    }    
    
    /// Find the real roots of the cubic equation
    /// x^3 + c(2) * x^2 + c(1) * x + c(0) = 0
    template <typename T>  
    void cubicRoots(const Vec3<T>& c, std::vector<T>& roots) {
        const T eps = T(1e3) * std::numeric_limits<T>::epsilon();
        
        // Solve the cubic for lambda, from Numerical Recipes in C, 2nd ed.
        const T Q = (c(2)*c(2) - 3*c(1)) / 9;
        const T R = (2*c(2)*c(2)*c(2) - 9*c(2)*c(1) + 27*c(0)) / 54;
        
        const T R_2 = R*R;
        const T Q_3 = Q*Q*Q;
        
        // If R^2 < Q^3, then we have three real roots
        if (R_2 < Q_3) {
            const T theta = acos(R / sqrt(Q_3));
            
            // The three roots
            roots.resize(3);
            roots[0] = T(-2 * sqrt(Q) * cos(theta/3)            - c(2) / 3);
            roots[1] = T(-2 * sqrt(Q) * cos((theta + 2*M_PI)/3) - c(2) / 3);
            roots[2] = T(-2 * sqrt(Q) * cos((theta - 2*M_PI)/3) - c(2) / 3);
            
        } else {
            
            const T sgn_R = (R < 0 ? T(-1) : T(1));
            const T A = -sgn_R * pow(T(std::abs(R) + sqrt(R_2 - Q_3)), T(1/3.0));
            const T B = (std::abs(A) < eps ? 0 : Q / A);
            
            // Okay, maybe two roots.
            if (std::abs(A - B) < eps) {
                roots.resize(2);
                roots[0] = (A + B) - c(2)/3;
                roots[1] = T(-0.5) * (A + B) - c(2) / 3;
                
                // Only one root
            } else {
                roots.resize(1);
                roots[0] = (A + B) - c(2)/3;
            }
        }
    }
        
    template <typename T, bool IntegerType = std::numeric_limits<T>::is_integer>
    struct Mat3FloatingPointHelper {
        static void eigs(const Mat3<T>& A, 
                         Vec3<T>& values, Mat3<T>& vectors) {
            assert(!"Taking the eigenvalues of an integer array doesn't make sense.  Convert to floating point if required.");
        }
    };
    
    template <typename T>
    struct Mat3FloatingPointHelper<T, false> {
        static void eigs(const Mat3<T>& A, Vec3<T>& values, Mat3<T>& vectors) {
            
            // The coefficients of the cubic characteristic equation
            // lambda^3 + c2*lambda^2 + c1*lambda + c0 = 0
            static Vec3<T> c;
            c(2) = - A(0,0) - A(1,1) - A(2,2);
            c(1) =   A(0,0) * A(1,1) + A(0,0) * A(2,2) + A(1,1) * A(2,2) 
                   -(A(1,0) * A(1,0) + A(2,1) * A(2,1) + A(2,0) * A(2,0));
            c(0) = - A(0,0) * A(1,1) * A(2,2) 
                   - T(2.0) * A(1,0) * A(2,1) * A(2,0) 
                   + A(0,0) * A(2,1) * A(2,1) 
                   + A(1,1) * A(2,0) * A(2,0) 
                   + A(2,2) * A(1,0) * A(1,0);
            
            std::vector<T> roots;
            cubicRoots(c, roots);
            //sort3(values);
            
            static Vec3<T> vec1, vec2, vec3;
            if (roots.size() == 1) {
                // A = lambda * I  
                values(0) = roots[0];
                values(1) = roots[0];
                values(2) = roots[0];
                vectors = Mat3<T>::ident();
                
            } else if (roots.size() == 2) {
                // First one is of multiplicity one.
                findEigVect1(A, roots[0], vec1);
                
                // Second one is of multiplicity two.
                findEigVect2(A, roots[1], vec2, vec3);
                
                values(0) = roots[0];
                values(1) = roots[1];
                values(2) = roots[1];
                
                vectors.setCol(0, vec1);
                vectors.setCol(1, vec2);
                vectors.setCol(2, vec3);
                
            } else if (roots.size() == 3) {
                for (int i = 0; i < 3; ++i) {
                    findEigVect1(A, values(i), vec1);
                    vectors.setCol(i, vec1);
                }
                
            } else {
                assert(false);
            }
        }
    };
    
#else
    
    template <typename T>
    struct MatFloatingPointHelper {
        static void eigs(const Mat$D<T>& A, Vec$D<T>& values, Mat$D<T>& vectors) {
            LapackHelper$D<T>::eigs(A, values, vectors);
        }
    };
    
#endif
    
} // private namespace


    
//
// External functions
//
    
// Note: if you add other *free* functions (not member functions) here, make sure they are also added to 
// LINALG_INSTANTIATE_MATRIX macro at the end of this file.

template <typename T>
Mat$D<T> operator*(const Mat$D<T>& u, const Mat$D<T>& v) {
    Mat$D<T> w;
    T x;
    for (unsigned i = 0; i < $D; ++i) {
        for (unsigned j = 0; j < $D; ++j) {
            x = 0;
            for (unsigned k = 0; k < $D; ++k) {
                x += u(i,k) * v(k,j);
            }
            w(i,j) = x;
        }
    }

    return w;
}

template <typename T>
Vec$D<T> operator*(const Mat$D<T>& m, const Vec$D<T>& v) {
    Vec$D<T> w;
    double x;
    for (unsigned i = 0; i < $D; ++i) {
        x = 0;
        for (unsigned j = 0; j < $D; ++j) {
            x += m.m_data[j * $D + i] * v(j);
        }
        w(i) = T(x);
    }

    return w;
}

template <typename T>
std::ostream& operator<<(std::ostream& o, const Mat$D<T>& m) {
    for (unsigned i = 0; i < $D; ++i) {
        for (unsigned j = 0; j < $D; ++j) {
            o << m(i,j);
            if (j != $D - 1)
                o << ' ';
        }
        if (i != $D - 1)
            o << '\n';
    }
    return o;
}

template <typename T>
std::istream& operator>>(std::istream& in, Mat$D<T>& m) {
    for (unsigned i = 0; i < $D; ++i) 
        for (unsigned j = 0; j < $D; ++j) 
            in >> m(i,j);
    
    return in;
}


//
// Member functions
//

template <typename T>
Vec$D<T> Mat$D<T>::getCol(const unsigned c) const {
    const unsigned offset = c * $D;
    Vec$D<T> v;
    for (unsigned i = 0; i < $D; ++i)
        v(i) = m_data[offset + i];
    return v;
}

template <typename T>
void Mat$D<T>::setCol(const unsigned c, const Vec$D<T>& v) {
    const unsigned offset = c * $D;
    for (unsigned i = 0; i < $D; ++i)
        m_data[offset + i] = v(i);
}

template <typename T>
Vec$D<T> Mat$D<T>::getRow(const unsigned r) const {
    Vec$D<T> v;
    for (unsigned i = 0; i < $D; ++i)
        v(i) = m_data[i * $D + r];
    return v;
}

template <typename T>
void Mat$D<T>::setRow(const unsigned r, const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        m_data[i * $D + r] = v(i);
}

template <typename T>
Mat$D<T> Mat$D<T>::trans() const {
    Mat$D<T> w;
    for (unsigned i = 0; i < $D; ++i) {
        for (unsigned j = 0; j < $D; ++j) {
            w.m_data[j * $D + i] = m_data[i * $D + j];
        }
    }

    return w;
}


// 
// Static member functions
//

template <typename T>
const Mat$D<T>& Mat$D<T>::zero() {
    static bool needInit = true;
    static Mat$D<T> zeroMat;
    if (needInit) {
        T* d = zeroMat.data();
        for (unsigned i = 0; i < $D * $D; ++i)
            d[i] = 0;
        needInit = false;
    }
    return zeroMat;
}

template <typename T>
const Mat$D<T>& Mat$D<T>::ident() {
    static bool needInit = true;
    static Mat$D<T> identMat;
    if (needInit) {
        for (unsigned j = 0; j < $D; ++j)
            for (unsigned i = 0; i < $D; ++i)
                identMat(i,j) = (i == j);
        needInit = false;
    }
    return identMat;
}

template <typename T>
bool Mat$D<T>::equal(const Mat$D<T>& other, const T eps) const {
    const const_iterator iEnd = end();
    const_iterator i, j;
    for (i = begin(), j = other.begin(); i != iEnd; ++i, ++j) 
        if (!SignedHelper<T>::equal(*i, *j, eps))
            return false;
    assert(j == other.end());

    return true;
}

template <typename T>
Mat$D<T> Mat$D<T>::inverse(const Vec$D<T>& values, const Mat$D<T>& vectors) {
    // Invert the diagonals.
    Vec$D<T> d;
    for (unsigned i = 0; i < $D; ++i)
        d(i) = T(1.0) / values(i);

    // Form the matrix Q*D*Q^T, where Q is the eigenvector matrix, and
    // $D = diag(d).
    Mat$D<T> inv;
    for (unsigned j = 0; j < $D; ++j) {
        for (unsigned i = 0; i < $D; ++i) {
            inv(i,j) = 0;
            for (unsigned k = 0; k < $D; ++k) {
                inv(i,j) += vectors(i,k) * d(k) * vectors(j,k);
            }
        }
    }
    
    return inv;
}

template <typename T>
void Mat$D<T>::eigs(Vec$D<T>& values, Mat$D<T>& vectors) const {
#if $D == 2
    Mat2FloatingPointHelper<T>::eigs(*this, values, vectors);
#elif $D == 3 
    Mat3FloatingPointHelper<T>::eigs(*this, values, vectors);
#else
    MatFloatingPointHelper<T>::eigs(*this, values, vectors);
#endif
}

template <typename T>
std::string Mat$D<T>::description() {
    static std::string d;
    if (d.empty()) {
        std::ostringstream o;
        o << $D << "D " << type_description<T>() << " matrix";
        d = o.str();
    }

    return d;
}


/// Make the compiler generate code for the above template definitions so the 
/// linker will have something to link to.  The alternative is to place these
/// in the header, but then we have to pull in <ostream> and other headers 
/// that are big and unwieldy and then your compile times suck.
/// Only external functions that are not inlined are included here.
#define LINALG_INSTANTIATE_MATRIX(type) \
template class Mat$D<type>; \
template std::ostream& operator<<(std::ostream& o, const Mat$D<type>& v); \
template std::istream& operator>>(std::istream& i, Mat$D<type>& v); \
template Mat$D<type> operator*(const Mat$D<type>& u, const Mat$D<type>& v); \
template Vec$D<type> operator*(const Mat$D<type>& m, const Vec$D<type>& v); \
template Mat$D<type> min(const Mat$D<type>& u, const Mat$D<type>& v); \
template Mat$D<type> max(const Mat$D<type>& u, const Mat$D<type>& v);
    
#if !defined(LINALG_SKIP_DEFAULT_INSTANTIATIONS) && !defined(DOXYGEN_SHOULD_SKIP_THIS)
    LINALG_INSTANTIATE_MATRIX(float);
    LINALG_INSTANTIATE_MATRIX(double);
    LINALG_INSTANTIATE_MATRIX(int);
    LINALG_INSTANTIATE_MATRIX(unsigned);
#endif

#if defined(LINALG_INSTANTIATE_USER_TYPE)
    LINALG_INSTANTIATE_MATRIX(LINALG_INSTANTIATE_USER_TYPE);
#endif
    
#undef LINALG_INSTANTIATE_MATRIX


} // namespace
