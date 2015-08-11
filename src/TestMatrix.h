#ifndef LINALG_TEST_MATRIX_H
#define LINALG_TEST_MATRIX_H

#include "Test.h"
#include "TestVector.h"

#include <string>
#include <cstdlib>
#include <cmath>
#include <cassert>

namespace LinAlg { namespace Test { namespace Matrix {

    template <typename Mat>
    class RandomMatrix: public Mat {
        public:
        typedef Mat Base;
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;
        
        RandomMatrix() {
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j) 
                    (*this)(i,j) = (T)rand();
        }
        
    };

    template <typename Mat, 
        bool IntegerType = std::numeric_limits<typename Mat::value_type>::is_integer,
        bool SignedType = std::numeric_limits<typename Mat::value_type>::is_signed>
    struct ApproxEquality {
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;

        static void test() {
            const std::string error = "Approximate equality";
            Mat v;
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    v(i,j) = 1;

            if (!v.equal(v, 0))
                fail(error, v);
            
            // Get the smallest increment you can add to one and get a different value.
            // Does not work for integer types, will return zero.
            const T eps = std::numeric_limits<T>::epsilon();
            assert(eps > 0);
            Mat copy(v);
            
            // Add or subtract 2*epsilon from each component
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    copy(i,j) += 2 * (rand() % 2 ? eps : -eps);
            
            // The matrix should still be within two epsilons
            if (!v.equal(copy, 2 * eps))
                fail(error, v);
            
            // But not within one epsilon
            if (v.equal(copy, eps))
                fail(error, v);
        }
    };

    // Signed integers
    template <typename Mat>
    struct ApproxEquality<Mat, true, true> {
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;

        static void test() {
            const std::string error = "Approximate equality";
            Mat v;
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    v(i,j) = 1;

            if (!v.equal(v, 0))
                fail(error, v);
            
            // Get the smallest increment you can add to one and get a different value.
            const T eps = 1;
            Mat copy(v);
            
            // Add or subtract 2*epsilon from each component
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    copy(i,j) += 2 * (rand() % 2 ? eps : -eps);
            
            // The vector should still be within two epsilons
            if (!v.equal(copy, 2 * eps))
                fail(error, v);
            
            // But not within one epsilon
            if (v.equal(copy, eps))
                fail(error, v);
        }
    };

    // Unsigned integers
    template <typename Mat>
    struct ApproxEquality<Mat, true, false> {
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;
        
        static void test() {
            const std::string error = "Approximate equality";
            Mat v;
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    v(i,j) = 10;

            if (!v.equal(v, 0))
                fail(error, v);
            
            // Add or subtract 2*epsilon from each component
            Mat copy(v);
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    copy(i,j) += 2 * (rand() % 2 ? 1 : -1);
            
            // The vector should still be within two epsilons
            if (!v.equal(copy, 2))
                fail(error, v);
            
            // But not within one epsilon
            if (v.equal(copy, 1))
                fail(error, v);
        }
    };

    template <typename Mat,
        bool IntegerType = std::numeric_limits<typename Mat::value_type>::is_integer,
        bool SignedType = std::numeric_limits<typename Mat::value_type>::is_signed>
    struct UnaryMinus {
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;

        static void test() {
            const std::string& error = "Unary minus";
            RandomMatrix<Mat> v;
            Mat copy(v);
            copy = -copy;
            for (unsigned i = 0; i < D; ++i)
                for (unsigned j = 0; j < D; ++j)
                    if (-v(i,j) != copy(i,j))
                        fail(error, copy);
        }
    };

    // Unsigned quantities -- skip test
    template <typename Mat, bool IntegerType>
    struct UnaryMinus<Mat, IntegerType, false> {
        static void test() {}
    };

    /// Run a series of tests on a matrix.
    template <typename Mat>
    class Tester {
        public:
        typedef typename Mat::value_type T;
        static const unsigned D = Mat::num_dims;

        typedef typename Mat::vec_type Vec;
        typedef RandomMatrix<Mat> Rand;
        
#if 0
        /// Conversion test
        template <typename OtherT>
            static void testConvert(const T eps) {
                typedef LinAlg::Matrix<OtherT,D> OtherMat;
                Rand v;
                OtherMat w = OtherMat::convert(v);
                for (unsigned j = 0; j < D; ++j) {
                    for (unsigned i = 0; i < D; ++i) {
                        //if (!SignedHelper<T>::equal((OtherT) v(i), w(i), eps))
                        if ((OtherT) v(i,j) != w(i,j))     // Need to allow for the generation of Inf, etc.
                            fail("Convert", (OtherT) v(i,j), w(i,j));        
                    }
                }
            }
#endif
        
        static void test() {
            //const T eps = std::numeric_limits<T>::epsilon();

            BasicTests<Mat>::test();
          
            std::string error;
            
            error = "Matrix size";
            { 
                Mat v;
                if (v.size1() != D || v.size2() != D || v.num_dims != D)
                    fail(error, v);
            }
            
            error = "Storage size";
            {
                Mat v;
                if (sizeof(v) != D * D * sizeof(T)) // This may not be valid -- compilers are allowed to add padding.
                    fail(error, v);
            }
            
            error = "Identity matrix";
            {
                const Mat& v = Mat::ident();
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j) 
                        if (v(i,j) != T(i == j))
                            fail(error, v);
            }
            
#if 0
            error = "Constructor from std::vector";
            {
                std::vector<T> other(D*D);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        other[i + D * j] = (T)rand();
                
                Mat v(other);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (v(i,j) != other[i + D * j])
                            fail(error, v, other);
            }
            
            error = "Constructor from Matrix<float,D>";
            {
                LinAlg::Matrix<float,D> other;
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        other(i,j) = (float) rand();
                
                Mat v = Mat::convert(other);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (v(i,j) != (T) other(i,j))
                            fail(error, v, other);
            }
            
            error = "Constructor from Matrix<unsigned,D>";
            {
                LinAlg::Matrix<unsigned,D> other;
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        other(i,j) = (unsigned) rand();
                
                Mat v = Mat::convert(other);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (v(i,j) != (T) other(i,j))
                            fail(error, v, other);
            }
#endif
            
            error = "Data access";
            { 
                Rand v;
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j) 
                        if (v.data() + (j * D + i) != &v(i,j))
                            fail(error, v);
            }
            
            error = "Transpose";
            {
                Rand v;
                Mat copy(v);
                copy = copy.trans();
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (v(i,j) != copy(j,i))
                            fail(error, v, copy);
            }
            
            error = "Column access";
            {
                Rand mat;
                Vec v;
                for (unsigned j = 0; j < D; ++j) {
                    v = mat.getCol(j);
                    for (unsigned i = 0; i < D; ++i) 
                        if (v(i) != mat(i,j))
                            fail(error, mat, j, v);
                }

                for (unsigned j = 0; j < D; ++j) {
                    Random<Vec> randVec;
                    mat.setCol(j, randVec);
                    for (unsigned i = 0; i < D; ++i) 
                        if (mat(i,j) != randVec(i))
                            fail(error, mat, j, v);
                }
            }
            
            error = "Row access";
            {
                Rand mat;
                Vec v;
                for (unsigned i = 0; i < D; ++i) {
                    v = mat.getRow(i);
                    for (unsigned j = 0; j < D; ++j) 
                        if (v(j) != mat(i,j))
                            fail(error, mat, i, v);
                }

                for (unsigned i = 0; i < D; ++i) {
                    Random<Vec> randVec;
                    mat.setRow(i, randVec);
                    for (unsigned j = 0; j < D; ++j) 
                        if (mat(i,j) != randVec(j))
                            fail(error, mat, i, v);
                }
            }
            
            error = "Minimum";
            {
                Rand v, w;
                Mat result = min(v, w);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (result(i,j) != (v(i,j) < w(i,j) ? v(i,j) : w(i,j)))
                            fail(error, v, w, result);
            }
            
            error = "Maximum";
            {
                Rand v, w;
                Mat result = max(v, w);
                for (unsigned i = 0; i < D; ++i)
                    for (unsigned j = 0; j < D; ++j)
                        if (result(i,j) != (v(i,j) > w(i,j) ? v(i,j) : w(i,j)))
                            fail(error, v, w, result);
            }
            
            
    #if 0
            // Not for integer types
            error = "Eigenvalue decomposition";
            {
                assert(eps > 0);
                Rand vals;
                
                // Construct an orthogonal basis
                Rand mat;
                
                Vec c = mat.getCol(0);
                c /= c.length();
                mat.setCol(0, c);

                for (unsigned i = 1; i < D; ++i) {
                    Vec c = mat.getCol(i);
                    for (unsigned j = 0; j < i; ++i) {
                        c -= proj(c, mat.getCol(j));
                        assert(std::abs((double)inner(c, mat.getCol(j))) < eps);
                    }
                    assert(c.lengthSqr() > 0);
                    mat.setCol(i, c / c.length());
                }
                
                for (unsigned i = 0; i < D; ++i)
                    mat.setCol(i, mat.getCol(i) * vals(i));
                
                Vec compVals;
                Mat compVecs;
                mat.eigs(compVals, compVecs);
            }
    #endif
            
            // Type-specific tests
            ApproxEquality<Mat>::test();
            UnaryMinus<Mat>::test();

#if 0
            testConvert<double>(eps);
            testConvert<float>(eps);
            testConvert<int>(eps);
            testConvert<unsigned>(eps);
#endif       
        }
    };
}}} // namespaces


#endif
