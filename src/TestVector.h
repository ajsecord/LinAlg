#ifndef LINALG_TEST_VECTOR_H
#define LINALG_TEST_VECTOR_H

#include "Test.h"
#include <string>
#include <cstdlib>
#include <limits>

namespace LinAlg { namespace Test { namespace Vector {
    
    /// Run a series of tests on a vector.
    /// The extra template parameters are used to specialize tests for a particular type.
    template <typename Vec, 
            bool IntegerType = std::numeric_limits<typename Vec::value_type>::is_integer, 
            bool SignedType = std::numeric_limits<typename Vec::value_type>::is_signed>
    class Tester {
        public:
        typedef Random<Vec> Rand;
        typedef typename Vec::value_type T;
        static const unsigned D = Vec::num_dims;
        
#if 0
        /// Conversion test
        template <typename OtherT>
        static void testConvert(const T eps) {
            typedef LinAlg::Vector<OtherT,D> OtherVec;
            Rand v;
            OtherVec w = OtherVec::convert(v);
            for (unsigned i = 0; i < D; ++i) {
                if ((OtherT) v(i) != w(i))     // Need to allow for the generation of Inf, etc.
                    fail("Convert", (OtherT) v(i), w(i));            
            }
        }
#endif
        
        /// Generic tests that work for all types.
        static void test() {
            BasicTests<Vec>::test();
            const T eps = std::numeric_limits<T>::epsilon();

            std::string error;
            
            error = "Vector size";
            { 
                Vec v;
                if (v.size() != D || v.num_dims != D)
                    fail(error, v);
            }
            
            error = "Storage size";
            {
                Vec v;
                if (sizeof(v) != D * sizeof(T)) // This may not be valid -- compilers are allowed to add padding.
                    fail(error, v);
            }
            
#if 0
            error = "Constructor from std::vector";
            {
                std::vector<T> other(D);
                for (unsigned i = 0; i < D; ++i)
                    other[i] = (T)rand();
                
                Vec v(other);
                for (unsigned i = 0; i < D; ++i)
                    if (v(i) != other[i])
                        fail(error, v, other);
            }

            error = "Conversion from Vector<float,D>";
            {
                LinAlg::Vector<float,D> other;
                for (unsigned i = 0; i < D; ++i)
                    other[i] = (float) rand();
                
                Vec v = Vec::convert(other);
                for (unsigned i = 0; i < D; ++i)
                    if (v(i) != (T) other[i])
                        fail(error, v, other);
            }
            
            error = "Conversion from Vector<unsigned,D>";
            {
                LinAlg::Vector<unsigned,D> other;
                for (unsigned i = 0; i < D; ++i)
                    other[i] = (unsigned) rand();
                
                Vec v = Vec::convert(other);
                for (unsigned i = 0; i < D; ++i)
                    if (v(i) != (T) other[i])
                        fail(error, v, other);
            }
#endif
            
            error = "Indexing operators";
            {   
                Rand v;
                for (unsigned i = 0; i < D; ++i)
                    if (v(i) != v(i))
                        fail(error, v);
            }
            
            error = "Data access";
            { 
                Rand v;
                for (unsigned i = 0; i < D; ++i) 
                    if (v.data() + i != &v(i))
                        fail(error, v);
            }
            
            error = "Length";
            { 
                Rand v;
                T result = 0;
                for (unsigned i = 0; i < D; ++i)
                    result += v(i) * v(i);
                result = (T)sqrt((double)result);
                
                const T l = v.length();
                if (!SignedHelper<T>::equal(result, l, eps))
                    fail(error, v, result, l);
            }
            
            error = "Length squared";
            { 
                Rand v;
                T result = 0;
                for (unsigned i = 0; i < D; ++i)
                    result += v(i) * v(i);
                
                if (!SignedHelper<T>::equal(result, v.lengthSqr(), eps))
                    fail(error, v, result);
            }  
                                    
            error = "Inner product";
            {
                Rand v, w;
                T result = 0;
                for (unsigned i = 0; i < D; ++i)
                    result += v(i) * w(i);
                if (!SignedHelper<T>::equal(result, inner(v,w), eps))
                    fail(error, v, w);
            }
        
            error = "Minimum";
            {
                Rand v, w;
                Vec result = min(v, w);
                for (unsigned i = 0; i < D; ++i)
                    if (result(i) != (v(i) < w(i) ? v(i) : w(i)))
                        fail(error, v, w, result);
            }
            
            error = "Maximum";
            {
                Rand v, w;
                Vec result = max(v, w);
                for (unsigned i = 0; i < D; ++i)
                    if (result(i) != (v(i) > w(i) ? v(i) : w(i)))
                        fail(error, v, w, result);
            }
            
#if 0
            testConvert<float>(eps);
            testConvert<double>(eps);
            testConvert<int>(eps);
            testConvert<unsigned>(eps);
#endif
        }
    };

}}}
          


#endif
