#ifndef LINALG_TEST_H
#define LINALG_TEST_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <cassert>

#include "LinAlg/TypeHelpers.h"  // Just grabbing the SignedHelper::equal() method for convenience.

namespace LinAlg { namespace Test {

template <typename T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
    o << '(';
    for (unsigned i = 0; i < v.size(); ++i)
        o << v[i] << (i < v.size() - 1 ? ", " : "");
    o << ')';
    return o;
}

/// Throw an exception because a test failed with some vector.
template <typename Vec>
static void fail(const std::string& s, const Vec& v) {
    std::ostringstream str;
    str.precision(16);
    str << s << ". Value: " << v;
    throw std::runtime_error(str.str());
}

/// Throw an exception because a test failed with some vector and a 
/// second thing.
template <typename Vec, typename Other>
static void fail(const std::string& s, const Vec& v, const Other& w) {
    std::ostringstream str;
    str.precision(16);
    str << s << ". Value: " << v << ", value: " << w;
    throw std::runtime_error(str.str());
}


/// Throw an exception because a test failed with some vector and a 
/// second thing.
template <typename Vec, typename Other1, typename Other2>
static void fail(const std::string& s, const Vec& v, const Other1& w, const Other2& x) {
    std::ostringstream str;
    str.precision(16);
    str << s << ". Value: " << v << ", value: " << w << ", value: " << x;
    throw std::runtime_error(str.str());
}

/// Throw an exception because a test failed.
template <typename Vec, typename Other1, typename Other2, typename Other3>
static void fail(const std::string& s, const Vec& v, const Other1& w, const Other2& x, const Other3& y) {
    std::ostringstream str;
    str.precision(16);
    str << s << ". Value: " << v << ", value: " << w << ", value: " << x << ", value: " << y;
    throw std::runtime_error(str.str());
}

template <typename Container>
class Random: public Container {
public:
    typedef typename Container::iterator Iter;
    
    Random() {
        // Keep the size of the elements somewhat reasonable to avoid 
        // overflowing to infinity in cases like taking the product of all elements.
        typedef typename Container::value_type T;
        const double maxValue = pow((double) std::numeric_limits<T>::max(), 1.0 / Container::num_elements);
        for (Iter i = Container::begin(); i != Container::end(); ++i)
            *i = (T)(maxValue * (rand() / ((double)RAND_MAX - 1) - 0.5));
    }
};



template <typename Container, 
          bool IntegerType = std::numeric_limits<typename Container::value_type>::is_integer, 
          bool SignedType = std::numeric_limits<typename Container::value_type>::is_signed>
struct ApproxEquality {
    static void test() {
        typedef typename Container::iterator Iter;
        typedef typename Container::value_type T;
        const std::string error = "Approximate equality";
        Container v;
        
        for (Iter i = v.begin(); i != v.end(); ++i)
            *i = 1;
        
        if (!v.equal(v, 0))
            fail(error, v);
        
        // Get the smallest increment you can add to one and get a different value.
        // Does not work for integer types, will return zero.
        const T eps = std::numeric_limits<T>::epsilon();
        assert(eps > 0);
        
        Container copy(v);
        // Add or subtract 2*epsilon from each component
        for (Iter i = copy.begin(); i != copy.end(); ++i)
            *i += 2 * (rand() % 2 ? eps : -eps);
        
        // The vector should still be within two epsilons
        if (!v.equal(copy, 2 * eps))
            fail(error, v);
        
        // But not within one epsilon
        if (v.equal(copy, eps))
            fail(error, v);
    }
};

// Signed integers
template<typename Container>
struct ApproxEquality<Container, true, true> {
    static void test() {
        typedef typename Container::iterator Iter;
        typedef typename Container::value_type T;
        const std::string error = "Approximate equality";
        
        Container v;
        for (Iter i = v.begin(); i != v.end(); ++i)
            *i = 1;
        
        if (!v.equal(v, 0))
            fail(error, v);
        
        // Get the smallest increment you can add to one and get a different value.
        const T eps = 1;
        Container copy(v);
        
        // Add or subtract 2*epsilon from each component
        for (Iter i = copy.begin(); i != copy.end(); ++i)
            *i += 2 * (rand() % 2 ? eps : -eps);
        
        // The vector should still be within two epsilons
        if (!v.equal(copy, 2 * eps))
            fail(error, v);
        
        // But not within one epsilon
        if (v.equal(copy, eps))
            fail(error, v);
    }
};

// Unsigned integers
template <typename Container>
struct ApproxEquality<Container, true, false> {
    static void test() {
        typedef typename Container::iterator Iter;
        typedef typename Container::value_type T;
        const std::string error = "Approximate equality";
        
        Container v;
        for (Iter i = v.begin(); i != v.end(); ++i)
            *i = 10;
        
        if (!v.equal(v, 0))
            fail(error, v);
        
        // Add or subtract 2*epsilon from each component
        Container copy(v);
        for (Iter i = copy.begin(); i != copy.end(); ++i)
            *i += 2 * (rand() % 2 ? 1 : -1);
        
        // The vector should still be within two epsilons
        if (!v.equal(copy, 2))
            fail(error, v);
        
        // But not within one epsilon
        if (v.equal(copy, 1))
            fail(error, v);
    }
};

// Generic version
template <typename Container, 
          bool SignedType = std::numeric_limits<typename Container::value_type>::is_signed>
struct UnaryMinus {
    static void test() {
        typedef typename Container::iterator Iter;
        typedef typename Container::value_type T;
        const std::string error = "Unary minus";
        Random<Container> v;
        Container copy(v);
        copy = -copy;
        for (Iter i = v.begin(), j = copy.begin(); i != v.end(); ++i, ++j) 
            if (-(*i) != *j)
                fail(error, copy);
    }
};

// Unsigned quantities -- skip test
template <typename Container>
struct UnaryMinus<Container, false> {
    static void test() {}
};

/// Tests that work for any container, ie. element-only tests.
template <typename Container>
struct BasicTests {
    static void test() {
        typedef typename Container::const_iterator CIter;
        typedef typename Container::iterator Iter;
        typedef typename Container::value_type T;
        typedef Random<Container> Rand;
        const T eps = std::numeric_limits<T>::epsilon();

        std::string error;
        
        error = "Zero vector/matrix";
        {
            const Container& v = Container::zero();
            for (CIter i = v.begin(); i != v.end(); ++i)
                if (*i != 0)
                    fail(error, v);
        }
        
        error = "Copy constructor";
        {
            Rand v;
            Container c(v);
            for (Iter i = v.begin(), j = c.begin(); i != v.end(); ++i, ++j)
                if (*i != *j)
                    fail(error, v, c);
        }
        
        error = "Assignment operator";
        {  
            Rand v;
            Container c;
            c = v;
            for (Iter i = v.begin(), j = c.begin(); i != v.end(); ++i, ++j)
                if (*i != *j)
                    fail(error, v, c);
        }
        
        error = "Component limits";
        {
            Container lower, upper;
            for (Iter i = lower.begin(), j = upper.begin(); i != lower.end(); ++i, ++j) {
                *i = (T)rand() / 2;
                *j = (T)rand() / 2 + RAND_MAX / 2;
            }
            
            // Using the first elements of lower, upper as a convenient single set of limits
            { 
                Rand v;
                v.limit(*lower.begin(), *upper.begin());
                for (Iter i = v.begin(); i != v.end(); ++i) 
                    if (*i < *lower.begin() || *i > *upper.begin())
                        fail(error, v);
            }
            
            // Using lower[0], upper[0] as a convenient single set of limits
            { 
                Rand v;
                v.limit(std::make_pair(*lower.begin(), *upper.begin()));
                for (Iter i = v.begin(); i != v.end(); ++i) 
                    if (*i < *lower.begin() || *i > *upper.begin())
                        fail(error, v);
            }
            
            { 
                Rand v;
                v.limit(lower, upper);
                for (Iter i = v.begin(), l = lower.begin(), u = upper.begin(); i != v.end(); ++i, ++l, ++u) 
                    if (*i < *l || *i > *u)
                        fail(error, v);
            }
        }
        
        error = "Sum";
        { 
            Rand v;
            T result = 0;
            for (Iter i = v.begin(); i != v.end(); ++i) 
                result += *i;
            
            if (!SignedHelper<T>::equal(result, v.sum(), eps))
                fail(error, v, result);
        }  
        
        error = "Product";
        { 
            Rand v;
            T result = 1;
            for (Iter i = v.begin(); i != v.end(); ++i) 
                result *= *i;
            
            if (!SignedHelper<T>::equal(result, v.prod(), eps))
                fail(error, v, v.prod(), result);
        }  
        
        error = "Operator +=";
        {
            Rand v, w;
            Container copy(v);
            v += w;
            for (Iter vi = v.begin(), wi = w.begin(), ci = copy.begin(); vi != v.end(); ++vi, ++wi, ++ci)
                if (!SignedHelper<T>::equal(*vi, *ci + *wi, eps))
                    fail(error, v);
        }
        
        error = "Operator -=";
        {
            Rand v, w;
            Container copy(v);
            v -= w;
            for (Iter vi = v.begin(), wi = w.begin(), ci = copy.begin(); vi != v.end(); ++vi, ++wi, ++ci)
                if (!SignedHelper<T>::equal(*vi, *ci - *wi, eps))
                    fail(error, v);
        }
        
        error = "Operator *=";
        {
            Rand v;
            T s = (T)rand();
            Container copy(v);
            v *= s;
            for (Iter vi = v.begin(), ci = copy.begin(); vi != v.end(); ++vi, ++ci)
                if (!SignedHelper<T>::equal(*vi, *ci * s, eps)) 
                    fail(error, v, copy, s, (*vi - *ci * s));
        }
        
        error = "Operator /=";
        {
            Rand v;
            T s;
            for (;;) {
                s = (T)rand();
                if (s != 0)
                    break;
            }
            Container copy(v);
            v /= s;
            for (Iter vi = v.begin(), ci = copy.begin(); vi != v.end(); ++vi, ++ci)
                if (!SignedHelper<T>::equal(*vi, *ci / s, eps)) 
                    fail(error, v, copy, s);
        }
        
        error = "Binary addition";
        {
            Rand v, w;
            Container result = v + w;
            for (Iter vi = v.begin(), wi = w.begin(), ri = result.begin(); vi != v.end(); ++vi, ++wi, ++ri)
                if (!SignedHelper<T>::equal(*ri, *vi + *wi, eps))
                    fail(error, v);
        }
        
        error = "Binary subtraction";
        {
            Rand v, w;
            Container result = v - w;
            for (Iter vi = v.begin(), wi = w.begin(), ri = result.begin(); vi != v.end(); ++vi, ++wi, ++ri)
                if (!SignedHelper<T>::equal(*ri, *vi - *wi, eps))
                    fail(error, v);
        }
        
        error = "Scalar multiplication";
        {
            Rand v;
            T s = (T)rand();
            Container result = v * s;
            for (Iter vi = v.begin(), ri = result.begin(); vi != v.end(); ++vi, ++ri)
                if (!SignedHelper<T>::equal(*ri, *vi * s, eps))
                    fail(error, result);
            
            result = s * v;
            for (Iter vi = v.begin(), ri = result.begin(); vi != v.end(); ++vi, ++ri)
                if (!SignedHelper<T>::equal(*ri, *vi * s, eps))
                    fail(error, result);
        }
        
        error = "Scalar division";
        {
            Rand v;

            T s = 0;
            while (s == 0)
                s = (T)rand();
    
            Container result = v / s;
            for (Iter vi = v.begin(), ri = result.begin(); vi != v.end(); ++vi, ++ri)
                if (!SignedHelper<T>::equal(*ri, *vi / s, eps))
                    fail(error, result);
        }
        
        error = "Maximum";
        {
            Rand v, w;
            Container result = max(v, w);
            for (Iter vi = v.begin(), wi = w.begin(), ri = result.begin(); vi != v.end(); ++vi, ++wi, ++ri)
                if (!SignedHelper<T>::equal(*ri, std::max(*vi, *wi), eps))
                    fail(error, v);
        }
        
        error = "Minimum";
        {
            Rand v, w;
            Container result = min(v, w);
            for (Iter vi = v.begin(), wi = w.begin(), ri = result.begin(); vi != v.end(); ++vi, ++wi, ++ri)
                if (!SignedHelper<T>::equal(*ri, std::min(*vi, *wi), eps))
                    fail(error, v);
        }
        
        
        error = "Equality";
        {
            Rand v;
            if (!(v == v))
                fail(error, v);
        }
        
        error = "Inequality";
        {
            Rand v;
            if (v != v)
                fail(error, v);
        }    
        
        // Type-specific tests
        ApproxEquality<Container>::test();
        UnaryMinus<Container>::test();
    }
};

}} // namespaces




#endif