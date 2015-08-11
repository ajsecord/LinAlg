/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * Vec$D.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_VEC$D_H
#define LINALG_VEC$D_H

#if defined(min) || defined(max)
#   error min/max macros defined -- probably contamination from scabby old windows.h.
#endif

#include <iosfwd>
#include <string>

#include <LinAlg/Vec$DFwd.h>

/// Linear algebra types and operations
namespace LinAlg {

/// Negate a vector.
template <typename T>
    Vec$D<T> operator-(const Vec$D<T>& u);

/// Add two vectors.
template <typename T> 
    Vec$D<T> operator+(const Vec$D<T>& u, const Vec$D<T>& v);

/// Subtract two vectors.
template <typename T> 
    Vec$D<T> operator-(const Vec$D<T>& u, const Vec$D<T>& v);
    
/// Right-multiply a vector by a scalar.
template <typename T>
    Vec$D<T> operator*(const Vec$D<T>& u, const T s);

/// Left-multiply a vector by a scalar.
template <typename T>
    Vec$D<T> operator*(const T s, const Vec$D<T>& u);

/// Right-divide a vector by a scalar.
template <typename T>
    Vec$D<T> operator/(const Vec$D<T>& u, const T s);
    
/// Exact equality -- use with caution on floating-point base types.
template <typename T> 
    bool operator==(const Vec$D<T>& u, const Vec$D<T>& v);

/// Exact inequality -- use with caution on floating-point base types.
template <typename T> 
    bool operator!=(const Vec$D<T>& u, const Vec$D<T>& v);

/// The inner (dot) product of \c u and \c v.
template <typename T> 
    T inner(const Vec$D<T>& u, const Vec$D<T>& v);

/// Return the projection of \c u onto \c v.
template <typename T>
    Vec$D<T> proj(const Vec$D<T>& u, const Vec$D<T>& v);

/// Per-element minimum.
template <typename T>
    Vec$D<T> min(const Vec$D<T>& u, const Vec$D<T>& v);

/// Per-element maximum.
template <typename T>
    Vec$D<T> max(const Vec$D<T>& u, const Vec$D<T>& v);

/// Unit vector in direction of \c u.  Does not check for u.length() == 0.
template <typename T>
    Vec$D<T> normalized(const Vec$D<T>& u);


/// \class Vec$D Vec$D.h <LinAlg/Vec$D.h>
/// A linear algebra vector class of fixed type and size $D.
template <typename T>
class Vec$D {
    public:
    /// This type
    typedef Vec$D<T> self_type;

    /// The underlying m_data type
    typedef T value_type;
    
    /// The number of dimensions
    enum { num_dims = $D };  

    /// The number of elements
    enum { num_elements = $D };
    
    /// An iterator
    typedef T* iterator;
    
    /// A constant iterator.
    typedef const T* const_iterator;
    
    /// Create an uninitialized vector
    Vec$D();

    /// Create a vector from arguments
    Vec$D($VECTOR_CONSTRUCTOR_ARGUMENTS);
    
    /// Copy a vector of another type.
    /// \note If the source and destination types are not the same, then the
    /// elements of \c v are cast into the appropriate type.
    template <typename Other>
    Vec$D(const Vec$D<Other>& v);

    /// Assign another vector to this one.
    Vec$D& operator=(const Vec$D& v);

    /// Access an element read-only.
    const T& operator()(const unsigned i) const;
    
    /// Access an element read-write.
    T& operator()(const unsigned i);

    /// Get a read-only pointer to the first element.
    const T* data() const;
    
    /// Get a read-write pointer to the first element.
    T* data();
    
    /// Get a read-only iterator to the first element.
    const_iterator begin() const;
    
    /// Get a read-write pointer to the first element.
    iterator begin();
    
    /// Get a read-only iterator to past the last element.
    const_iterator end() const;
    
    /// Get a read-write iterator to past the last element.
    iterator end();
    
    /// Return the number of elements in the vector.
    unsigned size() const;

    /// Limit each element of the vector to the range [l.first, l.second]
    void limit(const std::pair<T,T>& l);

    /// Limit each element of the vector to the range [lower, upper]
    void limit(const T lower, const T upper);
    
    /// Limit element i of the vector to the range [lower[i], upper[i]]
    void limit(const Vec$D& lower, const Vec$D& upper);

    /// Return the sum of all the elements.
    T sum() const;
    
    /// Return the product of all the elements.
    T prod() const;
    
    /// Return the geometric length of the vector.
    T length() const;
    
    /// Return the geometric length squared of the vector.
    T lengthSqr() const;

    /// Is this vector the same as another, within epsilon?
    bool equal(const Vec$D& v, const T eps) const;
   
    /// Add another vector to this one.
    Vec$D& operator+=(const Vec$D& v);
    
    /// Subtract another vector from this one.
    Vec$D& operator-=(const Vec$D& v);

    /// Multiply this vector by a scalar.
    template <typename Scalar>
    Vec$D& operator*=(const Scalar d);
    
    /// Divide this vector by a scalar.  Does not check for zero.
    template <typename Scalar>
    Vec$D& operator/=(const Scalar d);
    
    // Dear god the things you have to do to get a templated friend function
    friend Vec$D operator- <>(const Vec$D& u);
    friend Vec$D operator+ <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D operator- <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D operator* <>(const Vec$D<T>& u, const T s);
    friend Vec$D operator* <>(const T s, const Vec$D<T>& u);
    friend Vec$D operator/ <>(const Vec$D<T>& u, const T s);
    friend bool operator== <>(const Vec$D& u, const Vec$D& v);
    friend bool operator!= <>(const Vec$D& u, const Vec$D& v);
    friend T inner <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D proj <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D min <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D max <>(const Vec$D& u, const Vec$D& v);
    friend Vec$D normalized <>(const Vec$D& u);
    
    /// Return the zero vector.
    static const Vec$D& zero();
    
    /// Return the unit vector with v(i) == 0 and v(index) == 1.
    static const Vec$D unit(const unsigned index);

    /// Return a description like "2D float vector", useful for debugging and messages.
    static std::string description();
    
    protected:
    T m_data[$D];                  ///< The element storage
};

    
//
// Global function declarations
//

/// Print a vector to the output stream.
template <typename T>
std::ostream& operator<<(std::ostream& o, const Vec$D<T>& v);
    
/// Read a vector from the input stream.
template <typename T>
std::istream& operator>>(std::istream& i, Vec$D<T>& v);

} // namespace


//
// Inline function definitions
//
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace LinAlg {

template <typename T>
inline
Vec$D<T>::Vec$D() {
#if defined(LINALG_DEFAULT_INITIALIZE_TO_ZERO)
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] = 0;
#endif
}
    
template <typename T>
inline
Vec$D<T>::Vec$D($VECTOR_CONSTRUCTOR_ARGUMENTS) {
    $VECTOR_CONSTRUCTOR_BODY
}

template <typename T>
template <typename Other>
inline
Vec$D<T>::Vec$D(const Vec$D<Other>& v) {
    for (int i = 0; i < $D; ++i)
        m_data[i] = (T) v(i);
}

template <typename T>
inline
Vec$D<T>& Vec$D<T>::operator=(const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] = v.m_data[i];
    return *this;
}

template <typename T>
inline
const T& Vec$D<T>::operator()(const unsigned i) const {
    assert(i < $D);
    return m_data[i];
}

template <typename T>
inline
T& Vec$D<T>::operator()(const unsigned i) {
    assert(i < $D);
    return m_data[i];
}

template <typename T>
inline
const T* Vec$D<T>::data() const {
    return m_data;
}

template <typename T>
inline
T* Vec$D<T>::data() {
    return m_data;
}

template <typename T>
inline 
typename Vec$D<T>::const_iterator Vec$D<T>::begin() const {
    return m_data;
}

template <typename T>
inline 
typename Vec$D<T>::iterator Vec$D<T>::begin() {
    return m_data;
}

template <typename T>
inline 
typename Vec$D<T>::const_iterator Vec$D<T>::end() const {
    return m_data + $D;
}

template <typename T>
inline 
typename Vec$D<T>::iterator Vec$D<T>::end() {
    return m_data + $D;
}

template <typename T>
inline
unsigned Vec$D<T>::size() const {
    return $D;
}

template <typename T>
inline
void Vec$D<T>::limit(const std::pair<T,T>& l) {
    for (unsigned i = 0; i < $D; ++i) {
        if (m_data[i] < l.first) 
            m_data[i] = l.first;
        if (m_data[i] > l.second) 
            m_data[i] = l.second;
    }
}

template <typename T>
inline
void Vec$D<T>::limit(const T lower, const T upper) {
    for (unsigned i = 0; i < $D; ++i) {
        if (m_data[i] < lower) 
            m_data[i] = lower;
        if (m_data[i] > upper)
            m_data[i] = upper;
    }
}

template <typename T>
inline
void Vec$D<T>::limit(const Vec$D<T>& lower, const Vec$D<T>& upper) {
    for (unsigned i = 0; i < $D; ++i) {
        if (m_data[i] < lower(i)) 
            m_data[i] = lower(i);
        if (m_data[i] > upper(i))
            m_data[i] = upper(i);
    }
}

template <typename T>
inline
Vec$D<T>& Vec$D<T>::operator+=(const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] += v.m_data[i];
    return *this;
}

template <typename T>
inline
Vec$D<T>& Vec$D<T>::operator-=(const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] -= v.m_data[i];
    return *this;
}

template <typename T>
template <typename Scalar>
inline
Vec$D<T>& Vec$D<T>::operator*=(const Scalar d) {
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] *= d;
    return *this;
}

template <typename T>
template <typename Scalar>
inline
Vec$D<T>& Vec$D<T>::operator/=(const Scalar d) {
    // Don't get too tricky -- think of T == int and d = 2.
    for (unsigned i = 0; i < $D; ++i)
        m_data[i] /= d;
    return *this;
}

template <typename T>
inline 
T Vec$D<T>::sum() const {
    T w = 0;
    for (unsigned i = 0; i < $D; ++i)
        w += m_data[i];
    return w;
}

template <typename T>
inline 
T Vec$D<T>::prod() const {
    T w = 1;
    for (unsigned i = 0; i < $D; ++i)
        w *= m_data[i];
    return w;
}

template <typename T>
inline 
T Vec$D<T>::length() const {
    T w = 0;
    for (unsigned i = 0; i < $D; ++i)
        w += m_data[i] * m_data[i];
    return (T)sqrt((double)w);      // Need this to resolve ambiguity when T=int
}

template <typename T>
inline 
T Vec$D<T>::lengthSqr() const {
    T w = 0;
    for (unsigned i = 0; i < $D; ++i)
        w += m_data[i] * m_data[i];
    return w;
}

template <typename T>
inline
Vec$D<T> operator-(const Vec$D<T>& u) {
    Vec$D<T> v;
    for (unsigned i = 0; i < $D; ++i)
        v.m_data[i] = -u.m_data[i];
    return v;
}

template <typename T>
inline
Vec$D<T> operator+(const Vec$D<T>& u, const Vec$D<T>& v) {
    Vec$D<T> w;
    for (unsigned i = 0; i < $D; ++i) 
        w.m_data[i] = u.m_data[i] + v.m_data[i];
    return w;
}

template <typename T>
inline
Vec$D<T> operator-(const Vec$D<T>& u, const Vec$D<T>& v) {
    Vec$D<T> w;
    for (unsigned i = 0; i < $D; ++i) 
        w.m_data[i] = u.m_data[i] - v.m_data[i];
    return w;
}
    
template <typename T>
inline
Vec$D<T> operator*(const Vec$D<T>& u, const T s) {
    Vec$D<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Vec$D<T> operator*(const T s, const Vec$D<T>& u) {
    Vec$D<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Vec$D<T> operator/(const Vec$D<T>& u, const T s) {
    Vec$D<T> w(u);
    w /= s;
    return w;
}
    
template <typename T>
inline
bool operator==(const Vec$D<T>& u, const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        if (u.m_data[i] != v.m_data[i])
            return false;
    return true;
}

template <typename T>
inline
bool operator!=(const Vec$D<T>& u, const Vec$D<T>& v) {
    for (unsigned i = 0; i < $D; ++i)
        if (u.m_data[i] != v.m_data[i])
            return true;
    return false;
}

template <typename T>
inline
T inner(const Vec$D<T>& u, const Vec$D<T>& v) {
    T w = 0;
    for (unsigned i = 0; i < $D; ++i)
        w += u.m_data[i] * v.m_data[i];
    return w;
}

template <typename T>
inline
Vec$D<T> proj(const Vec$D<T>& u, const Vec$D<T>& v) {
    return v * inner(u, v) / v.lengthSqr();
}

template <typename T>
inline
Vec$D<T> min(const Vec$D<T>& u, const Vec$D<T>& v) {
    Vec$D<T> result;
    for (unsigned i = 0; i < $D; ++i)
        result.m_data[i] = (u.m_data[i] < v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Vec$D<T> max(const Vec$D<T>& u, const Vec$D<T>& v) {
    Vec$D<T> result;
    for (unsigned i = 0; i < $D; ++i)
        result.m_data[i] = (u.m_data[i] > v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Vec$D<T> normalized(const Vec$D<T>& u) {
    return u / u.length();
}



} // namespace

#endif
