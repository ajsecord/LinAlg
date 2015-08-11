// Warning: this file was generated from src/LinAlg/VecTemplate.h.  Do not edit if you value your changes!

/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * Vec4.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_VEC4_H
#define LINALG_VEC4_H

#if defined(min) || defined(max)
#   error min/max macros defined -- probably contamination from scabby old windows.h.
#endif

#include <iosfwd>
#include <string>

#include <LinAlg/Vec4Fwd.h>

/// Linear algebra types and operations
namespace LinAlg {

/// Negate a vector.
template <typename T>
    Vec4<T> operator-(const Vec4<T>& u);

/// Add two vectors.
template <typename T> 
    Vec4<T> operator+(const Vec4<T>& u, const Vec4<T>& v);

/// Subtract two vectors.
template <typename T> 
    Vec4<T> operator-(const Vec4<T>& u, const Vec4<T>& v);
    
/// Right-multiply a vector by a scalar.
template <typename T>
    Vec4<T> operator*(const Vec4<T>& u, const T s);

/// Left-multiply a vector by a scalar.
template <typename T>
    Vec4<T> operator*(const T s, const Vec4<T>& u);

/// Right-divide a vector by a scalar.
template <typename T>
    Vec4<T> operator/(const Vec4<T>& u, const T s);
    
/// Exact equality -- use with caution on floating-point base types.
template <typename T> 
    bool operator==(const Vec4<T>& u, const Vec4<T>& v);

/// Exact inequality -- use with caution on floating-point base types.
template <typename T> 
    bool operator!=(const Vec4<T>& u, const Vec4<T>& v);

/// The inner (dot) product of \c u and \c v.
template <typename T> 
    T inner(const Vec4<T>& u, const Vec4<T>& v);

/// Return the projection of \c u onto \c v.
template <typename T>
    Vec4<T> proj(const Vec4<T>& u, const Vec4<T>& v);

/// Per-element minimum.
template <typename T>
    Vec4<T> min(const Vec4<T>& u, const Vec4<T>& v);

/// Per-element maximum.
template <typename T>
    Vec4<T> max(const Vec4<T>& u, const Vec4<T>& v);

/// Unit vector in direction of \c u.  Does not check for u.length() == 0.
template <typename T>
    Vec4<T> normalized(const Vec4<T>& u);


/// \class Vec4 Vec4.h <LinAlg/Vec4.h>
/// A linear algebra vector class of fixed type and size 4.
template <typename T>
class Vec4 {
    public:
    /// This type
    typedef Vec4<T> self_type;

    /// The underlying m_data type
    typedef T value_type;
    
    /// The number of dimensions
    enum { num_dims = 4 };  

    /// The number of elements
    enum { num_elements = 4 };
    
    /// An iterator
    typedef T* iterator;
    
    /// A constant iterator.
    typedef const T* const_iterator;
    
    /// Create an uninitialized vector
    Vec4();

    /// Create a vector from arguments
    Vec4(const T v0, const T v1, const T v2, const T v3);
    
    /// Copy a vector of another type.
    /// \note If the source and destination types are not the same, then the
    /// elements of \c v are cast into the appropriate type.
    template <typename Other>
    Vec4(const Vec4<Other>& v);

    /// Assign another vector to this one.
    Vec4& operator=(const Vec4& v);

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
    void limit(const Vec4& lower, const Vec4& upper);

    /// Return the sum of all the elements.
    T sum() const;
    
    /// Return the product of all the elements.
    T prod() const;
    
    /// Return the geometric length of the vector.
    T length() const;
    
    /// Return the geometric length squared of the vector.
    T lengthSqr() const;

    /// Is this vector the same as another, within epsilon?
    bool equal(const Vec4& v, const T eps) const;
   
    /// Add another vector to this one.
    Vec4& operator+=(const Vec4& v);
    
    /// Subtract another vector from this one.
    Vec4& operator-=(const Vec4& v);

    /// Multiply this vector by a scalar.
    template <typename Scalar>
    Vec4& operator*=(const Scalar d);
    
    /// Divide this vector by a scalar.  Does not check for zero.
    template <typename Scalar>
    Vec4& operator/=(const Scalar d);
    
    // Dear god the things you have to do to get a templated friend function
    friend Vec4 operator- <>(const Vec4& u);
    friend Vec4 operator+ <>(const Vec4& u, const Vec4& v);
    friend Vec4 operator- <>(const Vec4& u, const Vec4& v);
    friend Vec4 operator* <>(const Vec4<T>& u, const T s);
    friend Vec4 operator* <>(const T s, const Vec4<T>& u);
    friend Vec4 operator/ <>(const Vec4<T>& u, const T s);
    friend bool operator== <>(const Vec4& u, const Vec4& v);
    friend bool operator!= <>(const Vec4& u, const Vec4& v);
    friend T inner <>(const Vec4& u, const Vec4& v);
    friend Vec4 proj <>(const Vec4& u, const Vec4& v);
    friend Vec4 min <>(const Vec4& u, const Vec4& v);
    friend Vec4 max <>(const Vec4& u, const Vec4& v);
    friend Vec4 normalized <>(const Vec4& u);
    
    /// Return the zero vector.
    static const Vec4& zero();
    
    /// Return the unit vector with v(i) == 0 and v(index) == 1.
    static const Vec4 unit(const unsigned index);

    /// Return a description like "2D float vector", useful for debugging and messages.
    static std::string description();
    
    protected:
    T m_data[4];                  ///< The element storage
};

    
//
// Global function declarations
//

/// Print a vector to the output stream.
template <typename T>
std::ostream& operator<<(std::ostream& o, const Vec4<T>& v);
    
/// Read a vector from the input stream.
template <typename T>
std::istream& operator>>(std::istream& i, Vec4<T>& v);

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
Vec4<T>::Vec4() {
#if defined(LINALG_DEFAULT_INITIALIZE_TO_ZERO)
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] = 0;
#endif
}
    
template <typename T>
inline
Vec4<T>::Vec4(const T v0, const T v1, const T v2, const T v3) {
    m_data[0] = v0; m_data[1] = v1; m_data[2] = v2; m_data[3] = v3; 
}

template <typename T>
template <typename Other>
inline
Vec4<T>::Vec4(const Vec4<Other>& v) {
    for (int i = 0; i < 4; ++i)
        m_data[i] = (T) v(i);
}

template <typename T>
inline
Vec4<T>& Vec4<T>::operator=(const Vec4<T>& v) {
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] = v.m_data[i];
    return *this;
}

template <typename T>
inline
const T& Vec4<T>::operator()(const unsigned i) const {
    assert(i < 4);
    return m_data[i];
}

template <typename T>
inline
T& Vec4<T>::operator()(const unsigned i) {
    assert(i < 4);
    return m_data[i];
}

template <typename T>
inline
const T* Vec4<T>::data() const {
    return m_data;
}

template <typename T>
inline
T* Vec4<T>::data() {
    return m_data;
}

template <typename T>
inline 
typename Vec4<T>::const_iterator Vec4<T>::begin() const {
    return m_data;
}

template <typename T>
inline 
typename Vec4<T>::iterator Vec4<T>::begin() {
    return m_data;
}

template <typename T>
inline 
typename Vec4<T>::const_iterator Vec4<T>::end() const {
    return m_data + 4;
}

template <typename T>
inline 
typename Vec4<T>::iterator Vec4<T>::end() {
    return m_data + 4;
}

template <typename T>
inline
unsigned Vec4<T>::size() const {
    return 4;
}

template <typename T>
inline
void Vec4<T>::limit(const std::pair<T,T>& l) {
    for (unsigned i = 0; i < 4; ++i) {
        if (m_data[i] < l.first) 
            m_data[i] = l.first;
        if (m_data[i] > l.second) 
            m_data[i] = l.second;
    }
}

template <typename T>
inline
void Vec4<T>::limit(const T lower, const T upper) {
    for (unsigned i = 0; i < 4; ++i) {
        if (m_data[i] < lower) 
            m_data[i] = lower;
        if (m_data[i] > upper)
            m_data[i] = upper;
    }
}

template <typename T>
inline
void Vec4<T>::limit(const Vec4<T>& lower, const Vec4<T>& upper) {
    for (unsigned i = 0; i < 4; ++i) {
        if (m_data[i] < lower(i)) 
            m_data[i] = lower(i);
        if (m_data[i] > upper(i))
            m_data[i] = upper(i);
    }
}

template <typename T>
inline
Vec4<T>& Vec4<T>::operator+=(const Vec4<T>& v) {
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] += v.m_data[i];
    return *this;
}

template <typename T>
inline
Vec4<T>& Vec4<T>::operator-=(const Vec4<T>& v) {
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] -= v.m_data[i];
    return *this;
}

template <typename T>
template <typename Scalar>
inline
Vec4<T>& Vec4<T>::operator*=(const Scalar d) {
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] *= d;
    return *this;
}

template <typename T>
template <typename Scalar>
inline
Vec4<T>& Vec4<T>::operator/=(const Scalar d) {
    // Don't get too tricky -- think of T == int and d = 2.
    for (unsigned i = 0; i < 4; ++i)
        m_data[i] /= d;
    return *this;
}

template <typename T>
inline 
T Vec4<T>::sum() const {
    T w = 0;
    for (unsigned i = 0; i < 4; ++i)
        w += m_data[i];
    return w;
}

template <typename T>
inline 
T Vec4<T>::prod() const {
    T w = 1;
    for (unsigned i = 0; i < 4; ++i)
        w *= m_data[i];
    return w;
}

template <typename T>
inline 
T Vec4<T>::length() const {
    T w = 0;
    for (unsigned i = 0; i < 4; ++i)
        w += m_data[i] * m_data[i];
    return (T)sqrt((double)w);      // Need this to resolve ambiguity when T=int
}

template <typename T>
inline 
T Vec4<T>::lengthSqr() const {
    T w = 0;
    for (unsigned i = 0; i < 4; ++i)
        w += m_data[i] * m_data[i];
    return w;
}

template <typename T>
inline
Vec4<T> operator-(const Vec4<T>& u) {
    Vec4<T> v;
    for (unsigned i = 0; i < 4; ++i)
        v.m_data[i] = -u.m_data[i];
    return v;
}

template <typename T>
inline
Vec4<T> operator+(const Vec4<T>& u, const Vec4<T>& v) {
    Vec4<T> w;
    for (unsigned i = 0; i < 4; ++i) 
        w.m_data[i] = u.m_data[i] + v.m_data[i];
    return w;
}

template <typename T>
inline
Vec4<T> operator-(const Vec4<T>& u, const Vec4<T>& v) {
    Vec4<T> w;
    for (unsigned i = 0; i < 4; ++i) 
        w.m_data[i] = u.m_data[i] - v.m_data[i];
    return w;
}
    
template <typename T>
inline
Vec4<T> operator*(const Vec4<T>& u, const T s) {
    Vec4<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Vec4<T> operator*(const T s, const Vec4<T>& u) {
    Vec4<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Vec4<T> operator/(const Vec4<T>& u, const T s) {
    Vec4<T> w(u);
    w /= s;
    return w;
}
    
template <typename T>
inline
bool operator==(const Vec4<T>& u, const Vec4<T>& v) {
    for (unsigned i = 0; i < 4; ++i)
        if (u.m_data[i] != v.m_data[i])
            return false;
    return true;
}

template <typename T>
inline
bool operator!=(const Vec4<T>& u, const Vec4<T>& v) {
    for (unsigned i = 0; i < 4; ++i)
        if (u.m_data[i] != v.m_data[i])
            return true;
    return false;
}

template <typename T>
inline
T inner(const Vec4<T>& u, const Vec4<T>& v) {
    T w = 0;
    for (unsigned i = 0; i < 4; ++i)
        w += u.m_data[i] * v.m_data[i];
    return w;
}

template <typename T>
inline
Vec4<T> proj(const Vec4<T>& u, const Vec4<T>& v) {
    return v * inner(u, v) / v.lengthSqr();
}

template <typename T>
inline
Vec4<T> min(const Vec4<T>& u, const Vec4<T>& v) {
    Vec4<T> result;
    for (unsigned i = 0; i < 4; ++i)
        result.m_data[i] = (u.m_data[i] < v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Vec4<T> max(const Vec4<T>& u, const Vec4<T>& v) {
    Vec4<T> result;
    for (unsigned i = 0; i < 4; ++i)
        result.m_data[i] = (u.m_data[i] > v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Vec4<T> normalized(const Vec4<T>& u) {
    return u / u.length();
}



} // namespace

#endif