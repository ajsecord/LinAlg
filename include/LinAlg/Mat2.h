// Warning: this file was generated from src/LinAlg/MatTemplate.h.  Do not edit if you value your changes!

/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * Mat2.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_MAT2_H
#define LINALG_MAT2_H

#include <iosfwd>
#include <string>

#include <LinAlg/Vec2Fwd.h>
#include <LinAlg/Mat2Fwd.h>

namespace LinAlg {

//
// Pre-declare template friends
//

/// Return the negation of a matrix.
template <typename T>
Mat2<T> operator-(const Mat2<T>& u);

/// Add two matrices together.
template <typename T>
Mat2<T> operator+(const Mat2<T>& u, const Mat2<T>& v);

/// Subtract two matrices.
template <typename T>
Mat2<T> operator-(const Mat2<T>& u, const Mat2<T>& v);

/// Multiply two matrices.
template <typename T>
Mat2<T> operator*(const Mat2<T>& u, const Mat2<T>& v);
    
/// Right-multiply a matrix by a scalar.
template <typename T>
Mat2<T> operator*(const Mat2<T>& u, const T s);

/// Left-multiply a matrix by a scalar.
template <typename T>
Mat2<T> operator*(const T s, const Mat2<T>& u);

/// Divide a matrix by a scalar.
template <typename T>
Mat2<T> operator/(const Mat2<T>& u, const T s);
    
/// Multiply a vector by a matrix.
template <typename T>
Vec2<T> operator*(const Mat2<T>& u, const Vec2<T>& v);

/// Exact equality -- use with caution on floating-point base types.
template <typename T> 
bool operator==(const Mat2<T>& u, const Mat2<T>& v);

/// Exact inequality -- use with caution on floating-point base types.
template <typename T> 
bool operator!=(const Mat2<T>& u, const Mat2<T>& v);

/// Per-element minimum.
template <typename T>
Mat2<T> min(const Mat2<T>& u, const Mat2<T>& v);

/// Per-element maximum.
template <typename T>
Mat2<T> max(const Mat2<T>& u, const Mat2<T>& v);


/// \class Mat2 Mat2.h <LinAlg/Mat2.h>
/// A linear algebra matrix class of fixed type and size 2.
template <typename T>
class Mat2 {
    public:
    /// This type
    typedef Mat2<T> self_type;
    
    /// The vector type of the same dimension.
    typedef Vec2<T> vec_type;

    /// The underlying data type
    typedef T value_type;
    
    /// The number of dimensions
    enum { num_dims = 2 };

    /// The number of elements
    enum { num_elements = 2 * 2 };
    
    /// An iterator
    typedef T* iterator;
    
    /// A constant iterator.
    typedef const T* const_iterator;
    
    
    /// Create an uninitialized matrix.
    Mat2();

    /// Create a matrix from named values.
    Mat2(const T v00, const T v10, const T v01, const T v11);
    
    /// Copy another matrix of possibly a different type.
    template <typename Other>
    Mat2(const Mat2<Other>& m);
    
    /// Assign another matrix to this matrix.
    Mat2& operator=(const Mat2& m);

    /// Read-only element access.
    const T& operator()(const unsigned i, const unsigned j) const;
    
    /// Read-write element access.
    T& operator()(const unsigned i, const unsigned j);

    /// Direct read-only access to the underlying m_data.  The elements are stored densely
    /// in column-major order.
    const T* data() const;
    
    /// Direct read-write access to the underlying m_data.  The elements are stored densely
    /// in column-major order.
    T* data();

    /// Get a read-only iterator to the first element.
    const_iterator begin() const;
    
    /// Get a read-write pointer to the first element.
    iterator begin();
    
    /// Get a read-only iterator to past the last element.
    const_iterator end() const;
    
    /// Get a read-write iterator to past the last element.
    iterator end();
    
    /// The number of rows.
    unsigned size1() const;
    
    /// The number of columns.
    unsigned size2() const;

    /// Limit each element to the range [l.first, l.second]
    void limit(const std::pair<T,T>& l);
    
    /// Limit each element to the range [lower, upper]
    void limit(const T lower, const T upper);
    
    /// Limit element (i,j) to the range [lower(i,j), upper(i,j)]
    void limit(const Mat2& lower, const Mat2& upper);
    
    /// Return the sum of all the elements.
    T sum() const;
    
    /// Return the product of all the elements.
    T prod() const;
    
    /// Return the transpose of this matrix.
    Mat2 trans() const;
    
    /// Get the \c i'th column of the matrix.
    Vec2<T> getCol(const unsigned i) const;
    
    /// Set the \c i'th column of the matrix.
    void setCol(const unsigned i, const Vec2<T>& v);

    /// Get the \c i'th row of the matrix.
    Vec2<T> getRow(const unsigned i) const;

    /// Set the \c i'th row of the matrix.
    void setRow(const unsigned i, const Vec2<T>& v);

    /// Compute the eigenvalues and eigenvectors.
    /// The eigenvectors are stored as the columns of \c vectors.
    void eigs(Vec2<T>& values, Mat2& vectors) const;

    /// Compare entries with another matrix, within eps
    bool equal(const Mat2<T>& v, const T eps) const;

    /// Add another matrix to this one.
    Mat2& operator+=(const Mat2& m);

    /// Subtract another matrix from this one.
    Mat2& operator-=(const Mat2& m);
    
    /// Multiply this matrix by a scalar value.
    Mat2& operator*=(const T d);
    
    /// Divided this matrix by a scalar value.
    Mat2& operator/=(const T d);

    // Dear god the things you have to do to get a templated friend function
    friend Mat2 operator-<>(const Mat2& u);
    friend Mat2 operator+<>(const Mat2& u, const Mat2& v);
    friend Mat2 operator-<>(const Mat2& u, const Mat2& v);
    friend Mat2 operator*<>(const Mat2& u, const Mat2& v);
    friend Mat2<T> operator* <>(const Mat2<T>& u, const T s);
    friend Mat2<T> operator* <>(const T s, const Mat2<T>& u);
    friend Mat2<T> operator/ <>(const Mat2<T>& u, const T s);
    friend Vec2<T> operator*<>(const Mat2& m, const Vec2<T>& v);
    friend bool operator== <>(const Mat2& u, const Mat2& v);
    friend bool operator!= <>(const Mat2& u, const Mat2& v);
    friend Mat2 min<>(const Mat2& u, const Mat2& v);
    friend Mat2 max<>(const Mat2& u, const Mat2& v);

    /// Return the zero matrix.
    static const Mat2& zero();
    
    /// Return the identity matrix.
    static const Mat2& ident();

    /// Calculate the inverse of a matrix from its eigenvalues and 
    /// eigenvectors.
    static Mat2 inverse(const Vec2<T>& values, const Mat2& vectors);

    /// Return a description like "2D float", useful for debugging and messages.
    static std::string description();
    
    protected:
    T m_data[2 * 2];                                      ///< The element storage.
};

    
//
// Global functions
//

/// Print a matrix to an output stream.
template <typename T>
std::ostream& operator<<(std::ostream& o, const Mat2<T>& m);

template <typename T>
std::istream& operator>>(std::istream& i, Mat2<T>& m);

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
Mat2<T>::Mat2() {
#if defined(LINALG_DEFAULT_INITIALIZE_TO_ZERO)
    for (unsigned i = 0; i < 2 * 2; ++i)
        m_data[i] = 0;
#endif    
}
    
template <typename T>
inline
Mat2<T>::Mat2(const T v00, const T v10, const T v01, const T v11) {
    m_data[0] = v00; m_data[1] = v10; m_data[2] = v01; m_data[3] = v11; 
}

template <typename T>
template <typename Other>
inline 
Mat2<T>::Mat2(const Mat2<Other>& m) {
    for (unsigned i = 0; i < 2 * 2; ++i) {
        m_data[i] = (T) m.m_data[i];
    }
}

template <typename T>
inline 
Mat2<T>& Mat2<T>::operator=(const Mat2<T>& m) {
    for (unsigned i = 0; i < 2 * 2; ++i) {
        m_data[i] = m.m_data[i];
    }
    return *this;
}

template <typename T>
inline
const T& Mat2<T>::operator()(const unsigned i, const unsigned j) const {
    assert(i < 2 && j < 2);
    return m_data[j * 2 + i];
}

template <typename T>
inline
T& Mat2<T>::operator()(const unsigned i, const unsigned j) {
    assert(i < 2 && j < 2);
    return m_data[j * 2 + i];
}

template <typename T>
inline
const T* Mat2<T>::data() const {
    return m_data;
}

template <typename T>
inline
T* Mat2<T>::data() {
    return m_data;
}

template <typename T>
inline 
typename Mat2<T>::const_iterator Mat2<T>::begin() const {
    return m_data;
}

template <typename T>
inline 
typename Mat2<T>::iterator Mat2<T>::begin() {
    return m_data;
}

template <typename T>
inline 
typename Mat2<T>::const_iterator Mat2<T>::end() const {
    return m_data + 2 * 2;
}

template <typename T>
inline 
typename Mat2<T>::iterator Mat2<T>::end() {
    return m_data + 2 * 2;
}

template <typename T>
inline
unsigned Mat2<T>::size1() const {
    return 2;
}

template <typename T>
inline
unsigned Mat2<T>::size2() const {
    return 2;
}

template <typename T>
inline
void Mat2<T>::limit(const std::pair<T,T>& l) {
    for (unsigned i = 0; i < 2* 2; ++i) {
        if (m_data[i] < l.first) 
            m_data[i] = l.first;
        if (m_data[i] > l.second) 
            m_data[i] = l.second;
    }
}

template <typename T>
inline
void Mat2<T>::limit(const T lower, const T upper) {
    for (unsigned i = 0; i < 2* 2; ++i) {
        if (m_data[i] < lower) 
            m_data[i] = lower;
        if (m_data[i] > upper)
            m_data[i] = upper;
    }
}

template <typename T>
inline
void Mat2<T>::limit(const Mat2<T>& lower, const Mat2<T>& upper) {
    for (unsigned j = 0; j < 2; ++j) {
        for (unsigned i = 0; i < 2; ++i) {
            if ((*this)(i,j) < lower(i,j)) 
                (*this)(i,j) = lower(i,j);
            if ((*this)(i,j) > upper(i,j))
                (*this)(i,j) = upper(i,j);
        }
    }
}

template <typename T>
inline 
T Mat2<T>::sum() const {
    T w = 0;
    for (unsigned i = 0; i < 2* 2; ++i)
        w += m_data[i];
    return w;
}

template <typename T>
inline 
T Mat2<T>::prod() const {
    T w = 1;
    for (unsigned i = 0; i < 2* 2; ++i)
        w *= m_data[i];
    return w;
}

template <typename T>
inline 
Mat2<T>& Mat2<T>::operator+=(const Mat2<T>& m) {
    for (unsigned i = 0; i < 2 * 2; ++i)
        m_data[i] += m.m_data[i];
    return *this;
}

template <typename T>
inline 
Mat2<T>& Mat2<T>::operator-=(const Mat2<T>& m) {
    for (unsigned i = 0; i < 2 * 2; ++i)
        m_data[i] -= m.m_data[i];
    return *this;
}

template <typename T>
inline 
Mat2<T>& Mat2<T>::operator*=(const T d) {
    for (unsigned i = 0; i < 2 * 2; ++i)
        m_data[i] *= d;
    return *this;
}

template <typename T>
inline 
Mat2<T>& Mat2<T>::operator/=(const T s) {
    // Don't get too tricky -- think of T == int and s = 2.
    for (unsigned i = 0; i < 2 * 2; ++i)
        m_data[i] /= s;
    return *this;
}

template <typename T>
inline 
Mat2<T> operator-(const Mat2<T>& u) {
    Mat2<T> w(u);
    for (unsigned i = 0; i < 2 * 2; ++i)
        w.m_data[i] = -u.m_data[i];
    return w;
}

template <typename T>
inline 
Mat2<T> operator+(const Mat2<T>& u, const Mat2<T>& v) {
    Mat2<T> w(u);
    for (unsigned i = 0; i < 2 * 2; ++i)
        w.m_data[i] = u.m_data[i] + v.m_data[i];
    return w;
}

template <typename T>
inline 
Mat2<T> operator-(const Mat2<T>& u, const Mat2<T>& v) {
    Mat2<T> w(u);
    for (unsigned i = 0; i < 2 * 2; ++i)
        w.m_data[i] = u.m_data[i] - v.m_data[i];
    return w;
}
    
template <typename T>
inline
Mat2<T> operator*(const Mat2<T>& u, const T s) {
    Mat2<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Mat2<T> operator*(const T s, const Mat2<T>& u) {
    Mat2<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Mat2<T> operator/(const Mat2<T>& u, const T s) {
    Mat2<T> w(u);
    w /= s;
    return w;
}
    
template <typename T>
inline
bool operator==(const Mat2<T>& u, const Mat2<T>& v) {
    for (unsigned i = 0; i < 2 * 2; ++i)
        if (u.m_data[i] != v.m_data[i])
            return false;
    return true;
}

template <typename T>
inline
bool operator!=(const Mat2<T>& u, const Mat2<T>& v) {
    for (unsigned i = 0; i < 2 * 2; ++i)
        if (u.m_data[i] != v.m_data[i])
            return true;
    return false;
}

template <typename T>
inline
Mat2<T> min(const Mat2<T>& u, const Mat2<T>& v) {
    Mat2<T> result;
    for (unsigned i = 0; i < 2 * 2; ++i)
        result.m_data[i] = (u.m_data[i] < v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Mat2<T> max(const Mat2<T>& u, const Mat2<T>& v) {
    Mat2<T> result;
    for (unsigned i = 0; i < 2 * 2; ++i)
        result.m_data[i] = (u.m_data[i] > v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}


} // namespace

#endif
