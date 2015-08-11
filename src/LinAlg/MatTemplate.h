/*
 * LinAlg 3.0: fixed-sized vector and matrix library for small dimensions with optional LAPACK bindings.
 * Mat$D.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_MAT$D_H
#define LINALG_MAT$D_H

#include <iosfwd>
#include <string>

#include <LinAlg/Vec$DFwd.h>
#include <LinAlg/Mat$DFwd.h>

namespace LinAlg {

//
// Pre-declare template friends
//

/// Return the negation of a matrix.
template <typename T>
Mat$D<T> operator-(const Mat$D<T>& u);

/// Add two matrices together.
template <typename T>
Mat$D<T> operator+(const Mat$D<T>& u, const Mat$D<T>& v);

/// Subtract two matrices.
template <typename T>
Mat$D<T> operator-(const Mat$D<T>& u, const Mat$D<T>& v);

/// Multiply two matrices.
template <typename T>
Mat$D<T> operator*(const Mat$D<T>& u, const Mat$D<T>& v);
    
/// Right-multiply a matrix by a scalar.
template <typename T>
Mat$D<T> operator*(const Mat$D<T>& u, const T s);

/// Left-multiply a matrix by a scalar.
template <typename T>
Mat$D<T> operator*(const T s, const Mat$D<T>& u);

/// Divide a matrix by a scalar.
template <typename T>
Mat$D<T> operator/(const Mat$D<T>& u, const T s);
    
/// Multiply a vector by a matrix.
template <typename T>
Vec$D<T> operator*(const Mat$D<T>& u, const Vec$D<T>& v);

/// Exact equality -- use with caution on floating-point base types.
template <typename T> 
bool operator==(const Mat$D<T>& u, const Mat$D<T>& v);

/// Exact inequality -- use with caution on floating-point base types.
template <typename T> 
bool operator!=(const Mat$D<T>& u, const Mat$D<T>& v);

/// Per-element minimum.
template <typename T>
Mat$D<T> min(const Mat$D<T>& u, const Mat$D<T>& v);

/// Per-element maximum.
template <typename T>
Mat$D<T> max(const Mat$D<T>& u, const Mat$D<T>& v);


/// \class Mat$D Mat$D.h <LinAlg/Mat$D.h>
/// A linear algebra matrix class of fixed type and size $D.
template <typename T>
class Mat$D {
    public:
    /// This type
    typedef Mat$D<T> self_type;
    
    /// The vector type of the same dimension.
    typedef Vec$D<T> vec_type;

    /// The underlying data type
    typedef T value_type;
    
    /// The number of dimensions
    enum { num_dims = $D };

    /// The number of elements
    enum { num_elements = $D * $D };
    
    /// An iterator
    typedef T* iterator;
    
    /// A constant iterator.
    typedef const T* const_iterator;
    
    
    /// Create an uninitialized matrix.
    Mat$D();

    /// Create a matrix from named values.
    Mat$D($MATRIX_CONSTRUCTOR_ARGUMENTS);
    
    /// Copy another matrix of possibly a different type.
    template <typename Other>
    Mat$D(const Mat$D<Other>& m);
    
    /// Assign another matrix to this matrix.
    Mat$D& operator=(const Mat$D& m);

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
    void limit(const Mat$D& lower, const Mat$D& upper);
    
    /// Return the sum of all the elements.
    T sum() const;
    
    /// Return the product of all the elements.
    T prod() const;
    
    /// Return the transpose of this matrix.
    Mat$D trans() const;
    
    /// Get the \c i'th column of the matrix.
    Vec$D<T> getCol(const unsigned i) const;
    
    /// Set the \c i'th column of the matrix.
    void setCol(const unsigned i, const Vec$D<T>& v);

    /// Get the \c i'th row of the matrix.
    Vec$D<T> getRow(const unsigned i) const;

    /// Set the \c i'th row of the matrix.
    void setRow(const unsigned i, const Vec$D<T>& v);

    /// Compute the eigenvalues and eigenvectors.
    /// The eigenvectors are stored as the columns of \c vectors.
    void eigs(Vec$D<T>& values, Mat$D& vectors) const;

    /// Compare entries with another matrix, within eps
    bool equal(const Mat$D<T>& v, const T eps) const;

    /// Add another matrix to this one.
    Mat$D& operator+=(const Mat$D& m);

    /// Subtract another matrix from this one.
    Mat$D& operator-=(const Mat$D& m);
    
    /// Multiply this matrix by a scalar value.
    Mat$D& operator*=(const T d);
    
    /// Divided this matrix by a scalar value.
    Mat$D& operator/=(const T d);

    // Dear god the things you have to do to get a templated friend function
    friend Mat$D operator-<>(const Mat$D& u);
    friend Mat$D operator+<>(const Mat$D& u, const Mat$D& v);
    friend Mat$D operator-<>(const Mat$D& u, const Mat$D& v);
    friend Mat$D operator*<>(const Mat$D& u, const Mat$D& v);
    friend Mat$D<T> operator* <>(const Mat$D<T>& u, const T s);
    friend Mat$D<T> operator* <>(const T s, const Mat$D<T>& u);
    friend Mat$D<T> operator/ <>(const Mat$D<T>& u, const T s);
    friend Vec$D<T> operator*<>(const Mat$D& m, const Vec$D<T>& v);
    friend bool operator== <>(const Mat$D& u, const Mat$D& v);
    friend bool operator!= <>(const Mat$D& u, const Mat$D& v);
    friend Mat$D min<>(const Mat$D& u, const Mat$D& v);
    friend Mat$D max<>(const Mat$D& u, const Mat$D& v);

    /// Return the zero matrix.
    static const Mat$D& zero();
    
    /// Return the identity matrix.
    static const Mat$D& ident();

    /// Calculate the inverse of a matrix from its eigenvalues and 
    /// eigenvectors.
    static Mat$D inverse(const Vec$D<T>& values, const Mat$D& vectors);

    /// Return a description like "2D float", useful for debugging and messages.
    static std::string description();
    
    protected:
    T m_data[$D * $D];                                      ///< The element storage.
};

    
//
// Global functions
//

/// Print a matrix to an output stream.
template <typename T>
std::ostream& operator<<(std::ostream& o, const Mat$D<T>& m);

template <typename T>
std::istream& operator>>(std::istream& i, Mat$D<T>& m);

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
Mat$D<T>::Mat$D() {
#if defined(LINALG_DEFAULT_INITIALIZE_TO_ZERO)
    for (unsigned i = 0; i < $D * $D; ++i)
        m_data[i] = 0;
#endif    
}
    
template <typename T>
inline
Mat$D<T>::Mat$D($MATRIX_CONSTRUCTOR_ARGUMENTS) {
    $MATRIX_CONSTRUCTOR_BODY
}

template <typename T>
template <typename Other>
inline 
Mat$D<T>::Mat$D(const Mat$D<Other>& m) {
    for (unsigned i = 0; i < $D * $D; ++i) {
        m_data[i] = (T) m.m_data[i];
    }
}

template <typename T>
inline 
Mat$D<T>& Mat$D<T>::operator=(const Mat$D<T>& m) {
    for (unsigned i = 0; i < $D * $D; ++i) {
        m_data[i] = m.m_data[i];
    }
    return *this;
}

template <typename T>
inline
const T& Mat$D<T>::operator()(const unsigned i, const unsigned j) const {
    assert(i < $D && j < $D);
    return m_data[j * $D + i];
}

template <typename T>
inline
T& Mat$D<T>::operator()(const unsigned i, const unsigned j) {
    assert(i < $D && j < $D);
    return m_data[j * $D + i];
}

template <typename T>
inline
const T* Mat$D<T>::data() const {
    return m_data;
}

template <typename T>
inline
T* Mat$D<T>::data() {
    return m_data;
}

template <typename T>
inline 
typename Mat$D<T>::const_iterator Mat$D<T>::begin() const {
    return m_data;
}

template <typename T>
inline 
typename Mat$D<T>::iterator Mat$D<T>::begin() {
    return m_data;
}

template <typename T>
inline 
typename Mat$D<T>::const_iterator Mat$D<T>::end() const {
    return m_data + $D * $D;
}

template <typename T>
inline 
typename Mat$D<T>::iterator Mat$D<T>::end() {
    return m_data + $D * $D;
}

template <typename T>
inline
unsigned Mat$D<T>::size1() const {
    return $D;
}

template <typename T>
inline
unsigned Mat$D<T>::size2() const {
    return $D;
}

template <typename T>
inline
void Mat$D<T>::limit(const std::pair<T,T>& l) {
    for (unsigned i = 0; i < $D* $D; ++i) {
        if (m_data[i] < l.first) 
            m_data[i] = l.first;
        if (m_data[i] > l.second) 
            m_data[i] = l.second;
    }
}

template <typename T>
inline
void Mat$D<T>::limit(const T lower, const T upper) {
    for (unsigned i = 0; i < $D* $D; ++i) {
        if (m_data[i] < lower) 
            m_data[i] = lower;
        if (m_data[i] > upper)
            m_data[i] = upper;
    }
}

template <typename T>
inline
void Mat$D<T>::limit(const Mat$D<T>& lower, const Mat$D<T>& upper) {
    for (unsigned j = 0; j < $D; ++j) {
        for (unsigned i = 0; i < $D; ++i) {
            if ((*this)(i,j) < lower(i,j)) 
                (*this)(i,j) = lower(i,j);
            if ((*this)(i,j) > upper(i,j))
                (*this)(i,j) = upper(i,j);
        }
    }
}

template <typename T>
inline 
T Mat$D<T>::sum() const {
    T w = 0;
    for (unsigned i = 0; i < $D* $D; ++i)
        w += m_data[i];
    return w;
}

template <typename T>
inline 
T Mat$D<T>::prod() const {
    T w = 1;
    for (unsigned i = 0; i < $D* $D; ++i)
        w *= m_data[i];
    return w;
}

template <typename T>
inline 
Mat$D<T>& Mat$D<T>::operator+=(const Mat$D<T>& m) {
    for (unsigned i = 0; i < $D * $D; ++i)
        m_data[i] += m.m_data[i];
    return *this;
}

template <typename T>
inline 
Mat$D<T>& Mat$D<T>::operator-=(const Mat$D<T>& m) {
    for (unsigned i = 0; i < $D * $D; ++i)
        m_data[i] -= m.m_data[i];
    return *this;
}

template <typename T>
inline 
Mat$D<T>& Mat$D<T>::operator*=(const T d) {
    for (unsigned i = 0; i < $D * $D; ++i)
        m_data[i] *= d;
    return *this;
}

template <typename T>
inline 
Mat$D<T>& Mat$D<T>::operator/=(const T s) {
    // Don't get too tricky -- think of T == int and s = 2.
    for (unsigned i = 0; i < $D * $D; ++i)
        m_data[i] /= s;
    return *this;
}

template <typename T>
inline 
Mat$D<T> operator-(const Mat$D<T>& u) {
    Mat$D<T> w(u);
    for (unsigned i = 0; i < $D * $D; ++i)
        w.m_data[i] = -u.m_data[i];
    return w;
}

template <typename T>
inline 
Mat$D<T> operator+(const Mat$D<T>& u, const Mat$D<T>& v) {
    Mat$D<T> w(u);
    for (unsigned i = 0; i < $D * $D; ++i)
        w.m_data[i] = u.m_data[i] + v.m_data[i];
    return w;
}

template <typename T>
inline 
Mat$D<T> operator-(const Mat$D<T>& u, const Mat$D<T>& v) {
    Mat$D<T> w(u);
    for (unsigned i = 0; i < $D * $D; ++i)
        w.m_data[i] = u.m_data[i] - v.m_data[i];
    return w;
}
    
template <typename T>
inline
Mat$D<T> operator*(const Mat$D<T>& u, const T s) {
    Mat$D<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Mat$D<T> operator*(const T s, const Mat$D<T>& u) {
    Mat$D<T> w(u);
    w *= s;
    return w;
}

template <typename T>
inline
Mat$D<T> operator/(const Mat$D<T>& u, const T s) {
    Mat$D<T> w(u);
    w /= s;
    return w;
}
    
template <typename T>
inline
bool operator==(const Mat$D<T>& u, const Mat$D<T>& v) {
    for (unsigned i = 0; i < $D * $D; ++i)
        if (u.m_data[i] != v.m_data[i])
            return false;
    return true;
}

template <typename T>
inline
bool operator!=(const Mat$D<T>& u, const Mat$D<T>& v) {
    for (unsigned i = 0; i < $D * $D; ++i)
        if (u.m_data[i] != v.m_data[i])
            return true;
    return false;
}

template <typename T>
inline
Mat$D<T> min(const Mat$D<T>& u, const Mat$D<T>& v) {
    Mat$D<T> result;
    for (unsigned i = 0; i < $D * $D; ++i)
        result.m_data[i] = (u.m_data[i] < v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}

template <typename T>
inline
Mat$D<T> max(const Mat$D<T>& u, const Mat$D<T>& v) {
    Mat$D<T> result;
    for (unsigned i = 0; i < $D * $D; ++i)
        result.m_data[i] = (u.m_data[i] > v.m_data[i] ? u.m_data[i] : v.m_data[i]);
    return result;
}


} // namespace

#endif
