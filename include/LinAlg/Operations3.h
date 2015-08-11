/*
 * LinAlg: fixed-size vector and matrix library with LAPACK bindings.
 * Operations3.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_OPERATIONS3_H
#define LINALG_OPERATIONS3_H

namespace LinAlg {
    // Pre-declare template friends
    template <typename T> class Vec3;

    /// \name 3D-specific operators.
    /// @{
    
    /// Return the dot product of two 3D vectors.
    template <typename T>
    T dot(const Vec3<T>& v, const Vec3<T>& w) {
        T sum = 0;
        for (unsigned i = 0; i < 3; ++i)
            sum += v(i) * w(i);
        return sum;
    }

    /// Return the cross product of two 3D vectors.
    template <typename T>
    Vec3<T> cross(const Vec3<T>& v, const Vec3<T>& w) {
        Vec3<T> result;

        result(0) = v(1) * w(2) - v(2) * w(1);
        result(1) = v(2) * w(0) - v(0) * w(2);
        result(2) = v(0) * w(1) - v(1) * w(0);

        return result;
    }
    
    /// @}

} // namespace

#endif
