/*
 * LinAlg: fixed-size vector and matrix library with LAPACK bindings.
 * Transform3.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_AFFINE_TRANSFORM3_H
#define LINALG_AFFINE_TRANSFORM3_H

#include <LinAlg/Transform3Fwd.h>

#include <LinAlg/Mat4.h>
#include <LinAlg/Vec4.h>
#include <LinAlg/Mat3.h>
#include <LinAlg/Vec3.h>

namespace LinAlg {

    // If you want support for non-uniform scales or support for shears, 
    // you need to write the correct generators below and change the inverse()
    // function appropriately.

    /// A 3D affine transformation represented as a 4x4 matrix.
    /// This version supports translation, rotation, and uniform scales.
    template <typename T>
    class Transform3: public Mat4<T> {
        public:

        /// The identity transform.
        Transform3(): Mat4<T>(1, 0, 0, 0, 
                              0, 1, 0, 0, 
                              0, 0, 1, 0, 
                              0, 0, 0, 1) {}

        /// Initialize from a matrix of values.
        Transform3(const Mat4<T>& values): Mat4<T>(values) {}

        /// Initialize from a rotation/scaling matrix and a translation vector.
        Transform3(const Mat3<T>& R, const Vec3<T>& t):
            Mat4<T>(R(0,0), R(1,0), R(2,0), 0, 
                    R(0,1), R(1,1), R(2,1), 0, 
                    R(0,2), R(1,2), R(2,2), 0, 
                    t(0), t(1), t(2), 1) {}

        /// Copy an affine transform.
        Transform3(const Transform3& other): Mat4<T>(other) {}

        /// Assign an affine transform to this one.
        Transform3& operator=(const Transform3& other) {
            Mat4<T>::operator=(other);
            return *this;
        }

        /// \name Various affine transforms creators.
        /// @{

        /// Return the identity transform.
        static Mat4<T> identity() {
            return Mat4<T>(1, 0, 0, 0, 
                           0, 1, 0, 0, 
                           0, 0, 1, 0, 
                           0, 0, 0, 1);
        }

        /// Return a translation transform.
        static Mat4<T> translation(const Vec3<T>& v) {
            return Mat4<T>(1, 0, 0, 0, 
                           0, 1, 0, 0, 
                           0, 0, 1, 0, 
                           v(0), v(1), v(2), 1);
        }

        /// Return a uniform scaling transform.
        static Mat4<T> scale(const T& s) {
            return Mat4<T>(s, 0, 0, 0, 
                           0, s, 0, 0, 
                           0, 0, s, 0, 
                           0, 0, 0, 1);
        }

        /// Return a rotation transform of \c angle radians about the unit vector \c axis.
        /// \note This clearly won't work very well for integer types \c T.
        static Mat4<T> rotation(const T angle, const Vec3<T>& v) {
            const T sa = sin(angle), ca = cos(angle);
            return Mat4<T>(v[0] * v[0] * (1 - ca) +        ca, v[0] * v[1] * (1 - ca) - v[2] * sa, v[0] * v[2] * (1 - ca) + v[1] * sa, 0,
                           v[1] * v[0] * (1 - ca) + v[2] * sa, v[1] * v[1] * (1 - ca) +        ca, v[1] * v[2] * (1 - ca) - v[0] * sa, 0,
                           v[2] * v[0] * (1 - ca) - v[1] * sa, v[2] * v[1] * (1 - ca) + v[0] * sa, v[2] * v[2] * (1 - ca) +        ca, 0, 
                           0,                               0,                                  0,                                  1);
        }

        /// Return a rotation transform of <code>axis.length()</code> radians about the unit 
        /// vector in the direction of \c axis.
        /// \note This clearly won't work very well for integer types \c T.
        static Mat4<T> rotation(const Vec3<T>& axisAngle) {
            const T len = axisAngle.length();
            assert(len > 0);

            return rotation(len, axis / len);
        }

        // Add shear?

        /// @}

        /// Return the uniform scaling factor.
        T getScale() const {
            return sqrt(getScaleSqr());
        }

        /// Return the uniform scaling factor squared.
        T getScaleSqr() const {
            const Transform3& M = *this;
            return M(0,0) * M(0,0) + M(1,0) * M(1,0) + M(2,0) * M(2,0);
        }

        /// Return the translation.
        Vec3<T> getTranslation() const {
           const Transform3& M = *this;
           return Vec3<T>(M(0,3), M(1,3), M(2,3));
        }

        /// Calculate the inverse of the transfromation.
        Transform3 inverse() const {
            const Transform3& M = *this;

            // The following only works for translations, rotations and uniform scales!
            assert(M(3,0) == 0 && M(3,1) == 0 && M(3,2) == 0 && M(3,3) == 1);

            // Find the scale factor -- length of any of the row or column vectors of the rotation part.
            const T scaleSqr = getScaleSqr();
            assert(scaleSqr > 0);

            // Construct the inverse rotation/scale matrix
            // The 3x3 upper-left block is sR, where s is the scale and R is the rotation matrix.
            // We need (1/s) R^T, so the following computes (1/s^2) * (sR)^T.
            const T s = 1 / scaleSqr;
            const Mat3<T> invRotScale(s * M(0,0), s * M(0,1), s * M(0,2),
                                      s * M(1,0), s * M(1,1), s * M(1,2),
                                      s * M(2,0), s * M(2,1), s * M(2,2));

            // Construct the inverse translation vector
            const Vec3<T> invTranslation(-M(0,3), -M(1,3), -M(2,3));

            // Build the final matrix
            return Transform3(invRotScale, invRotScale * invTranslation);
        }

        /// Transform a point
        Vec3<T> transPoint(const Vec3<T>& v) const {
            const Transform3<T>& M = *this;
            return Vec3<T>(M(0,0) * v(0) + M(0,1) * v(1) + M(0,2) * v(2) + M(0,3),
                           M(1,0) * v(0) + M(1,1) * v(1) + M(1,2) * v(2) + M(1,3),
                           M(2,0) * v(0) + M(2,1) * v(1) + M(2,2) * v(2) + M(2,3));
        }

        /// Transform a vector
        Vec3<T> transVector(const Vec3<T>& v) const {
            const Transform3<T>& M = *this;
            return Vec3<T>(M(0,0) * v(0) + M(0,1) * v(1) + M(0,2) * v(2),
                           M(1,0) * v(0) + M(1,1) * v(1) + M(1,2) * v(2),
                           M(2,0) * v(0) + M(2,1) * v(1) + M(2,2) * v(2));
        }

    };
}

#endif


