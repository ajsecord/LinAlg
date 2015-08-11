/*
 * LinAlg: fixed-size vector and matrix library with LAPACK bindings.
 * Transform2.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_AFFINE_TRANSFORM2_H
#define LINALG_AFFINE_TRANSFORM2_H

#include <LinAlg/Transform2Fwd.h>

#include <LinAlg/Mat3.h>
#include <LinAlg/Vec3.h>
#include <LinAlg/Mat2.h>
#include <LinAlg/Vec2.h>

namespace LinAlg {

    // If you want support for non-uniform scales or support for shears, 
    // you need to write the correct generators below and change the inverse()
    // function appropriately.

    /// A 2D affine transformation represented as a 3x3 matrix.
    /// This version supports translation, rotation, and uniform scales.
    template <typename T>
    class Transform2: public Mat3<T> {
        public:

        /// The identity transform.
        Transform2(): Mat3<T>(1, 0, 0, 
                              0, 1, 0, 
                              0, 0, 1) {}

        /// Initialize from a matrix of values.
        Transform2(const Mat3<T>& values): Mat3<T>(values) {}

        /// Initialize from a rotation/scaling matrix and a translation vector.
        Transform2(const Mat2<T>& R, const Vec2<T>& t):
            Mat3<T>(R(0,0), R(1,0), 0, 
                    R(0,1), R(1,1), 0, 
                    t(0),   t(1), 1) {}

        /// Copy an affine transform.
        Transform2(const Transform2& other): Mat3<T>(other) {}

        /// Assign an affine transform to this one.
        Transform2& operator=(const Transform2& other) {
            Mat3<T>::operator=(other);
            return *this;
        }

        /// \name Various affine transforms creators.
        /// @{

        /// Return the identity transform.
        static Mat3<T> identity() {
            return Mat3<T>(1, 0, 0, 
                           0, 1, 0, 
                           0, 0, 1);
        }

        /// Return a translation transform.
        static Mat3<T> translation(const Vec2<T>& v) {
            return Mat3<T>(1, 0, 0, 
                           0, 1, 0, 
                           v(0), v(1), 1);
        }

        /// Return a uniform scaling transform.
        static Mat3<T> scale(const T& s) {
            return Mat3<T>(s, 0, 0, 
                           0, s, 0, 
                           0, 0, 1);        
        }

        /// Return a rotation transform of \c angle radians about the origin.
        /// \note This clearly won't work very well for integer types \c T.
        static Mat3<T> rotation(const T angle) {
            const T sa = sin(angle), ca = cos(angle);
            return Mat3<T>( ca, sa, 0, 
                           -sa, ca, 0, 
                             0,  0, 1);
        }

        // Add shear?

        /// @}

        /// Return the uniform scaling factor.
        T getScale() const {
            return sqrt(getScaleSqr());
        }

        /// Return the uniform scaling factor squared.
        T getScaleSqr() const {
            const Transform2& M = *this;
            return M(0,0) * M(0,0) + M(1,0) * M(1,0);
        }

        /// Return the translation.
        Vec2<T> getTranslation() const {
           const Transform2& M = *this;
           return Vec2<T>(M(0,2), M(1,2));
        }

        /// Calculate the inverse of the transfromation.
        Transform2 inverse() const {
            const Transform2& M = *this;

            // The following only works for translations, rotations and uniform scales!
            assert(M(2,0) == 0 && M(2,1) == 0 && M(2,2) == 1);

            // Find the scale factor -- length of any of the row or column vectors of the rotation part.
            const T scaleSqr = getScaleSqr();
            assert(scaleSqr > 0);

            // Construct the inverse rotation/scale matrix
            // The 2x2 upper-left block is sR, where s is the scale and R is the rotation matrix.
            // We need (1/s) R^T, so the following computes (1/s^2) * (sR)^T.
            const T s = 1 / scaleSqr;
            const Mat2<T> invRotScale(s * M(0,0), s * M(0,1), 
                                      s * M(1,0), s * M(1,1));

            // Construct the inverse translation vector
            const Vec2<T> invTranslation(-M(0,2), -M(1,2));

            // Build the final matrix
            return Transform2(invRotScale, invRotScale * invTranslation);
        }

        /// Transform a point
        Vec2<T> transPoint(const Vec2<T>& v) const {
            const Transform2<T>& M = *this;
            return Vec2<T>(M(0,0) * v(0) + M(0,1) * v(1) + M(0,2), 
                              M(1,0) * v(0) + M(1,1) * v(1) + M(1,2));
        }

        /// Transform a vector
        Vec2<T> transVector(const Vec2<T>& v) const {
            const Transform2<T>& M = *this;
            return Vec2<T>(M(0,0) * v(0) + M(0,1) * v(1), 
                              M(1,0) * v(0) + M(1,1) * v(1));
        }
    };


}

#endif


