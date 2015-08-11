/*
 * LinAlg: A fixed-length vector and matrix library.
 * TypeHelpers.h
 * Copyright 2005-2008 Adrian Secord.
 * Contact: <a href="http://www.google.com/search?q=%22Adrian%20Secord%22">Adrian Secord</a>
 */

#ifndef LINALG_TYPE_HELPERS_H
#define LINALG_TYPE_HELPERS_H

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)

#include <cmath>

// Code to help handle the low-level details of the different types.
// Everything in this file deals with "plain old data" types.

namespace LinAlg {

    /// Helper class that defines functions that have different implementations
    /// depending on whether T is signed or not.
    /// The default implementation is for signed quantities.
    template <typename T, bool SignedType = std::numeric_limits<T>::is_signed>
    struct SignedHelper {
        /// Approximate equality
        static inline bool equal(const T first, const T second, const T eps) {
            //return std::abs(first - second) <= eps;
            return std::abs(first - second) <= eps * std::abs(first);
        }
    };
    
    /// Specialize for unsigned quantities.
    template <typename T>
    struct SignedHelper<T, false> {
        /// Approximate equality for unsigned types -- you cannot store 
        /// the difference of two unsigned numbers into the same unsigned number type.
        /// If we could ask for "the signed type that can store the difference between
        /// two numbers of type T" then we wouldn't have to do this.
        static inline bool equal(const T first, const T second, const T eps) {
            return (first <= second + eps) && (first + eps >= second);
        }
    };

    /// Return a string description of the parameter type.
    template <typename T> inline
    std::string type_description() {
        return "unknown";
    }

    /// Return a string description of the parameter type.
    template <> inline
    std::string type_description<float>() {
        return "float";
    }

    /// Return a string description of the parameter type.
    template <> inline
    std::string type_description<double>() {
        return "double";
    }

    /// Return a string description of the parameter type.
    template <> inline
    std::string type_description<int>() {
        return "int";
    }

    /// Return a string description of the parameter type.
    template <> inline
    std::string type_description<unsigned>() {
        return "unsigned";
    }
}

#endif
#endif