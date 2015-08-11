/*
 * LinAlg: A fixed-length vector and matrix library.
 * Vec$D.cpp
 * Copyright 2005-2008 Adrian Secord.
 * Contact: <a href="http://www.google.com/search?q=%22Adrian%20Secord%22">Adrian Secord</a>
 */

#include <iostream>
#include <limits>
#include <sstream>

#include "TypeHelpers.h"
#include <LinAlg/Vec$D.h>

namespace LinAlg {
    
// Note: if you add other *free* functions (not member functions) here, make sure they are also added to 
// LINALG_INSTANTIATE_VECTOR macro at the end of this file.

template <typename T>
std::ostream& operator<<(std::ostream& o, const Vec$D<T>& v) {
    o << '(';
    
    for (unsigned i = 0; i < $D; ++i) 
        o << v(i) << (i < $D - 1 ? ", " : "");
    
    return o << ')';
}
    
template <typename T>
std::istream& operator>>(std::istream& in, Vec$D<T>& v) {
    char c[$D + 1];
    
    in >> c[0];    
    for (unsigned i = 0; i < $D; ++i) 
        in >> v(i) >> c[i+1];
    
    bool good = (c[0] == '(' && c[$D] == ')');
    for (unsigned i = 0; i < $D-1; ++i)
        good = good && (c[i+1] == ',');
    
    if (!good)
        in.setstate(std::ios::failbit);  
    
    return in;
}

template <typename T>
const Vec$D<T>& Vec$D<T>::zero() {
    static bool needInit = true;
    static Vec$D<T> zeroVec;
    if (needInit) {
        for (unsigned i = 0; i < $D; ++i)
            zeroVec(i) = 0;
        needInit = false;
    }
    return zeroVec;
}

template <typename T>
const Vec$D<T> Vec$D<T>::unit(const unsigned index) {
    assert(index < $D);
    Vec$D<T> v = zero();
    v.m_data[index] = 1;
    return v;
}

template <typename T>
bool Vec$D<T>::equal(const Vec$D<T>& v, const T eps) const {
    for (unsigned i = 0; i < $D; ++i)
        if (!SignedHelper<T>::equal(m_data[i], v(i), eps))
            return false;
    return true;
}

template <typename T>
std::string Vec$D<T>::description() {
    static std::string d;
    if (d.empty()) {
        std::ostringstream o;
        o << $D << "D " << type_description<T>() << " vector";
        d = o.str();
    }

    return d;
}
    
// Make the compiler generate code for the above template definitions so the 
// linker will have something to link to.  The alternative is to place these
// in the header, but then we have to pull in <ostream> and other headers 
// that are big and unwieldy and then your compile times suck.
#define LINALG_INSTANTIATE_VECTOR(T) \
template class Vec$D<T>; \
template std::ostream& operator<<(std::ostream& o, const Vec$D<T>& v); \
template std::istream& operator>>(std::istream& i, Vec$D<T>& v);

#if !defined(LINALG_SKIP_DEFAULT_INSTANTIATIONS) 
#   if !defined(DOXYGEN_SHOULD_SKIP_THIS)
        LINALG_INSTANTIATE_VECTOR(float);
        LINALG_INSTANTIATE_VECTOR(double);
        LINALG_INSTANTIATE_VECTOR(int);
        LINALG_INSTANTIATE_VECTOR(unsigned);
#   endif
#endif

#if defined(LINALG_INSTANTIATE_USER_TYPE)
    LINALG_INSTANTIATE_VECTOR(LINALG_INSTANTIATE_USER_TYPE);
#endif

#undef LINALG_INSTANTIATE_VECTOR
    
} // namespace
