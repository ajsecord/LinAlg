/*
 * LinAlg: fixed-size vector and matrix library with LAPACK bindings.
 * Transform2Fwd.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_AFFINE_TRANSFORM2_FWD_H
#define LINALG_AFFINE_TRANSFORM2_FWD_H

namespace LinAlg {
    
    template <typename T> class Transform2;
    
    /// A 2D float transform 
    typedef Transform2<float> Transform2f;

    /// A 2D double transform 
    typedef Transform2<double> Transform2d;

    /// A 2D int transform 
    typedef Transform2<int> Transform2i;

    /// A 2D unsigned transform 
    typedef Transform2<unsigned> Transform2u;
}

#endif