/*
 * LinAlg: fixed-size vector and matrix library with LAPACK bindings.
 * Transform3Fwd.h
 * Contact: http://www.google.com/search?q=%22Adrian+Secord%22
 * Copyright 2005-2008 Adrian Secord.
 */

#ifndef LINALG_AFFINE_TRANSFORM3_FWD_H
#define LINALG_AFFINE_TRANSFORM3_FWD_H

namespace LinAlg {
    
    template <typename T> class Transform3;
    
    /// A 3D float transform
    typedef Transform3<float> Transform3f;
    
    /// A 3D double transform
    typedef Transform3<double> Transform3d;
    /// A 3D int transform
    typedef Transform3<int> Transform3i;
    
    /// A 3D unsigned transform
    typedef Transform3<unsigned> Transform3u;
}

#endif


