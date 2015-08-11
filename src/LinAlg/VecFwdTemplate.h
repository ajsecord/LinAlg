#ifndef LINALG_VEC$D_FWD_H
#define LINALG_VEC$D_FWD_H

namespace LinAlg {
    template <typename T> class Vec$D;

    /// A float vector
    typedef Vec$D<float> Vec$Df;
    
    /// A double vector
    typedef Vec$D<double> Vec$Dd;
    
    /// A int vector
    typedef Vec$D<int> Vec$Di;
    
    /// An unsigned int vector
    typedef Vec$D<unsigned> Vec$Du;
}

#endif
