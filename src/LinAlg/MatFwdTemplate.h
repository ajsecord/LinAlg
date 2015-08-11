#ifndef LINALG_MAT$D_FWD_H
#define LINALG_MAT$D_FWD_H

namespace LinAlg {
    template <typename T> class Mat$D;
    
    /// A float matrix
    typedef Mat$D<float> Mat$Df;
    
    /// A double matrix
    typedef Mat$D<double> Mat$Dd;
    
    /// A int matrix
    typedef Mat$D<int> Mat$Di;
    
    /// An unsigned int matrix
    typedef Mat$D<unsigned> Mat$Du;
}

#endif
