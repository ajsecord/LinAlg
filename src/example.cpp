/* A simple example to show the basic idea of the library. */

#include <LinAlg/Vec3.h>
#include <LinAlg/Mat3.h>
#include <iostream>

int main(int argc, char* argv[]) {
    using namespace LinAlg;
    
    // Create a vector
    Vec3f v(1, 2, 3);
    std::cout << v << std::endl;
    
    // Create a matrix
    Mat3f m = Mat3f::ident();
    m(0,0) = 1;
    m(1,1) = 2;
    m(2,2) = 3;

    // Matrix-vector multiply
    Vec3f result = m * v;
    
    std::cout << "Multiply " << v << " by \n" << m << "\nand you get " << result << '\n';
    
    return 0;
}