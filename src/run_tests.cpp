#include "TestVector.h"
#include "TestMatrix.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

#include <LinAlg/Vec2.h>
#include <LinAlg/Vec3.h>
#include <LinAlg/Vec4.h>

#include <LinAlg/Mat2.h>
#include <LinAlg/Mat3.h>
#include <LinAlg/Mat4.h>

#define TEST_VECTOR(VecType, trials) { \
    std::cout << "Testing " << VecType::description() << std::endl; \
    for (int i = 0; i < trials; ++i) \
        LinAlg::Test::Vector::Tester<VecType>::test(); \
}

#define TEST_MATRIX(MatType, trials) { \
    std::cout << "Testing " << MatType::description() << std::endl; \
    for (int i = 0; i < trials; ++i) \
        LinAlg::Test::Matrix::Tester<MatType>::test(); \
}

int main(int argc, char* argv[]) {
    using namespace LinAlg;
    
    //const unsigned seed = (unsigned)(argc > 1 ? atol(argv[1]) : std::time(NULL));
    const unsigned seed = 1211904253;
    std::cout << "Random seed: " << seed << std::endl;
    srand(seed);
    
    const int trials = (argc > 2 ? atol(argv[2]) : 100);
    std::cout << "Number of trials: " << trials << std::endl;
    
    try {
        TEST_VECTOR(Vec2f, trials);
        TEST_VECTOR(Vec2d, trials);
        TEST_VECTOR(Vec2i, trials);
        TEST_VECTOR(Vec2u, trials);
        
        TEST_VECTOR(Vec3f, trials);
        TEST_VECTOR(Vec3d, trials);
        TEST_VECTOR(Vec3i, trials);
        TEST_VECTOR(Vec3u, trials);
        
        TEST_VECTOR(Vec4f, trials);
        TEST_VECTOR(Vec4d, trials);
        TEST_VECTOR(Vec4i, trials);
        TEST_VECTOR(Vec4u, trials);
        
    } catch (std::runtime_error& e) {
        std::cout << "FAIL: " << e.what() << std::endl;
        return -1;
    }

    try {
        TEST_MATRIX(Mat2f, trials);
        TEST_MATRIX(Mat2d, trials);
        TEST_MATRIX(Mat2i, trials);
        TEST_MATRIX(Mat2u, trials);
        
        TEST_MATRIX(Mat3f, trials);
        TEST_MATRIX(Mat3d, trials);
        TEST_MATRIX(Mat3i, trials);
        TEST_MATRIX(Mat3u, trials);
        
        TEST_MATRIX(Mat4f, trials);
        TEST_MATRIX(Mat4d, trials);
        TEST_MATRIX(Mat4i, trials);
        TEST_MATRIX(Mat4u, trials);
        
    } catch (std::runtime_error& e) {
        std::cout << "FAIL: " << e.what() << std::endl;
        return -1;
    }
    
    std::cout << "All tests passed with " << trials << " sets of random data." << std::endl;
    return 0;
}