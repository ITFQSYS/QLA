#include <gtest/gtest.h>

#include"test_linalg.hpp"
#include"test_Matrix.hpp"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}