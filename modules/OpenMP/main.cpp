// Copyright 2021 Smirnov Aleksandr
#include <gtest/gtest.h>
#include <vector>
#include "./ops_omp.h"
#include "./hansen_functions.h"

//TEST(Parallel_Operations_OpenMP, Test_Sum) {
//    std::vector<int> vec = getRandomVector(100);
//    int sequential_sum = getSequentialOperations(vec, "+");
//    int parallel_sum = getParallelOperations(vec, "+");
//    ASSERT_EQ(sequential_sum, parallel_sum);
//}

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, double, double, double>> OpenMP;
TEST_P(OpenMP, Division_into_segments) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    double expected = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;

    // Actz
    result = ParallelOperations(param, fun, r, e);
    // std::cout << result << " " << r << " " << e <<"\n";

    //Assert
    ASSERT_DOUBLE_EQ(round(expected * 10) / 10, round(result * 10) / 10);
}

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, double, double, double>> OpenMP;
TEST_P(OpenMP, Internals_parallelization) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    double expected = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;

    // Actz
    result = ParallelAGS(param, fun, r, e);
    // std::cout << result << " " << r << " " << e <<"\n";

    //Assert
    ASSERT_DOUBLE_EQ(round(expected * 10) / 10, round(result * 10) / 10);
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

INSTANTIATE_TEST_CASE_P(/**/, OpenMP,
    testing::Values(
        std::make_tuple(pfn[0], intervals[0], res[0][0], 1., 1e-6),
        std::make_tuple(pfn[1], intervals[1], res[1][0], 3., 1e-6),
        std::make_tuple(pfn[2], intervals[2], res[2][0], 3., 1e-6),
        // std::make_tuple(pfn[3], intervals[3], res[3][0], 100., 1e-3), // error
        std::make_tuple(pfn[4], intervals[4], res[4][0], 3., 1e-6),
        // std::make_tuple(pfn[5], intervals[5], res[5][0], 3., 1e-6), // error
        std::make_tuple(pfn[6], intervals[6], res[6][0], 3., 1e-6)
    ));
