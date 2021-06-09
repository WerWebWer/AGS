// Copyright 2021 Smirnov Aleksandr
#include <gtest/gtest.h>

#include <vector>

#include "./ops_omp.h"
#include "./hansen_functions.h"
#include <omp.h>

::testing::AssertionResult resultInExpected(std::vector<double> expected, double result) {
    bool flag = false;
    for (size_t i = 0; i < expected.size(); i++) {
        if (round(expected[i] * 10) / 10 == round(result * 10) / 10) {
            return ::testing::AssertionSuccess();
        }
    }
    std::cout << "Actual:\n    " << result << "\nExpected:\n    { ";
    for (double i : expected) std::cout << i << " ";
    std::cout << "}" << std::endl;
    return ::testing::AssertionFailure();

}

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, std::vector<double>, double, double>> OpenMP;
TEST_P(OpenMP, Division_into_segments) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    std::vector<double> expected = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;
    double start = clock();

    // Actz
    result = ParallelGSA(param, fun, r, e);

    printf(">>> TIME = %.4lf sec\n", (clock() - start) / CLOCKS_PER_SEC);

    //Assert
    ASSERT_TRUE(resultInExpected(expected, result));
}

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, std::vector<double>, double, double>> OpenMP;
TEST_P(OpenMP, Internals_parallelization) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    std::vector<double> expected = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;
    double start = clock();

    // Actz
    result = ParallelOperations(param, fun, r, e);

    printf(">>> TIME = %.4lf sec\n", (clock() - start) / CLOCKS_PER_SEC);

    //Assert
    ASSERT_TRUE(resultInExpected(expected, result));
}

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, std::vector<double>, double, double>> OpenMP;
TEST_P(OpenMP, Parallel_search_for_new_points) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    std::vector<double> expected = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;
    double start = clock();

    // Actz
    result = ParallelNewPoints(param, fun, r, e);

    printf(">>> TIME = %.4lf sec\n", (clock() - start) / CLOCKS_PER_SEC);

    //Assert
    ASSERT_TRUE(resultInExpected(expected, result));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

INSTANTIATE_TEST_CASE_P(/**/, OpenMP, testing::Values(
    std::make_tuple(pfn[0], intervals[0], res[0], 2., 1e-6),
    std::make_tuple(pfn[1], intervals[1], res[1], 2., 1e-6),
    std::make_tuple(pfn[2], intervals[2], res[2], 2., 1e-6),
    std::make_tuple(pfn[3], intervals[3], res[3], 2., 1e-6),
    std::make_tuple(pfn[4], intervals[4], res[4], 2., 1e-6),
    std::make_tuple(pfn[5], intervals[5], res[5], 2., 1e-6),
    std::make_tuple(pfn[6], intervals[6], res[6], 2., 1e-6),
    std::make_tuple(pfn[7], intervals[7], res[7], 2., 1e-6),
    std::make_tuple(pfn[8], intervals[8], res[8], 2., 1e-6),
    std::make_tuple(pfn[9], intervals[9], res[9], 2., 1e-6),
    std::make_tuple(pfn[10], intervals[10], res[10], 2., 1e-6),
    std::make_tuple(pfn[11], intervals[11], res[11], 2., 1e-6),
    std::make_tuple(pfn[12], intervals[12], res[12], 2., 1e-6),
    std::make_tuple(pfn[13], intervals[13], res[13], 2., 1e-6),
    std::make_tuple(pfn[14], intervals[14], res[14], 2., 1e-6),
    std::make_tuple(pfn[15], intervals[15], res[15], 2., 1e-6),
    std::make_tuple(pfn[16], intervals[16], res[16], 2., 1e-6),
    std::make_tuple(pfn[17], intervals[17], res[17], 2., 1e-6),
    std::make_tuple(pfn[18], intervals[18], res[18], 2., 1e-6),
    std::make_tuple(pfn[19], intervals[19], res[19], 2., 1e-6)
));
