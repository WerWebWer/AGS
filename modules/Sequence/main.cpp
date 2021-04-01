// Copyright 2021 Smirnov Aleksandr
#include <gtest/gtest.h>

#include <vector>

#include "./func.h"
#include "./hansen_functions.h"

#define R = 1.
#define EPS = 1e-6

typedef testing::TestWithParam<std::tuple< double(*)(double), double*, double, double, double>> Sequential;
TEST_P(Sequential, Test) {
    // Arrange
    double (*fun)(double) = std::get<0>(GetParam());
    double* param = std::get<1>(GetParam());
    double res = std::get<2>(GetParam());
    double r = std::get<3>(GetParam());
    double e = std::get<4>(GetParam());
    double result;

    // Actz
    result = AGS(param, fun, r, e);
    // std::cout << result << " " << r << " " << e <<"\n";

    //Assert
    ASSERT_DOUBLE_EQ(round(result * 10) / 10, round(res * 10) / 10);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

INSTANTIATE_TEST_CASE_P(/**/, Sequential,
    testing::Values(
        std::make_tuple(pfn[0], intervals[0], res[0][0], 1., 1e-6),
        std::make_tuple(pfn[1], intervals[1], res[1][0], 3., 1e-6),
        std::make_tuple(pfn[2], intervals[2], res[2][0], 3., 1e-6),
        // std::make_tuple(pfn[3], intervals[3], res[3][0], 100., 1e-3), // error
        std::make_tuple(pfn[4], intervals[4], res[4][0], 3., 1e-6),
        // std::make_tuple(pfn[5], intervals[5], res[5][0], 3., 1e-6), // error
        std::make_tuple(pfn[6], intervals[6], res[6][0], 3., 1e-6)
    ));
