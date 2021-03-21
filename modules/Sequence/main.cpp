// Copyright 2021 Smirnov Aleksandr
#include <gtest/gtest.h>

#include <vector>

#include "./func.h"
#include "./hansen_functions.h"

TEST(Sequential, Test_hfunc1) {
    double result = AGS(intervals[0], hfunc1);
    ASSERT_DOUBLE_EQ(round(result * 10) / 10, round(res[0][0] * 10) / 10);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
