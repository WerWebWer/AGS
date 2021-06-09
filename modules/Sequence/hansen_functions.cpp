// Copyright 2021 Smirnov Aleksandr
#include "../../modules/Sequence/hansen_functions.h"

double intervals[20][2] = { {-1.5, 11}, {2.7, 7.5}, {-10.0, 10.0}, {1.9, 3.9}, {0.0, 1.2},
{-10.0, 10.0}, {2.7, 7.5}, {-10.0, 10.0}, {3.1, 20.4}, {0.0, 10.0},
{-1.57, 6.28}, {0.0, 6.28}, {0.001, 0.99}, {0.0, 4.0}, {-5.0, 5.0},
{-3.0, 3.0}, {-4.0, 4.0}, {0.0, 6.0}, {0.0, 6.5}, {-10.0, 10.0} };

double hfunc1(double x) {
    return pow(x, 6) / 6.0 - 52.0 / 25.0 * pow(x, 5) + 39.0 / 80.0 * pow(x, 4) +
        71.0 / 10.0 * pow(x, 3) - 79.0 / 20.0 * pow(x, 2) - x + 0.1;
}

double hfunc2(double x) {
    return sin(x) + sin(10 * x / 3);
}

double hfunc3(double x) {
    double res = 0;
    for (int i = 1; i < 6; i++)
        res += i * sin((i + 1) * x + i);
    return -res;
}

double hfunc4(double x) {
    return (-16 * x * x + 24 * x - 5) * exp(-x);
}

double hfunc5(double x) {
    return -(-3 * x + 1.4) * sin(18 * x);
}

double hfunc6(double x) {
    return -(x + sin(x)) * exp(-x * x);
}

double hfunc7(double x) {
    return sin(x) + sin(10 * x / 3) + log(x) - 0.84 * x + 3;
}

double hfunc8(double x) {
    double res = 0;
    for (int i = 1; i < 6; i++)
        res += i * cos((i + 1) * x + i);
    return -res;
}

double hfunc9(double x) {
    return sin(x) + sin(2.0 / 3.0 * x);
}

double hfunc10(double x) {
    return -x * sin(x);
}

double hfunc11(double x) {
    return 2 * cos(x) + cos(2 * x);
}

double hfunc12(double x) {
    return pow(sin(x), 3) + pow(cos(x), 3);
}

double hfunc13(double x) {
    double sgn = 0.0;
    if (x * x - 1 < 0)
        sgn = -1.0;
    else
        sgn = 1.0;
    return -pow(x * x, 1.0 / 3.0) + sgn * pow(sgn * (x * x - 1.0), 1.0 / 3.0);
}

double hfunc14(double x) {
    return -exp(-x) * sin(2 * acos(-1.0) * x);
}

double hfunc15(double x) {
    return (x * x - 5 * x + 6) / (x * x + 1);
}

double hfunc16(double x) {
    return 2 * (x - 3) * (x - 3) + exp(x * x / 2);
}

double hfunc17(double x) {
    return pow(x, 6) - 15 * pow(x, 4) + 27 * x * x + 250;
}

double hfunc18(double x) {
    if (x <= 3)
        return (x - 2) * (x - 2);
    else
        return 2 * log(x - 2) + 1;
}

double hfunc19(double x) {
    return -x + sin(3 * x) - 1;
}

double hfunc20(double x) {
    return -(x - sin(x)) * exp(-x * x);
}

double hpfunc1(double x) {
    return pow(x, 5) - 10.4 * pow(x, 4) + 1.95 * pow(x, 3) + 21.3 * x * x -
        7.9 * x - 1.0;
}

double hpfunc2(double x) {
    return cos(x) + 10.0 * cos(10.0 * x / 3.0) / 3.0;
}

double hpfunc3(double x) {
    double res = 0.0;
    for (int i = 1; i < 6; i++)
        res += i * (i + 1) * cos((i + 1) * x + i);
    return -res;
}

double hpfunc4(double x) {
    return (16.0 * x * x - 56.0 * x + 29.0) * exp(-x);
}

double hpfunc5(double x) {
    return 3.0 * sin(18.0 * x) - 18.0 * (-3.0 * x + 1.4) * cos(18.0 * x);
}

double hpfunc6(double x) {
    return (2.0 * x * (x + sin(x)) - cos(x) - 1) * exp(-x * x);
}

double hpfunc7(double x) {
    return cos(x) + 10.0 * cos(10.0 * x / 3.0) / 3.0 + 1 / x - 0.84;
}

double hpfunc8(double x) {
    double res = 0.0;
    for (int i = 1; i < 6; i++)
        res += i * (i + 1) * sin((i + 1) * x + i);
    return res;
}

double hpfunc9(double x) {
    return cos(x) + 2.0 * cos(2.0 * x / 3.0) / 3.0;
}

double hpfunc10(double x) {
    return -sin(x) - x * cos(x);
}

double hpfunc11(double x) {
    return -2.0 * (sin(x) + sin(2.0 * x));
}

double hpfunc12(double x) {
    return 3.0 * cos(x) * sin(x) * (sin(x) - cos(x));
}

double hpfunc13(double x) {
    double st = (1.0 / 3.0);
    if (x == 0.0)
        return 0.0;
    return (2.0 * x / pow((x * x - 1) * (x * x - 1), st) - 2.0 * pow(x, -st)) / 3.0;
}

double hpfunc14(double x) {
    double pi = acos(-1.0);
    return exp(-x) * (sin(2.0 * pi * x) - 2.0 * pi * cos(2 * pi * x));
}

double hpfunc15(double x) {
    return (2.0 * x * (-x * x + 5.0 * x - 6.0) - (x * x + 1.0) * (5.0 - 2.0 * x))
        / ((x * x + 1.0) * (x * x + 1.0));
}

double hpfunc16(double x) {
    return 4.0 * (x - 3.0) + x * exp(x * x / 2.0);
}

double hpfunc17(double x) {
    return 6.0 * pow(x, 5) - 60.0 * x * x * x + 54.0 * x;
}

double hpfunc18(double x) {
    if (x <= 3)
        return 2.0 * x - 4.0;
    else
        return 2.0 / (x - 2.0);
}

double hpfunc19(double x) {
    return 3.0 * cos(3.0 * x) - 1.0;
}

double hpfunc20(double x) {
    return exp(-x * x) * (2.0 * x * (x - sin(x)) - 1.0 + cos(x));
}

double(*pfn[])(double x) = { hfunc1, hfunc2, hfunc3, hfunc4, hfunc5,
hfunc6, hfunc7, hfunc8, hfunc9, hfunc10,
hfunc11, hfunc12, hfunc13, hfunc14, hfunc15,
hfunc16, hfunc17, hfunc18, hfunc19, hfunc20 };

std::vector<std::vector<double>> res = { {10},{5.14575},{-0.49139, -6.77458, 5.79179},{2.868},{0.966086},
{0.679578},{5.199776},{-0.80032, -7.08351, 5.48286},{17.0392},{7.97867},
{2.09444, 4.18879},{4.71239, 3.14159},{0.70711},{0.22488},{2.41421},
{1.590721},{-3, 3},{2},{5.8728656},{1.195137} };
