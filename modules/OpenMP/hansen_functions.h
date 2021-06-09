// Copyright 2021 Smirnov Aleksandr
#ifndef MODULES_SEQUENCE_HANSEN_FUNCTION_H_
#define MODULES_SEQUENCE_HANSEN_FUNCTION_H_

#include <cmath>
#include <vector>

extern  double intervals[20][2];

double hfunc1(double x);
double hfunc2(double x);
double hfunc3(double x);
double hfunc4(double x);
double hfunc5(double x);
double hfunc6(double x);
double hfunc7(double x);
double hfunc8(double x);
double hfunc9(double x);
double hfunc10(double x);
double hfunc11(double x);
double hfunc12(double x);
double hfunc13(double x);
double hfunc14(double x);
double hfunc15(double x);
double hfunc16(double x);
double hfunc17(double x);
double hfunc18(double x);
double hfunc19(double x);
double hfunc20(double x);
double hpfunc1(double x);
double hpfunc2(double x);
double hpfunc3(double x);
double hpfunc4(double x);
double hpfunc5(double x);
double hpfunc6(double x);
double hpfunc7(double x);
double hpfunc8(double x);
double hpfunc9(double x);
double hpfunc10(double x);
double hpfunc11(double x);
double hpfunc12(double x);
double hpfunc13(double x);
double hpfunc14(double x);
double hpfunc15(double x);
double hpfunc16(double x);
double hpfunc17(double x);
double hpfunc18(double x);
double hpfunc19(double x);
double hpfunc20(double x);

extern double(*pfn[])(double x);

extern std::vector<std::vector<double>> res;


#endif  // MODULES_SEQUENCE_HANSEN_FUNCTION_H_
