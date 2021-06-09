// Copyright 2021 Smirnov Aleksandr
#ifndef MODULES_THREAD_OPS_STD_H_
#define MODULES_THREAD_OPS_STD_H_

#include <vector>

double ParallelGSA(double inter[2], double (*fun)(double x), double r, double e);
double ParallelOperations(double inter[2], double (*fun)(double x), double r, double e);
double ParallelNewPoints(double inter[2], double (*fun)(double x), double r, double e);

#endif  // MODULES_THREAD_OPS_STD_H_
