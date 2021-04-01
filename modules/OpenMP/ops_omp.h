// Copyright 2021 Smirnov Aleksandr
#ifndef MODULES_OPENMP_OPS_OMP_H_
#define MODULES_OPENMP_OPS_OMP_H_

#include <vector>
#include <string>

double ParallelOperations(double inter[2], double fun(double x), double r, double e);
double ParallelAGS(double inter[2], double fun(double x), double r, double e);  // internals parallelization 

double Rmax(std::vector<double> a);
double AGS(double inter[2], double fun(double x), double r, double e);

#endif  // MODULES_OPENMP_OPS_OMP_H_
