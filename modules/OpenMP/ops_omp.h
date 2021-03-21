// Copyright 2021 Smirnov Aleksandr
#ifndef MODULES_OPENMP_OPS_OMP_H_
#define MODULES_OPENMP_OPS_OMP_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
int getParallelOperations(std::vector<int> vec, const std::string& ops);
int getSequentialOperations(std::vector<int> vec, const std::string& ops);

#endif  // MODULES_OPENMP_OPS_OMP_H_
