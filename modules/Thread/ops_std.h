// Copyright 2021 Smirnov Aleksandr
#ifndef MODULES_THREAD_OPS_STD_H_
#define MODULES_THREAD_OPS_STD_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int  sz);
int getParallelOperations(std::vector<int> vec, const std::string& ops);
int getSequentialOperations(std::vector<int> vec, const std::string& ops);

#endif  // MODULES_THREAD_OPS_STD_H_
