#ifndef FUNK
#define FUNK

#include <cmath>
#include <algorithm>
#include <vector>
#include "./Pos.h"

double function(double x);
double M_m(std::vector<Pos>& a, double r);
std::pair<int, double> func_R(std::vector<Pos>& a, double m);

#endif  // FUNK