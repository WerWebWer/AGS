#ifndef FUNK
#define FUNK
#include <iostream>
#include <algorithm>
#include "vector"
#include "iterator"
#include "stdio.h"
#include "math.h"
#include "time.h"

double function(double x) {
    return sin(x);
};

double Rmax(std::vector<double> a){
    double R = a[0];
    int count = 0;
    for (int i = 1; i < a.size(); i++) {
        if (a[i] > R) {
            R = a[i];
            count = i;
        }
    }
    return count;
};

double AGS(double left, double right){
    double r = 3; //
    double M = 0; //
    double e = 0.02; //
    double new_point;
    double m; 
    int k = 2;
    double r_max = 0;
    std::vector <std::pair<double, double>> point;
    std::vector <double> R;
    std::pair<double, double> p1(left, function(left));
    std::pair<double, double> p2(right, function(right));
    point.insert(point.begin(), p1);
    point.push_back(p2);
    while ((point[r_max + 1].first - point[r_max].first) > e) {
        for (int i = 0; i < k - 1; i++) {
            int tmpM = abs((point[i + 1].second - point[i].second) / (point[i + 1].first - point[i].first));
            if (M < tmpM) M = tmpM;
        }
        m = (M == 0) ? 1 : r * M;
        for (int i = 0; i < k - 1; i++) {
            R.push_back(m * (point[i + 1].first - point[i].first) + ((point[i + 1].second - point[i].second) * (point[i + 1].second - point[i].second)) / (m * (point[i + 1].first - point[i].first)) - 2 * (point[i + 1].second + point[i].second));
        }
        r_max = Rmax(R);
        R.clear();
        new_point = (point[r_max + 1].first + point[r_max].first) / 2 - (point[r_max + 1].second - point[r_max].second) / (2 * m);
        point.push_back(std::pair<double, double>(new_point, function(new_point)));
        sort(point.begin(), point.end());
        k++;
    }
    return point[r_max].second;
}

#endif  // FUNK