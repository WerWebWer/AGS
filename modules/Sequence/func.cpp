// Copyright 2021 Smirnov Aleksandr
#include "stdio.h"
#include "math.h"
#include "time.h"

#include <iostream>
#include <algorithm>

#include "vector"
#include "iterator"

#include "../../modules/Sequence/func.h"

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

double AGS(double inter[2], double fun(double x)){
    double left = inter[0];
    double right = inter[1];
    std::cout << "[ " << left << " ; " << right << " ]" << std::endl;
    double r = 3; //
    double M = 0; //
    double e = 0.000001; //
    double new_point;
    double m; 
    int k = 2;
    double r_max = 0;
    std::vector <std::pair<double, double>> point;
    std::vector <double> R;
    std::pair<double, double> p1(left, fun(left));
    std::pair<double, double> p2(right, fun(right));
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
        point.push_back(std::pair<double, double>(new_point, fun(new_point)));
        sort(point.begin(), point.end());
        k++;
    }
    std::cout << point[r_max].first << " " << point[r_max].second << std::endl;
    return point[r_max].first;
}