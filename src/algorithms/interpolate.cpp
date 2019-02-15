/*
 * interpolate.cpp
 *
 *  Created on: 08.02.2019
 *      Author: Manik
 */

#include "algorithms/algorithms.h"

double interpolate(double x, std::vector<double> xList, std::vector<double> fList)
{
    double c;
    // Assumes that "table" is sorted by .first
    // Check if x is out of bound
    if (x > xList.back()) {
        c = fList.back();
        return c;
    }
    if (x < xList[0]) {
        c = fList[0];
        return c;
    }
    size_t i = std::distance(xList.begin(), std::lower_bound(xList.begin(), xList.end(), x));
    c = fList[i-1] + (fList[i] - fList[i-1]) * (x - xList[i-1])/(xList[i]- xList[i-1]);
    return c;
}

