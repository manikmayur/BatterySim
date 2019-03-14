/*
 * Utilities.h
 *
 *  Created on: 07.02.2019
 *      Author: Manik
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

// Sort string by length
static bool sortbystrlength(std::string lhs, std::string rhs) {return lhs.length() > rhs.length();}

// Get message from stderr to log file
static std::string execCMD(const char* cmd);



#endif /* UTILITIES_H_ */
