//============================================================================
// Name        : Test.cpp
// Author      : Manik Mayur
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include "algorithms/algorithms.h"       // odeint function definitions


 double func(double t, double x)
 {
	 return exp(t);
 }


 // ------  Main
 int main()
 {
	 double res;
     double x0 = 10.0;
     // Integration parameters
     double t0 = 0.0;
     double t1 = 10.0;
     res = integrate(func, x0, t0, t1);

     //Print
     std::cout<<"out: "<<res<<"\n";

 }
/*int main() {
	std::string result;
	std::string res;
	std::string command = "D:\\HSO\\EESSoftwareLibrary\\McDENIS\\cantera\\ctml_writer.exe D:\\HSO\\EESSoftwareLibrary\\McDENIS\\cantera\\Final_KokamSiMET_2018_ModVal_KraftwerkBatterie.cti 2>&1 >&2";
	auto cti2xml_success = exec(command.c_str());
	std::cout <<"<>"<< result << std::endl;
//
	std::vector<std::string> strVec;
	std::string string_original;
	unsigned int pos;
	strVec.push_back("e_anode");
	strVec.push_back("e_anode_inactive");

	std::string cmpStr = "0.1*e_anode_inactive + e_anode_inactive";
	std::cout<<cmpStr<<std::endl;
	for(unsigned int n = 0; n < strVec.size(); n++)
	{
		string_original = strVec[n];
		if(cmpStr.find(string_original) != std::string::npos)
		{
			printf(
					"Before sorting: Found %s at n=%d\n",
					string_original.c_str(), n);
			break;
		}
	}
	std::sort(strVec.begin(), strVec.end(), sortbystrlength);
	for(unsigned int n = 0; n < strVec.size(); n++)
	{
		string_original = strVec[n];
		if(cmpStr.find(string_original) != std::string::npos)
		{
			printf(
					"After sorting: Found %s at n=%d\n",
					string_original.c_str(), n);
			break;
		}
	}
	std::cout << "!!!Hello World!!!" << std::endl; // prints !!!Hello World!!!
	return 0;
}*/
