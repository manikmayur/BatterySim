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
using namespace std;

bool sortbystrlength(string lhs, string rhs) {return lhs.length() > rhs.length();}

int main() {
	std::vector<std::string> strVec;
	std::string string_original;
	unsigned int pos;
	strVec.push_back("e_anode");
	strVec.push_back("e_anode_inactive");

	string cmpStr = "0.1*e_anode_inactive + e_anode_inactive";
	cout<<cmpStr<<endl;
	for(unsigned int n = 0; n < strVec.size(); n++)
	{
		string_original = strVec[n];
		//pos = cmpStr.find(string_original);
		if(cmpStr.find(string_original) != string::npos)
		{
			printf(
					"Before sorting: Found %s at n=%d\n",
					string_original.c_str(), n);
			break;
		}
	}
	sort(strVec.begin(), strVec.end(), sortbystrlength);
	for(unsigned int n = 0; n < strVec.size(); n++)
	{
		string_original = strVec[n];
		if(cmpStr.find(string_original) != string::npos)
		{
			printf(
					"After sorting: Found %s at n=%d\n",
					string_original.c_str(), n);
			break;
		}
	}

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
