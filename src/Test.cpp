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
using namespace std;

// Sort string by length
bool sortbystrlength(string lhs, string rhs) {return lhs.length() > rhs.length();}
// Get message from stderr to log file
std::string exec(const char* cmd);

std::string exec(const char* cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    auto pipe = popen(cmd, "r");

    if (!pipe) throw std::runtime_error("popen() failed!");

    while (!feof(pipe))
    {
        if (fgets(buffer.data(), 128, pipe) != nullptr)
            result += buffer.data();
    }

    auto rc = pclose(pipe);

    if (rc == EXIT_SUCCESS)
    {
        std::cout << "SUCCESS\n";
    }
    else
    {
        std::cout << "FAILED\n";
    }

    return result;
}

int main() {
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

	string cmpStr = "0.1*e_anode_inactive + e_anode_inactive";
	cout<<cmpStr<<endl;
	for(unsigned int n = 0; n < strVec.size(); n++)
	{
		string_original = strVec[n];
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

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!*/
	return 0;
}
