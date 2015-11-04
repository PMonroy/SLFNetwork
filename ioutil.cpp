#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream> // std::cout, std::endl
#include <fstream>
#include <sstream>
#include <iomanip> // std::setfill, std::setw

using namespace std;

string DoubletoString(int ndigits, int ndecimals, double number)
{
  ostringstream stream;// it needs to include  <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>

  stream << fixed; //it needs to include <iostream>

  stream << setfill('0') << setw(ndigits);

  stream << setprecision(0) << (abs(number)*factor);
  
  return stream.str();
}


// This function counts the number of lines of a file
int CountLines(const char *filename)
  {
    ifstream ifile(filename);
    int numlines;
    string line;
    
    if (!ifile.is_open()) 
      {
	string me = "countlines()";
	cout << me <<": Skipping unreadable file \"" << filename << "\" "<<endl;
	return -1;
      }
    for ( numlines=0; getline(ifile, line); ++numlines)
      ;

    ifile.close();
    return numlines;
  }

bool FileExist(const char *filename)
{
    ifstream ifile(filename);
    return ifile.good();
}
