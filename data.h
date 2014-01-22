/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *  In data.h, the class Data is declared.
 *
 */

#ifndef DATA_H
#define DATA_H

#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

class Data {
public:
	string data_name;  // name of the data file
	
	// The following vector are used to store data.
	vector< vector<double> > data;
  
  // The following variables are related to the data
	int64_t data_size;   // size of the data
	int64_t num_row;     // number of rows, (Not all rows are data.)
	int64_t num_col;     // number of columns
	int64_t skip_rows;   // specify how many rows to skip

	int64_t count_rows();
	void collect_data();
	
	// constructor
	Data();
	Data(string);                 // argument string for data file name
	Data(string, int64_t);        // argument string for data file name, long for skip_rows
	Data(const Data& original);   // copy constructor
	Data& operator=(const Data& original);  // assignment operator
	~Data() {};
};


class DataException {   // How to throw an exception if the constructor above fails
private: 
	string  MessageText;
public:
	DataException(string ErrorMessage);  //  Set the exception message 
	string ExceptionMessage();           //  Retrieve the exception message outside the constructor
};

#endif
