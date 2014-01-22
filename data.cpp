/*
 *  Fengji Hou
 *  fh417@nyu.edu
 *  New York University
 *
 *  In data.cpp, Member functions of data class are defined. The functions
 *  include collecting data from data files and identify various sizes of the
 *  data.
 *
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <vector>
#include "data.h"

using namespace std;


// This routine counts how many rows there are in the data file, whether or not
// some of the rows are not data.
int64_t Data::count_rows() {
	string data_file_name("./data/" + data_name + ".dat"); 
	ifstream input(data_file_name.c_str());
	
	if(input.fail()) {
    throw DataException("Data constructor could not open file " + data_file_name + " !");
		return -1;
	}

	int64_t num_rows = 0;
	string temporary; // temporary string used to store the line read by getline()
	while ( !input.eof() ) {
		getline( input, temporary ); // get and count a line until the line is EOF
		if( temporary[0] != 0 ) ++ num_rows;
	}

	input.close();
	return num_rows;
}


// This function reads the data into 2-d vector 'data'.
void Data::collect_data() {
	string data_file_name("./data/" + data_name + ".dat"); 
	ifstream input(data_file_name.c_str());
 	
 	string temporary;
 	for (int64_t i = 0; i < skip_rows; ++i) {
 		getline( input, temporary );
 	}
 	
 	vector<double> input_array; // used to store the data in an array
 	double temp;
 	while(!input.eof()) {
 		input >> temp;
		input_array.push_back(temp);
	}
	
	num_col = input_array.size() / data_size; // determining number of columns
	data.resize(data_size, vector<double>(num_col, 0));
	for (int64_t i = 0; i < data_size; ++i) {
		for (int64_t j = 0; j < num_col; ++j) {
			data[i][j] = input_array[i*num_col + j];
		}
	}
  
	input.close();
}

Data::Data() {
}

// constructor
Data::Data(string name):          // name of the data
           data_name(name)
{
	skip_rows = 0;
	num_row   = count_rows();
	data_size = num_row - skip_rows;
	if (data_size > 0) {
		collect_data();
	}
	else {
		//throw DataException("Data constructor: data file is empty!");
	}
	
}

// constructor
Data::Data(string name,          // name of the data
           int64_t sr):          // how many rows to skip
           data_name(name),
           skip_rows(sr)
{
	num_row   = count_rows();
	data_size = num_row - skip_rows;
	if (data_size > 0) {
		collect_data();
	}
	else {
		//throw DataException("Data constructor: data file is empty!");
	}
	
}

// copy constructor
Data::Data(const Data& original):
           data_name(original.data_name),
           num_col(original.num_col),
           num_row(original.num_row),
           skip_rows(original.skip_rows),
           data_size(original.data_size)
{
	if (data_size > 0) {
		data.resize(data_size, vector<double>(num_col, 0));
		for (int64_t i = 0; i < data_size; ++i) {
			for (int64_t j = 0; j < num_col; ++j) {
				data[i][j] = original.data[i][j];
			}
		}
	}
}

Data& Data::operator=(const Data& original)
{
	data_name = original.data_name;
	num_col = original.num_col;
	num_row = original.num_row;
	skip_rows = original.skip_rows;
	data_size = original.data_size;
	if (data_size > 0) {
		data.resize(data_size, vector<double>(num_col, 0));
		for (int64_t i = 0; i < data_size; ++i) {
			for (int64_t j = 0; j < num_col; ++j) {
				data[i][j] = original.data[i][j];
			}
		}
	}
	return *this;
}

//    Code for throwing an exception with an error message and catching it outside

DataException::DataException(string ErrorMessage){  //  Set the exception message
	MessageText = ErrorMessage;
}

string DataException::ExceptionMessage(){
	return MessageText;
}
