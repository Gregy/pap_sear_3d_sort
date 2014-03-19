/* 
 * File:   main.cpp
 * Author: gregy
 *
 * otazky: musi se lokalne seradit pred kazdym krokem shearsortu?
 */

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string.h>
#include <omp.h>
#include <sstream>

#define VTYPE int

using namespace std;


class CLIArgumentsParser {
public:
    static char* getCmdOption(char ** begin, char ** end, const std::string & option)
	{
		char ** itr = std::find(begin, end, option);
		if (itr != end && ++itr != end)
		{
			return *itr;
		}
		return 0;
	}

	static bool cmdOptionExists(char** begin, char** end, const std::string& option)
	{
		return std::find(begin, end, option) != end;
	}
};

class CPU {
	vector<VTYPE>::iterator startdata;
	vector<VTYPE>::iterator enddata;
	vector<VTYPE> temp;

public:
	CPU(vector<VTYPE>::iterator startdata, vector<VTYPE>::iterator enddata, size_t maxbucket):temp(2*maxbucket) {
		this->startdata = startdata;
		this->enddata = enddata;
	}
	void sort(bool reverse) {
		if(reverse)
			std::sort(this->startdata, this->enddata, std::greater<VTYPE>());
		else
			std::sort(this->startdata, this->enddata);
	}
	void mergesplit(CPU * other, bool reverse) {
		if(reverse)
			std::merge(this->startdata, this->enddata, other->startdata, other->enddata, this->temp.begin(), greater<VTYPE>());
		else
			std::merge(this->startdata, this->enddata, other->startdata, other->enddata, this->temp.begin());
		
		vector<VTYPE>::iterator tempit = this->temp.begin();
		for(vector<VTYPE>::iterator i = this->startdata; i != this->enddata; ++i,++tempit) {
			*i = *tempit;
		}
		for(vector<VTYPE>::iterator i = other->startdata; i != other->enddata; ++i,++tempit) {
			*i = *tempit;
		}
	}
	void printData(char delim = ' ', ostream & out = cout) {
		for(vector<VTYPE>::iterator i = this->startdata;i!=this->enddata; i++) {
			out << *i << delim;
		}
	}
	void printDataReverse(char delim = ' ', ostream & out = cout) {
		if(this->startdata == this->enddata)
			return;
		for(vector<VTYPE>::iterator i = this->enddata-1;i!=this->startdata; i--) {
			out << *i << delim;
		}
		out << *(this->startdata) << delim;
	}
};

class Coordinates {
public:
	int x;
	int y;
	int z;

	Coordinates(int x, int y, int z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	Coordinates& operator+=(const Coordinates& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}
	Coordinates& operator*=(const int& rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}
	Coordinates operator+(const Coordinates& rhs) const{
		Coordinates n = (*this);
		n += rhs;
	    return n;
	}
	Coordinates operator*(const int& rhs) const{
		Coordinates n = (*this);
		n *= rhs;
	    return n;
	}
};

static const Coordinates XVector(1,0,0);
static const Coordinates YVector(0,1,0);
static const Coordinates ZVector(0,0,1);

class CPUMatrix {
	size_t x,y,z = 0;
    vector<CPU*> CPUs;
	size_t bucket = 0;

public:
	CPUMatrix(char * dimensions, vector<VTYPE> &numbers) {
		stringstream ss;
		ss<< dimensions;
		int x = 1,y = 1,z = 1;
		ss >> x;
		if(ss.peek() == 'x') {
			ss.ignore();
		}
		else {
			throw "spatny format dimenzi mantice ma byt NxN[xN]";
		}
		ss >> y;
		if(ss.good()) {
			if(ss.peek() == 'x') {
				ss.ignore();
			}
			else {
				throw "spatny format dimenzi mantice ma byt NxN[xN]";
			}
			ss >> z;
		}
		cout << "Dimenze jsou " << x << "x" << y <<"x"<<z<<endl;
        if(x<1 || y<1 || z<1) {
            throw "Dimenze musi byt kladne!";
        }
        this->CPUs.resize(x*y*z);
		init(x,y,z, numbers);
	}
	CPUMatrix(int x, int y, int z, vector<VTYPE> &numbers): CPUs(x*y*z) {
       init(x,y,z,numbers); 
    }
	void init(int x, int y, int z, vector<VTYPE> &numbers) {
		this->x = x;
		this->y=y;
		this->z = z;
		
		if(numbers.size() < CPUs.size()) {
			cout << "Grid procesoru je vetsi nez pocet cisel...fail" << endl;
			throw "grid procesoru je vetsi nez pocet cisel";
		}
		size_t bucket = numbers.size()/this->CPUs.size()+ (numbers.size() % this->CPUs.size() != 0);
		this->bucket = bucket;

		int pos = 0;
		for(vector<CPU*>::iterator i = this->CPUs.begin(); i!=this->CPUs.end(); ++i,++pos) {
			//handle crazy case where there are cpus with no numbers
			if(pos*bucket > numbers.size()) {
				cout << "Warning, cpus without work..." << endl;
				*i = new CPU(numbers.end(), numbers.end(), bucket);
			}
			else if((pos+1)*bucket > numbers.size()) {
				*i = new CPU(numbers.begin()+pos*bucket, numbers.end(), bucket);
			}
			else {
				*i = new CPU(numbers.begin()+pos*bucket, numbers.begin()+(pos+1)*bucket, bucket);
			}
		}
	}
    CPU*& get(size_t x, size_t y, size_t z) {
		if(!this->coorValid(x,y,z)) {
			throw "mimo rozmer procesoru";
		}
        return CPUs.at(x + y * this->x + z * this->x * this->y);
    }
    CPU*& get(Coordinates c) {
        return this->get(c.x,c.y,c.z);
    }
	size_t getBucketSize() {
		return this->bucket;
	}
	bool coorValid(size_t x, size_t y, size_t z) {
		if(x>=this->x || y>=this->y || z>=this->z) {
			return false;
		}
		return true;
	}
	CPU*& get(size_t i) {
		return CPUs.at(i);
	}
	size_t size(){
		return this->CPUs.size();
	}
	size_t sizeForVector(Coordinates vector) {
		return this->x*vector.x+this->y*vector.y+this->z*vector.z;
	}
	void printMatrix() {
		for(size_t z=0;z<this->z; z++) {
			cout << "-----------------------" << endl;
			for(int y = this->y-1;y>=0;y--) {
				cout << endl;
				for(size_t x=0;x<this->x;x++) {
					cout << "|";
					this->get(x,y,z)->printData();
					cout << "|";
				}
				cout << endl;
			}
			cout << "-----------------------" << endl;
		}
	}
};
size_t logbin(size_t number) {
	size_t ret = 0;
	while (number >>= 1) { ++ret; }
	return ret;
}
void eotSortOpenMP(Coordinates line, Coordinates dynvector, CPUMatrix * matrix, bool reverse = false, bool justOneRun = false) {

	size_t dimension_size = matrix->sizeForVector(dynvector);
	//locally sort
	#pragma omp parallel num_threads(dimension_size)
	{
		#pragma omp for
		for(size_t i=0; i<dimension_size; i++) {
			matrix->get(line+dynvector*i)->sort(reverse);
		}
	}
    size_t runtimes;
    if(justOneRun) {
        runtimes = 1;
    }
    else {
		runtimes = dimension_size/2+1;
    }
	//TODO: prozkoumat to +1
	for(size_t count=0;count<runtimes;count++) {
		#pragma omp parallel num_threads(dimension_size/2)
		{
			//even pairs
			#pragma omp for
			for(size_t i =0; i<dimension_size-1; i+=2) {
				matrix->get(line+dynvector*i)->mergesplit(matrix->get(line+dynvector*(i+1)), reverse);
			}
			#pragma omp barrier
			//odd pairs
			#pragma omp for
			for(size_t i=1; i<dimension_size-1; i+=2) {
				matrix->get(line+dynvector*i)->mergesplit(matrix->get(line+dynvector*(i+1)), reverse);
			}
		}
	}
}
void shearSortOpenMP(Coordinates plane, Coordinates rowvector, Coordinates colvector, CPUMatrix * matrix, bool reverse = false) {

	size_t rowsize = matrix->sizeForVector(rowvector);
	size_t colsize = matrix->sizeForVector(colvector);
	for(size_t count=0;count<logbin(colsize)+2;count++) {

		//faze schvalne opacne, posledni musi vzdy byt snake faze
		#pragma omp parallel num_threads(rowsize)
		{
			//down phase
			#pragma omp for
			for(size_t i=0; i<rowsize; i++) {
				eotSortOpenMP(plane+rowvector*i, colvector, matrix, reverse);
			}
		}

		#pragma omp parallel num_threads(colsize)
		{
			//snake phase
			#pragma omp for
			for(size_t i =0; i<colsize; i++) {
				if(i%2==0)
					eotSortOpenMP(plane+colvector*i, rowvector, matrix, reverse);
				else
					eotSortOpenMP(plane+colvector*i, rowvector, matrix, !reverse);
			}
		}
	}
}
void dSortOpenMP(CPUMatrix * matrix, bool reverse = false) {

	size_t xsize = matrix->sizeForVector(XVector);
	size_t ysize = matrix->sizeForVector(YVector);
	size_t zsize = matrix->sizeForVector(ZVector);
    Coordinates zero(0,0,0);
	#pragma omp parallel num_threads(zsize)
    {
		#pragma omp for
		for(size_t xyplane=0;xyplane<zsize;xyplane++) {
			shearSortOpenMP(zero+ZVector*xyplane, XVector, YVector, matrix, reverse);
		}
    }
	#pragma omp parallel num_threads(xsize)
    {
		#pragma omp for
		for(size_t yzplane=0;yzplane<xsize;yzplane++) {
			shearSortOpenMP(zero+XVector*yzplane, ZVector, YVector, matrix, reverse);
		}
    }
	#pragma omp parallel num_threads(ysize)
    {
		#pragma omp for
		for(size_t xzplane=0;xzplane<ysize;xzplane++) {
            if(xzplane%2 == 0) {
				shearSortOpenMP(zero+YVector*xzplane, XVector, ZVector, matrix, reverse);
            }
			else {
				shearSortOpenMP(zero+YVector*xzplane, XVector, ZVector, matrix, !reverse);
            }
		}
    }
	#pragma omp parallel num_threads(xsize*zsize)
    {
		#pragma omp for
		for(size_t column=0;column<xsize*zsize;column++) {
			eotSortOpenMP(zero+XVector*(column%xsize)+ZVector*(column/xsize), YVector, matrix,reverse,true);
		}
    }
	#pragma omp parallel num_threads(ysize)
    {
		#pragma omp for
		for(size_t xzplane=0;xzplane<ysize;xzplane++) {
			shearSortOpenMP(zero+YVector*xzplane, XVector, ZVector, matrix, reverse);
		}
    }

}

void readFile(char * filename, vector<VTYPE> &data) {
	ifstream inputFile(filename, std::ifstream::in);

	if (inputFile) {        
    	VTYPE value;

		while ( inputFile >> value ) {
			data.push_back(value);
		}
        if(!inputFile.eof()) {
			throw "chyba nacitani dat ze souboru - fakt tam jsou jen cisla?";
        }
		inputFile.close();
    }
    else {
        throw "chyba nacitani dat ze souboru";
    }
}

void writeFile(char * filename, vector<VTYPE>& data) {
	streambuf * buf;
	ofstream of;

	if(filename) {
		of.open(filename);
		buf = of.rdbuf();
	} else {
		buf = std::cout.rdbuf();
	}

	ostream out(buf);
	
	if (out) {        
		for(vector<VTYPE>::iterator i = data.begin(); i!=data.end(); i++) {
			out << *i << endl;
		}
    }
    else {
        throw "chyba vypisu dat";
    }
}
void getOutStream(char * filename, ostream &ost){
}
void writeFileSnake(ostream &out, Coordinates plane, Coordinates rowvector, Coordinates colvector, CPUMatrix * matrix) {
	
	if (out) {
		size_t rows = matrix->sizeForVector(colvector);
		size_t cols = matrix->sizeForVector(rowvector);
		for(size_t row =0;row<rows;row++) {
			if(row%2==0) {
				for(size_t col = 0; col<cols;col++) {
					matrix->get(plane+rowvector*col+colvector*row)->printData('\n', out);
				}
			}
			else {
				for(size_t col = cols-1; col>0;col--) {
					matrix->get(plane+rowvector*col+colvector*row)->printDataReverse('\n', out);
				}
				matrix->get(plane+colvector*row)->printDataReverse('\n', out);
			}
		}

    }
    else {
        throw "chyba vypisu dat";
    }
}

int main(int argc, char** argv) {
	omp_set_nested(1);

	try {
		char * filename = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-f");
		if (!filename) {
			cout << "Chybi argument -f se souborem cisel" << endl;
			return -1;
		}
		char * output = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-o");
		streambuf * buf;
		ofstream of;

		if(output) {
			of.open(output);
			buf = of.rdbuf();
		} else {
			buf = std::cout.rdbuf();
		}

		ostream outputStream(buf);


		char * dimensions = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-S");
		if(dimensions) {
			//SHEAR
			vector<VTYPE> data;
			readFile(filename, data);

			CPUMatrix matrix(dimensions, data);
            if(matrix.sizeForVector(ZVector) != 1) {
                cout << "Shear funguje jen na 2d mrizce..." << endl;
                return -1;
            }
			shearSortOpenMP(Coordinates(0,0,0), XVector, YVector,&matrix);
			writeFileSnake(outputStream, Coordinates(0,0,0), XVector,YVector, &matrix);
			return 0;
		}

		
		dimensions = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-3");
		if(dimensions) {
			//3D
			vector<VTYPE> data;
			readFile(filename, data);


			CPUMatrix matrix(dimensions, data);
            dSortOpenMP(&matrix, false);
			for(size_t xzplanes=0;xzplanes<matrix.sizeForVector(YVector); xzplanes++) {
				writeFileSnake(outputStream, Coordinates(0,0,0)+YVector*xzplanes, XVector,ZVector, &matrix);
            }
			return 0;
		}

		cout << "Nespecifikovan ani shear (-S) ani 3dsort (-3). Syntax [-S|-3] NxN[xN] -- rozmer matice procesoru" << endl;
		return -1;

	}
	catch (char const * e) {
		cout << e << endl;
		return -1;
	}
} 