/* 
 * File:   main.cpp
 * Author: gregy
 *
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
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <sys/time.h>
#define VTYPE int

using namespace std;


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

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
public:
	unsigned int startindex;
	unsigned int endindex;
	unsigned int tempindex;

	CPU(unsigned int startindex=0, unsigned int endindex=0, unsigned int tempindex = 0) {
		this->startindex = startindex;
		this->endindex = endindex;
		this->tempindex = tempindex;
	}
	void printData( int * data,char delim = ' ', ostream & out = cout) {
		for(unsigned int i=startindex; i< endindex; i++) {
			out <<  data[i] << delim;
		}
	}
	void printDataReverse(int * data, char delim = ' ', ostream & out = cout) {
		for(unsigned int i = endindex-1;i>=this->startindex && i<endindex; i--) {
			out << data[i] << delim;
		}
	}
};

struct Coordinates {
	int x;
	int y;
	int z;

	__host__ __device__ Coordinates(int x, int y, int z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
	__host__ __device__ Coordinates& operator+=(const Coordinates& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}
	 __host__ __device__ Coordinates& operator*=(const int& rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}
	__host__ __device__ Coordinates operator+(const Coordinates& rhs) const{
		Coordinates n = (*this);
		n += rhs;
	    return n;
	}
	__host__ __device__ Coordinates operator*(const int& rhs) const{
		Coordinates n = (*this);
		n *= rhs;
	    return n;
	}
};

static const Coordinates XVector(1,0,0);
static const Coordinates YVector(0,1,0);
static const Coordinates ZVector(0,0,1);

class CPUMatrix {
	size_t x,y,z;

public:
	size_t bucket;
	vector<VTYPE> *numbers;
    vector<CPU> CPUs;
	CPUMatrix(char * dimensions, vector<VTYPE> &numbers) {
		this->numbers = &numbers;
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
		for(vector<CPU>::iterator i = this->CPUs.begin(); i!=this->CPUs.end(); ++i,++pos) {
			//handle crazy case where there are cpus with no numbers
			if(pos*bucket > numbers.size()) {
				cout << "Warning, cpus without work..." << endl;
				*i = *(new CPU(0, 0,pos*2*bucket));
			}
			else if((pos+1)*bucket > numbers.size()) {
				*i = *(new CPU(pos*bucket, numbers.size(), pos*2*bucket));
			}
			else {
				*i = *(new CPU(pos*bucket, (pos+1)*bucket, pos*2*bucket));
			}
		}
	}
    CPU* get(size_t x, size_t y, size_t z) {
		if(!this->coorValid(x,y,z)) {
			throw "mimo rozmer procesoru";
		}
        return &CPUs.at(x + y * this->x + z * this->x * this->y);
    }
    CPU* get(Coordinates c) {
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
	CPU* get(size_t i) {
		return (&CPUs.at(i));
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
					this->get(x,y,z)->printData(this->numbers->data());
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

__device__ int getCpuIndex(Coordinates override = Coordinates(0,0,0)) {
	if(threadIdx.x+override.x >= blockDim.x || blockIdx.x+override.y >= gridDim.x || blockIdx.y+override.z >= gridDim.y) {
		return -1;
	}
	if(((int)threadIdx.x+override.x) < 0 || ((int)blockIdx.x+override.y) < 0 || ((int)blockIdx.y+override.z) < 0) {
		return -1;
	}
	return (threadIdx.x + override.x) + (blockIdx.x+override.y) * blockDim.x + (blockIdx.y+override.z) * gridDim.x * blockDim.x;
}
__device__ int getIndexForVector(Coordinates vector) {
	if(vector.x != 0)
		return threadIdx.x;
	if(vector.y != 0)
		return blockIdx.x;
	if(vector.z !=0)
		return blockIdx.y;
	
	return -1;
}

__device__ inline void mergeAndSplit(VTYPE* d_data, VTYPE* d_temp, size_t tempindex,size_t list1min, size_t list1max, size_t list2min, size_t list2max) {
	VTYPE * list3 = d_temp+tempindex;

	int index1 = list1min, index2 = list2min, index3 = 0;

    // Loop untill both arrays have reached their upper bound.
    while (index1 < list1max || index2 < list2max) {

        // Make sure the first array hasn't reached 
        // its upper bound already and make sure we 
        // don't compare outside bounds of the second 
        // array.
        // order of conditions very important, otherwise outside of bounds access!
        if (index2 >= list2max || (index1 < list1max && d_data[index1] <= d_data[index2]) ) {
            list3[index3] = d_data[index1];
            index1++;
        }
        else {
            list3[index3] = d_data[index2];
            index2++;
        }
        index3++;
    }
	index1 = list1min;
	index2 = list2min;
	index3 = 0;
	while(index1<list1max) {
		d_data[index1] = list3[index3];
		index1++;
		index3++;
	}
	while(index2<list2max) {
		d_data[index2] = list3[index3];
		index2++;
		index3++;
	}
}
__global__ void eotSortIterationCuda(CPU * d_cpu, VTYPE* d_data, Coordinates direction, VTYPE*d_temp,bool even, Coordinates snake = Coordinates(0,0,0), bool reverse = false, Coordinates skipvector = Coordinates(0,0,0)) {
	int myindex = getCpuIndex();
	if(snake.x == 1 || snake.y ==1|| snake.z ==1) {
		if(getIndexForVector(snake)%2 ==1) {
			direction = direction * -1;
		}
	}
	if(skipvector.x == 1 || skipvector.y ==1|| skipvector.z ==1) {
		if(getIndexForVector(skipvector)%2 ==0) {
			return;
		}
	}
    if(reverse) {
        direction = direction * -1;
    }
	int oindex = getCpuIndex(direction);
	if(oindex < 0) {
		return;
	}
	if(getIndexForVector(direction)%2 == 0 && even) {
		mergeAndSplit(d_data, d_temp,d_cpu[myindex].tempindex,d_cpu[myindex].startindex, d_cpu[myindex].endindex, d_cpu[oindex].startindex, d_cpu[oindex].endindex);
	}
	if(getIndexForVector(direction)%2 == 1 && !even) {
		mergeAndSplit(d_data, d_temp,d_cpu[myindex].tempindex,d_cpu[myindex].startindex, d_cpu[myindex].endindex, d_cpu[oindex].startindex, d_cpu[oindex].endindex);
	}
}

//knihovni implementace heapsortu
__global__ void localSort(CPU * d_cpu, VTYPE * d_data) {
	int vSize = d_cpu[getCpuIndex()].endindex - d_cpu[getCpuIndex()].startindex;
	VTYPE * vec = d_data+d_cpu[getCpuIndex()].startindex;
	int i;
    int originalVal;
    int promoteIndx;
    int parentIndx;

    for (i = 1; i < vSize; i++) {
        originalVal = vec[i];
        promoteIndx = i;
        parentIndx = (promoteIndx-1) / 2;
        while (promoteIndx > 0 && vec[parentIndx] < originalVal) {
            vec[promoteIndx] = vec[parentIndx];
            promoteIndx = parentIndx;
            parentIndx = (promoteIndx-1) / 2;
        }
        vec[promoteIndx] = originalVal;
    }

    int bottom;
    int displacedVal;
    int vacantNodeIndx;
    int leftIndx;
    int rightIndx;
    int maxIndx;

    for (bottom = vSize-1; bottom > 0; bottom--) {
        displacedVal = vec[bottom];
        vec[bottom] = vec[0];
        // ASSERT: Value in root moved to current bottom of vec

        vacantNodeIndx = 0;
        while (true) {
            leftIndx = 2*vacantNodeIndx + 1;
            if (leftIndx >= bottom)
                break;
            rightIndx = 2*vacantNodeIndx + 2;
            if (rightIndx >= bottom || vec[leftIndx] > vec[rightIndx])
                maxIndx = leftIndx;
            else
                maxIndx = rightIndx;
            if (vec[maxIndx] <= displacedVal)
                break;
            vec[vacantNodeIndx] = vec[maxIndx];
            vacantNodeIndx = maxIndx;
       }
       vec[vacantNodeIndx] = displacedVal;
       // ASSERT: Heap has been recreated in vec[0..bottom-1]
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
void writeFileSnake(ostream &out, Coordinates plane, Coordinates rowvector, Coordinates colvector, CPUMatrix * matrix, bool reverse = false) {
	
	if (out) {
		size_t rows = matrix->sizeForVector(colvector);
		size_t cols = matrix->sizeForVector(rowvector);
		for(size_t row =0;row<rows;row++) {
			if(row%2==0) {
				for(size_t col = 0; col<cols;col++) {
                    if(reverse) {
						matrix->get(plane+rowvector*col+colvector*row)->printDataReverse(matrix->numbers->data(),'\n', out);
                    }
                    else {
						matrix->get(plane+rowvector*col+colvector*row)->printData(matrix->numbers->data(),'\n', out);
                    }
				}
			}
			else {
				for(size_t col = cols-1; col>0;col--) {
                    if(reverse)
						matrix->get(plane+rowvector*col+colvector*row)->printDataReverse(matrix->numbers->data(),'\n', out);
					else
						matrix->get(plane+rowvector*col+colvector*row)->printData(matrix->numbers->data(),'\n', out);
				}
                if(reverse)
					matrix->get(plane+colvector*row)->printDataReverse(matrix->numbers->data(),'\n', out);
                else
					matrix->get(plane+colvector*row)->printData(matrix->numbers->data(),'\n', out);
			}
		}

    }
    else {
        throw "chyba vypisu dat";
    }
}
void eotSortCuda(Coordinates dynvector, CPUMatrix &matrix, CPU * d_cpu, VTYPE * d_data, VTYPE * d_temp, Coordinates snake = Coordinates(0,0,0), bool reverse = false, Coordinates skipvector = Coordinates(0,0,0), size_t override_runs = 0) {

	size_t dimension_size = matrix.sizeForVector(dynvector);
    size_t runtimes = dimension_size/2+1;
	if(override_runs > 0) {
		runtimes = override_runs;
	}
	//TODO: prozkoumat to +1
	for(size_t count=0;count<runtimes;count++) {
		//even pairs
		eotSortIterationCuda<<<dim3(matrix.sizeForVector(YVector), matrix.sizeForVector(ZVector),1), matrix.sizeForVector(XVector)>>> (d_cpu,d_data, dynvector, d_temp,true, snake, reverse, skipvector);
		//odd pairs
		eotSortIterationCuda<<<dim3(matrix.sizeForVector(YVector), matrix.sizeForVector(ZVector),1), matrix.sizeForVector(XVector)>>> (d_cpu,d_data, dynvector, d_temp,false, snake, reverse, skipvector);
	}
}

void shearSortCuda(vector<VTYPE> &data, CPUMatrix &matrix,CPU*d_cpu,VTYPE* d_data,VTYPE* d_temp, Coordinates xvector = XVector, Coordinates yvector = YVector, bool reverse = false, Coordinates skipvector = Coordinates(0,0,0)) {
	size_t rowsize = matrix.sizeForVector(xvector);
	size_t colsize = matrix.sizeForVector(yvector);
	for(size_t count=0;count<logbin(colsize)+2;count++) {
		eotSortCuda(yvector,matrix,d_cpu,d_data, d_temp, Coordinates(0,0,0), reverse, skipvector);
		eotSortCuda(xvector,matrix,d_cpu,d_data, d_temp, yvector, reverse, skipvector);
	}

}
void dSortCuda(vector<VTYPE> &data, CPUMatrix &matrix,CPU*d_cpu,VTYPE* d_data,VTYPE* d_temp) {

	size_t xsize = matrix.sizeForVector(XVector);
	size_t ysize = matrix.sizeForVector(YVector);
	size_t zsize = matrix.sizeForVector(ZVector);
    Coordinates zero(0,0,0);
	shearSortCuda(data, matrix, d_cpu, d_data, d_temp, XVector, YVector);
	shearSortCuda(data, matrix, d_cpu, d_data, d_temp, ZVector, YVector);
	shearSortCuda(data, matrix, d_cpu, d_data, d_temp, XVector, ZVector, true);
	shearSortCuda(data, matrix, d_cpu, d_data, d_temp, XVector, ZVector, false, YVector);
	eotSortCuda(YVector,matrix,d_cpu,d_data, d_temp, Coordinates(0,0,0), false, Coordinates(0,0,0), 1);
	shearSortCuda(data, matrix, d_cpu, d_data, d_temp, XVector, ZVector, false);
}

int main(int argc, char** argv) {
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
            struct timeval time;
            gettimeofday(&time, NULL);
            double t1=time.tv_sec+(time.tv_usec/1000000.0);
			//CUDA INIT
			CPU * d_cpu;
			const size_t sz = size_t(matrix.CPUs.size()) * sizeof(CPU);
			gpuErrchk(cudaMalloc((void**)&d_cpu, sz));
			cudaMemcpy(d_cpu, matrix.CPUs.data(), sz, cudaMemcpyHostToDevice);
			//copy data to gpu
			VTYPE * d_data;
			const size_t sd = size_t(data.size()) * sizeof(VTYPE);
			cudaMalloc((void**)&d_data, sd);
			cudaMemcpy(d_data, data.data(), sd, cudaMemcpyHostToDevice);
			//docasny prostor pro merge
			VTYPE * d_temp;
			const size_t st = 2*matrix.bucket * sizeof(VTYPE)*matrix.CPUs.size();
			cudaMalloc((void**)&d_temp, st);
			localSort<<<dim3(matrix.sizeForVector(YVector), matrix.sizeForVector(ZVector),1), matrix.sizeForVector(XVector)>>> (d_cpu,d_data);
			shearSortCuda(data, matrix,d_cpu,d_data,d_temp, XVector, YVector);

			cudaMemcpy(data.data(), d_data, sd, cudaMemcpyDeviceToHost);
            gettimeofday(&time, NULL);
            double t2=time.tv_sec+(time.tv_usec/1000000.0);
            printf("Sorting took %.6lf seconds\n", t2-t1);
			writeFileSnake(outputStream, Coordinates(0,0,0), XVector,YVector, &matrix);

			return 0;
		}
        
		
		dimensions = CLIArgumentsParser::getCmdOption(argv, argv + argc, "-3");
		if(dimensions) {
			//3DSort
			vector<VTYPE> data;
			readFile(filename, data);

			CPUMatrix matrix(dimensions, data);
            if(matrix.sizeForVector(ZVector) < 2) {
                cout << "3d sort potrebuje z souradnici" << endl;
                return -1;
            }
            struct timeval time;
            gettimeofday(&time, NULL);
            double t1=time.tv_sec+(time.tv_usec/1000000.0);
			//CUDA INIT
			CPU * d_cpu;
			const size_t sz = size_t(matrix.CPUs.size()) * sizeof(CPU);
			gpuErrchk(cudaMalloc((void**)&d_cpu, sz));
			cudaMemcpy(d_cpu, matrix.CPUs.data(), sz, cudaMemcpyHostToDevice);
			//copy data to gpu
			VTYPE * d_data;
			const size_t sd = size_t(data.size()) * sizeof(VTYPE);
			cudaMalloc((void**)&d_data, sd);
			cudaMemcpy(d_data, data.data(), sd, cudaMemcpyHostToDevice);
			//docasny prostor pro merge
			VTYPE * d_temp;
			const size_t st = 2*matrix.bucket * sizeof(VTYPE)*matrix.CPUs.size();
			cudaMalloc((void**)&d_temp, st);
			localSort<<<dim3(matrix.sizeForVector(YVector), matrix.sizeForVector(ZVector),1), matrix.sizeForVector(XVector)>>> (d_cpu,d_data);
			dSortCuda(data, matrix,d_cpu,d_data,d_temp);

			cudaMemcpy(data.data(), d_data, sd, cudaMemcpyDeviceToHost);
            gettimeofday(&time, NULL);
            double t2=time.tv_sec+(time.tv_usec/1000000.0);
            printf("Sorting took %.6lf seconds\n", t2-t1);
            for(int y=0;y<matrix.sizeForVector(YVector);y++) {
				writeFileSnake(outputStream, Coordinates(0,y,0), XVector,ZVector, &matrix);
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