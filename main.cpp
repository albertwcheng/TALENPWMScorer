#include "CppUtilClasses/Matrix.h"
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class PosMatrix
{
	public:
	Matrix<double> matData;
	void readFile(const string&filename){
			
	}
};




int main(int argc,char** argv){
	Matrix<double> M;
	vector<double> row;
	row.push_back(0.0);
	row.push_back(0.5);
	row.push_back(0.7);
	row.push_back(0.8);
	
	for(int i=0;i<3;i++){
		for(int j=0;j<4;j++){
			row[j]+=1;	
		}	
		
		M.addRow(row);
	}
	
	M.print(cout);
	cerr<<"change 0,0"<<endl;
	M.itemAt(0,0)=-20;
	M.print(cout);
	M.finalize();
	cout<<"finalize"<<endl;
	M.print(cout);
	M.itemAt(1,2)=-10;
	cout<<"change 1,2"<<endl;
	M.print(cout);
	cout<<"read matrix"<<endl;
	Matrix<double> Mf;
	StringToFloat stf;
	Mf.readFile("matrix.dat",stf);
	Mf.print(cout);
	Mf.finalize();
	cout<<"finalize read matrix"<<endl;
	Mf.print(cout);
}