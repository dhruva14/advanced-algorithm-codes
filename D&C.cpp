#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

float order;

//sum two matrix and store the result
void sum(vector< vector<float> > &A, 
         vector< vector<float> > &B, 
         vector< vector<float> > &C, float size) {
   
	for (float i = 0; i < size; i++) {
	        for (float j = 0; j < size; j++) {
	            C[i][j] = A[i][j] + B[i][j];
	        }
	}
}
//divide and conquer takes the matrices and their starting,ending row and column indexs
void d_c(vector <vector <float> > &A,
	vector <vector <float> > &B,
	vector <vector <float> > &C,float size,
	float row1,float column1,float row11,float column11,
	float row2,float column2,float row22,float column22){
	
	if(size<=order/8){//base size is (order of matrix/8) for better and faster execution time(optimisation based on system)
		for(float i=0;i<size;i++){ //initialising result matrix to zero
			for(float j=0;j<size;j++){
				C[i][j]=0;	
			}
		}
	for(float i=row1,rindex=0;(i<=row11);i++,rindex++){//iterative matrix multiplication
		for(float k=column2,cindex=0;k<=column22;k++,cindex++){
			for(float j=row2,j1=column1;j<=row22;j++,j1++){
				C[rindex][cindex]+=A[i][j1]*B[j][k];
			}
			
		}
	}	

}
	else{
		float size1=size/2;//divide the size
		vector <float> row(size1);
		vector < vector <float> > c11(size1,row),c12(size1,row),c21(size1,row),c22(size1,row),res1(size1,row),res2(size1,row);
	
		d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);//A11*B11
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);//A12*B21
		
		sum(res1,res2,c11,size1);
	
		d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);//A11*B21
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);//A12*B22
	
		sum(res1,res2,c12,size1);
	
		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);//A21*B11
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);//A22*B21
	
		sum(res1,res2,c21,size1);

		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);//A21*B12
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);//A22*B22

		sum(res1,res2,c22,size1);
	
		for(float i=0;i<size1;i++){
			for(float j=0;j<size1;j++){
				C[i][j]=c11[i][j];
				C[i][j+size1]=c12[i][j];
				C[i+size1][j]=c21[i][j];
				C[i+size1][j+size1]=c22[i][j];
			
			}
		}
	

	}
}

int main(){
	float size;
	cin>>size;
	order=size;
	clock_t t;
	vector <float> row(size);
	vector < vector <float> > A(size,row),B(size,row),C(size,row);
	for(float i=0;i<size;i++){
		for(float j=0;j<size;j++){	
			cin>>A[i][j];
		}
	}
	for(float i=0;i<size;i++){
		for(float j=0;j<size;j++){	
			cin>>B[i][j];
		}
		
	}
	t=clock();
	d_c(A,B,C,size,0,0,size-1,size-1,0,0,size-1,size-1);
	t=clock()-t;
	
	
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	
	
	for(float i=0;i<size;i++){
		for(float j=0;j<size;j++){	
			cout<<C[i][j]<<" ";
		}
		cout<<endl;
	}
	printf("d&c time taken is %f seconds\n",time_taken);
	return 0;
}
