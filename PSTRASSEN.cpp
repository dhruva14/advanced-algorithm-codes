#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>

using namespace std;

float order;

void lsum(vector< vector<float> > &A, 
         vector< vector<float> > &B, 
         vector< vector<float> > &C, float size,float row1,float column1,float row2,float column2) {
	
	    for (float i = row1,i1=row2,rindex=0; rindex < size; i++,i1++,rindex++) {
	        for (float j = column1,j1=column2,cindex=0; cindex < size; j++,j1++,cindex++) {
	            C[rindex][cindex] = A[i][j] + B[i1][j1];
	        }
	    }
}
void sum(vector< vector<float> > &A, 
              vector< vector<float> > &B, 
              vector< vector<float> > &C, float size) {
    
	for (float i = 0;i < size; i++) {
	        for (float j = 0; j < size; j++) {
	            C[i][j] = A[i][j] + B[i][j];
        	}
    	}   
}
void lsubtract(vector< vector<float> > &A, 
         vector< vector<float> > &B, 
         vector< vector<float> > &C, float size,float row1,float column1,float row2,float column2) {
    for (float i = row1,i1=row2,rindex=0; rindex < size; i++,i1++,rindex++) {
        for (float j = column1,j1=column2,cindex=0; cindex< size; j++,j1++,cindex++) {
            C[rindex][cindex] = (A[i][j] - B[i1][j1]);
        }
    }
}
void subtract(vector< vector<float> > &A, 
              vector< vector<float> > &B, 
              vector< vector<float> > &C, float size) {
	
    float i, j;
	for (i = 0; i < size; i++) {
	        for (j = 0; j < size; j++) {
	            C[i][j] = (A[i][j] - B[i][j]);
        	}
    	}   
}
void strassen(vector< vector<float> > &A, 
              vector< vector<float> > &B, 
              vector< vector<float> > &C, float size,
	      float row1,float column1,float row11,float column11,
	      float row2,float column2,float row22,float column22) {
	if (size <= order/8) {
        	for(float i=0;i<size;i++){
			for(float j=0;j<size;j++){
				C[i][j]=0;	
			}
		}
		for(float i=row1,rindex=0;(i<row11);i++,rindex++){
			for(float k=column2,cindex=0;k<column22;k++,cindex++){
				for(float j=row2,j1=column1;j<row22;j++,j1++){
					C[rindex][cindex]+=A[i][j1]*B[j][k];
				}
				
			}
		}
		
        return;
}

    
	else{
        	float size1 = size/2;
        	vector<float> row (size1);
        	vector< vector<float> > c11(size1,row), c12(size1,row), c21(size1,row), c22	(size1,row),
            	p1(size1,row), p2(size1,row), p3(size1,row), p4(size1,row), 
            	p5(size1,row), p6(size1,row), p7(size1,row),
            	aResult(size1,row), bResult(size1,row);

       		if(size>32){//to make sure execution time using threads is not more than that without threads
	

	        	lsum(A,A, aResult, size1,0+row1,0+column1,size1+row1,size1+column1); // a11 + a22
			lsum(B,B, bResult, size1,0+row2,0+column2,size1+row2,size1+column2); // b11 + b22
		
			thread thread1(strassen,ref(aResult), ref(bResult), ref(p1), size1,0,0,size1,size1,0,0,size1,size1);
			thread1.join();

	        	lsum(A,A,aResult, size1,size1+row1,0+column1,size1+row1,size1+column1); // a21 + a22
	
			thread thread2(strassen,ref(aResult),ref(B), ref(p2), size1,0,0,size1,size1,0+row2,0+column2,size1+row2,size1+column2);
			thread2.join();

	        	lsubtract(B,B, bResult, size1,0+row2,size1+column2,size1+row2,size1+column2); // b12 - b22

			thread thread3(strassen,ref(A), ref(bResult), ref(p3), size1,0+row1,0+column1,size1+row1,size1+column1,0,0,size1,size1);
			thread3.join();
        
			lsubtract(B,B, bResult, size1,size1+row2,0+column2,0+row2,0+column2); // b21 - b11
	
			thread thread4(strassen,ref(A), ref(bResult), ref(p4), size1,size1+row1,size1+column1,size+row1,size+column1,0,0,size1,size1);
			thread4.join();
		
			lsum(A,A, aResult,size1,0+row1,0+column1,0+row1,size1+column1); // a11 + a12
	
			thread thread5(strassen,ref(aResult),ref(B), ref(p5), size1,0,0,size1,size1,size1+row2,size1+column2,size+row2,size+column2);
           		thread5.join();
        
			lsubtract(A,A, aResult,size1,row1+0,column1+0,row1+size1,0+column1); // a11 - a21
	        	lsum(B,B, bResult,size1,0+row2,0+column2,0+row2,size1+column2); // b11 + b12
	
			thread thread6(strassen,ref(aResult), ref(bResult), ref(p6), size1,0,0,size1,size1,0,0,size1,size1);
			thread6.join();

	        	lsubtract(A,A, aResult, size1,0+row1,size1+column1,size1+row1,size1+column1); // a12 - a22
	        	lsum(B,B, bResult, size1,size1+row2,0+column2,size1+row2,size1+column2); // b21 + b22

			thread thread7(strassen,ref(aResult), ref(bResult), ref(p7), size1,0,0,size1,size1,0,0,size1,size1);
			thread7.join();
        	}
	else{
			lsum(A,A, aResult, size1,0+row1,0+column1,size1+row1,size1+column1); // a11 + a22
			lsum(B,B, bResult, size1,0+row2,0+column2,size1+row2,size1+column2); // b11 + b22
			strassen(aResult, bResult, p1, size1,0,0,size1,size1,0,0,size1,size1); // p1 = (a11+a22) * (b11+b22)
	
	        	lsum(A,A,aResult, size1,size1+row1,0+column1,size1+row1,size1+column1); // a21 + a22
			strassen(aResult,B, p2, size1,0,0,size1,size1,0+row2,0+column2,size1+row2,size1+column2); // p2 = (a21+a22) * (b11)
	
	        	lsubtract(B,B, bResult, size1,0+row2,size1+column2,size1+row2,size1+column2); // b12 - b22
			strassen(A, bResult, p3, size1,0+row1,0+column1,size1+row1,size1+column1,0,0,size1,size1); // p3 = (a11) * (b12 - b22)

		        lsubtract(B,B, bResult, size1,size1+row2,0+column2,0+row2,0+column2); // b21 - b11
			strassen(A, bResult, p4, size1,size1+row1,size1+column1,size+row1,size+column1,0,0,size1,size1); // p4 = (a22) * (b21 - b11)
		        lsum(A,A, aResult,size1,0+row1,0+column1,0+row1,size1+column1); // a11 + a12
		        strassen(aResult,B, p5, size1,0,0,size1,size1,size1+row2,size1+column2,size+row2,size+column2); // p5 = (a11+a12) * (b22)   
		        lsubtract(A,A, aResult,size1,row1+0,column1+0,row1+size1,0+column1); // a11 - a21
			lsum(B,B, bResult,size1,0+row2,0+column2,0+row2,size1+column2); // b11 + b12
			strassen(aResult, bResult, p6, size1,0,0,size1,size1,0,0,size1,size1); // p6 = (a11-a21) * (b11+b12)
	
        		lsubtract(A,A, aResult, size1,0+row1,size1+column1,size1+row1,size1+column1); // a12 - a22
        		lsum(B,B, bResult, size1,size1+row2,0+column2,size1+row2,size1+column2); // b21 + b22
       		        strassen(aResult, bResult, p7, size1,0,0,size1,size1,0,0,size1,size1); // p7 = (a12-a22) * (b21+b22)

	}

        // c21, c21, c11 e c22:

        sum(p3, p5, c12, size1); // c12 = p3 + p5
        sum(p2, p4, c21, size1); // c21 = p2 + p4

        sum(p1, p4, aResult, size1); // p1 + p4
        sum(aResult, p7, bResult, size1); // p1 + p4 + p7
        subtract(bResult, p5, c11, size1); // c11 = p1 + p4 - p5 + p7

        sum(p1, p3, aResult, size1); // p1 + p3
        subtract(aResult, p6, bResult, size1); // p1 + p3 - p6
        subtract(bResult, p2, c22, size1); // c22 = p1 + p3 - p2 - p6

        
        for (float i = 0; i < size1 ; i++) {
            for (float j = 0 ; j < size1 ; j++) {
                C[i][j] = c11[i][j];
                C[i][j + size1] = c12[i][j];
                C[i + size1][j] = c21[i][j];
                C[i + size1][j + size1] = c22[i][j];
            }
        }
    }
}
int main () {
	clock_t t;
	float n;
	cin>>n;
	order=n;
    	vector<float> row (n);
    	vector< vector<float> > A(n, row), B(n, row), C(n, row);
	for(float i=0;i<n;i++){
		for(float j=0;j<n;j++){
			cin>>A[i][j];
		}
	}
	for(float i=0;i<n;i++){
		for(float j=0;j<n;j++){
			cin>>B[i][j];
		}
	}
	t = clock();
	strassen(A,B,C,n,0,0,n,n,0,0,n,n);
	t = clock() - t;
	
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	
	for(float i=0;i<n;i++){
		for(float j=0;j<n;j++){	
			cout<<C[i][j]<<" ";
		}
		cout<<endl;
	}
	
	printf("PStrassen took %f seconds to execute \n", time_taken);
    return 0;
}
