#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>


using namespace std;

int order;
//Sum two matrix and store
void sum(vector< vector<int> > &A, 
         vector< vector<int> > &B, 
         vector< vector<int> > &C, int size) {
    
	for (int i = 0; i < size; i++) {
	        for (int j = 0; j < size; j++) {
	            C[i][j] = A[i][j] + B[i][j];
	        }
	}
}
//divide and conquer with starting and ending index
void d_c(vector <vector <int> > &A,
	vector <vector <int> > &B,
	vector <vector <int> > &C,int size,
	int row1,int column1,int row11,int column11,
	int row2,int column2,int row22,int column22){

	if(size<=order/8){
		for(int i=0;i<size;i++){
			for(int j=0;j<size;j++){
				C[i][j]=0;	
			}
		}
		for(int i=row1,rindex=0;(i<=row11);i++,rindex++){
			for(int k=column2,cindex=0;k<=column22;k++,cindex++){
				for(int j=row2,j1=column1;j<=row22;j++,j1++){
					C[rindex][cindex]+=A[i][j1]*B[j][k];
				}		
			}
		}	

	}
	else{	
		
		int size1=size/2;
		vector <int> row(size1);
		vector < vector <int> > c11(size1,row),c12(size1,row),c21(size1,row),c22(size1,row),res1(size1,row),res2(size1,row);

		if((size>=order/2)&&(size>32)){
		d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);
		thread thread1(sum,ref(res1),ref(res2),ref(c11),size1);
		//sum(res1,res2,c11,size1);
		thread1.join();
		d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);
		thread thread2(sum,ref(res1),ref(res2),ref(c12),size1);
		//sum(res1,res2,c12,size1);
		thread2.join();
		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);
		thread thread3(sum,ref(res1),ref(res2),ref(c21),size1);	
		//sum(res1,res2,c21,size1);
		thread3.join();
		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);
		thread thread4(sum,ref(res1),ref(res2),ref(c22),size1);
		//sum(res1,res2,c22,size1);
		thread4.join();
		}
		else{
			d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);
		//thread thread1(sum,ref(res1),ref(res2),ref(c11),size1);
		sum(res1,res2,c11,size1);
		//thread1.join();
		d_c(A,B,res1,size1,row1,column1,size1-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);
		d_c(A,B,res2,size1,row1,size1+column1,size1-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);
		//thread thread2(sum,ref(res1),ref(res2),ref(c12),size1);
		sum(res1,res2,c12,size1);
		//thread2.join();
		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,column2,size1-1+row2,size1-1+column2);
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,column2,size-1+row2,size1-1+column2);
		//thread thread3(sum,ref(res1),ref(res2),ref(c21),size1);	
		sum(res1,res2,c21,size1);
		//thread3.join();
		d_c(A,B,res1,size1,size1+row1,column1,size-1+row1,size1-1+column1,row2,size1+column2,size1-1+row2,size-1+column2);
		d_c(A,B,res2,size1,size1+row1,size1+column1,size-1+row1,size-1+column1,size1+row2,size1+column2,size-1+row2,size-1+column2);
		//thread thread4(sum,ref(res1),ref(res2),ref(c22),size1);
		sum(res1,res2,c22,size1);
		//thread4.join();


		}
		for(int i=0;i<size1;i++){
			for(int j=0;j<size1;j++){
				C[i][j]=c11[i][j];
				C[i][j+size1]=c12[i][j];
				C[i+size1][j]=c21[i][j];
				C[i+size1][j+size1]=c22[i][j];
			
			}
		}
	
	
	
	

	}

}

int main(){
	int size;
	cin>>size;
	order=size;
	clock_t t;
	vector <int> row(size);
	vector < vector <int> > A(size,row),B(size,row),C(size,row);
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){	
			cin>>A[i][j];
		}
	}
	for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){	
			cin>>B[i][j];
		}
		
	}
	t=clock();
	d_c(A,B,C,size,0,0,size-1,size-1,0,0,size-1,size-1);
	t=clock()-t;
	
	
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	
	printf("pdc time taken is %f seconds\n",time_taken);
	/*for(int i=0;i<size;i++){
		for(int j=0;j<size;j++){	
			cout<<C[i][j]<<" ";
		}
		cout<<endl;
	}*/
	return 0;
}
