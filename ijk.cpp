#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

void ijk1(vector< vector<float> > &A, 
                                   vector< vector<float> > &B,
                                   vector< vector<float> > &C, float n) {
    for (float i = 0; i < n; i++) {
        for (float j = 0; j < n; j++) {
            for (float k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

int main(){
	float size;
	cin>>size;
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
	ijk1(A,B,C,size);
	t=clock()-t;
	
	
	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
	printf("ijk time taken is %f seconds\n",time_taken);
	for(float i=0;i<size;i++){
		for(float j=0;j<size;j++){	
			cout<<C[i][j]<<" ";
		}
		cout<<endl;
	}
	return 0;
}
