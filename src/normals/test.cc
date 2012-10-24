/*******************************************************
 A simple program that demonstrates NewMat10 library.
 The program defines a random symmetric matrix
 and computes its eigendecomposition.
 For further details read the NewMat10 Reference Manual
********************************************************/

#include <stdlib.h>
#include <time.h>
#include <string.h>

// the following two are needed for printing
#include <iostream.h>
#include <iomanip.h>
/**************************************
/* The NewMat10 include files         */
#include <include.h>
#include <newmat.h>
#include <newmatap.h>
#include <newmatio.h>
/***************************************/


int main(int argc, char **argv) {
	int M = 3, N = 5;
	Matrix X(M,N); // Define an M x N general matrix

	// Fill X by random numbers between 0 and 9
	// Note that indexing into matrices in NewMat is 1-based!
	srand(time(NULL));
	for (int i = 1; i <= M; ++i) {
		for (int j = 1; j <= N; ++j) {
			X(i,j) = rand() % 10;
		}
	}

	SymmetricMatrix C;
	C << X * X.t(); // fill in C by X * X^t. 
					// Works because we *know* that the result is symmetric

	cout << "The symmetrix matrix C" << endl;
	cout << setw(5) << setprecision(0) << C << endl;
	

	// compute eigendecomposition of C
	Matrix			V(3,3); // for eigenvectors
	DiagonalMatrix	D(3);   // for eigenvalues

	// the decomposition
	Jacobi(C, D, V);
	
	// Print the result
	cout << "The eigenvalues matrix:" << endl;
	cout << setw(10) << setprecision(5) << D << endl;
	cout << "The eigenvectors matrix:" << endl;
	cout << setw(10) << setprecision(5) << V << endl;

	// Check that the first eigenvector indeed has the eigenvector property
	ColumnVector v1(3);
	v1(1) = V(1,1);
	v1(2) = V(2,1);
	v1(3) = V(3,1);

	ColumnVector Cv1 = C * v1;
	ColumnVector lambda1_v1 = D(1) * v1;

	cout << "The max-norm of the difference between C*v1 and lambda1*v1 is " <<
		   NormInfinity(Cv1 - lambda1_v1) << endl << endl;

	// Build the inverse and check the result
	Matrix Ci = C.i();
	Matrix I  = Ci * C;

	cout << "The inverse of C is" << endl;
	cout << setw(10) << setprecision(5) << Ci << endl;
	cout << "And the inverse times C is identity" << endl;
	cout << setw(10) << setprecision(5) << I << endl;

	// Example for multiple solves (see NewMat documentation)
	ColumnVector r1(3), r2(3);
	for (i = 1; i <= 3; ++i) {
		r1(i) = rand() % 10;
		r2(i) = rand() % 10;
	}
	LinearEquationSolver CLU = C; // decomposes C
	ColumnVector s1 = CLU.i() * r1;
	ColumnVector s2 = CLU.i() * r2;

	cout << "solution for right hand side r1" << endl;
	cout << setw(10) << setprecision(5) << s1 << endl;
	cout << "solution for right hand side r2" << endl;
	cout << setw(10) << setprecision(5) << s2 << endl;
	
	return 0;
}
