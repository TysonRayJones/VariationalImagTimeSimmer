/** @file 
 * Builds diagonal Hamiltonians from 3SAT problems
 */

#include <QuEST.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include "sat_generator.h"
#include "hamiltonian_builder.h"
 

/**
 * Builds a diagonal hamiltonian from the given 3SAT equation,
 * where the energy of a state (bitstring) is the total number of
 * clauses that it fails. The returned hamiltonian must be freed!
 */
double* getDiagHamilFrom3SAT(int* equ, int numBools, int numClauses) {
	
	long long int length = pow(2LL, numBools);
	double* hamiltonian = malloc(length * sizeof *hamiltonian);
	
	// for every basis vector
	for (long long int basis=0LL; basis < length; basis++) {
		
		double energy = 0;
		
		// check the energy contribution of each clause
		for (int term=0; term < 3*numClauses; term+=3) {
			
			// bool indices in the clause
			int inds[3];
			for (int i=0; i<3; i++)
				inds[i] = abs(equ[term+i])-1;
			
			// bad bool values
			int vals[3];
			for (int i=0; i<3; i++)
				vals[i] = equ[term+i] < 0;
			
			// bits of basis at indices (indexed from LEFT)
			int bits[3];
			for (int i=0; i<3; i++)
				bits[i] = (basis >> (numBools-inds[i]-1)) & 1;
				
			// increase basis energy if basis fails clause
			if (bits[0] == vals[0] && 
				bits[1] == vals[1] &&
				bits[2] == vals[2])
					energy += 1;
		}
		
		hamiltonian[basis] = energy;
	}
	
	// must be freed!
	return hamiltonian;
}




/**
 * Given a filename, loads a 3SAT equation (setting equ, numBools and numClauses),
 * finds its solution (setting sol) by brute-force, and builds the corresponding
 * hamiltonian (setting hamil) where the energy of a state (bitstring) is the number
 * of clauses in equ that it fails. equ, sol and hamil must be freed!
 */
Hamiltonian load3SATAndHamilFromFile(
	char *filename, 
	int **equ, int **sol,
	int *numBools, int *numClauses
) {
	
	// load the equation, infer numBools and numClauses
	*equ = loadEquation(filename, numBools, numClauses);
	
	// generate the solution
	*sol = malloc(*numBools * sizeof(int));
	int *candidate = malloc(*numBools * sizeof *candidate);
	findSingleSolution(*equ, *numBools, *numClauses, *sol, candidate, pow(2, *numBools));
	free(candidate);
	
	// build the hamiltonian
	Hamiltonian hamil;
	hamil.type = DIAGONAL;
	hamil.numQubits = *numBools;
	hamil.numAmps = pow(2LL, *numBools);
	hamil.diagHamil = getDiagHamilFrom3SAT(*equ, *numBools, *numClauses);
	
	hamil.terms = NULL;
	hamil.termCoeffs = NULL;
	hamil.numTerms = -1;
	
	return hamil;
	
	// equ, sol, hamiltonian must be freed!
}


/**
 * Creates a random 3SAT equation (setting equ and numClauses) of the given size
 * with a unique solution (and other constraints in Simon's paper) (setting sol), 
 * and builds the  corresponding hamiltonian (setting hamil) where the energy of a 
 * state (bitstring) is the number of clauses in equ that it fails. 
 * equ, sol and hamil must be freed!
 */
Hamiltonian getRandom3SATAndHamil(
	int numBools, 
	int **equ, int **sol, int *numClauses
) {
	
	getRandomEquAndSol(numBools, equ, sol, numClauses);
	
	Hamiltonian hamil;
	hamil.type = DIAGONAL;
	hamil.numQubits = numBools;
	hamil.numAmps = pow(2LL, numBools);
	hamil.diagHamil = getDiagHamilFrom3SAT(*equ, numBools, *numClauses);
	
	hamil.terms = NULL;
	hamil.termCoeffs = NULL;
	hamil.numTerms = -1;
	
	return hamil;
	
	// equ, sol hamiltonian must be freed!
}



void print3SATEquSol(int *equ, int *sol, int numBools, int numClauses) {
	
	printf("equ:\n");
	printEquation(equ, numClauses);
	
	printf("\nsol:\n");
	for (int i=0; i<numBools; i++)
		printf("%d", sol[i]);
	printf("\n");
}


/**
 * Prints a diagonal hamiltonian, with format bitstring: energy (as int)
 */
void printDiagHamil(double* hamil, long long int length) {
		
	// for every basis vector
	for (long long int i=0LL; i < length; i++) {
		
		// print its bitstring
		for (unsigned int mask = length >> 1; mask; mask >>= 1)
			printf("%d", !!(mask & i));
		
		// print energy
		printf(": %d\n", (int) hamil[i]);
	}
}



int getPauliHamilFromFile(char *filename, double** coeffs, int*** terms, int *numTerms, int *numQubits) {
	
	/*
	 * file format: coeff {term} \n where {term} is #numQubits values of 
	 * 0 1 2 3 signifying I X Y Z acting on that qubit index
	 */
	FILE* file = fopen(filename, "r");
	
	if (file == NULL) {
		printf("ERROR: hamiltonian file (%s) not found!\n", filename);
		return 1;
	}
	
	// count the number of qubits
	*numQubits = -1;
	char ch;
	while ((ch=getc(file)) != '\n')
		if (ch == ' ')
			*numQubits += 1;
			
	// count the number of terms
	rewind(file);
	*numTerms = 0;
	while ((ch=getc(file)) != EOF)
		if (ch == '\n')
			*numTerms += 1;
	
	// collect coefficients and terms
	rewind(file);
	*coeffs = malloc(*numTerms * sizeof **coeffs);
	*terms = malloc(*numTerms * sizeof **terms);
	for (int t=0; t < *numTerms; t++) {
		
		// record coefficient
		if (fscanf(file, "%lf ", &((*coeffs)[t])) != 1) {
			printf("ERROR: hamiltonian file (%s) has a line with no coefficient!\n", filename);
			return 1;
		}
		
		// record #numQubits operations in the term
		(*terms)[t] = malloc(*numQubits * sizeof ***terms);
		for (int q=0; q < *numQubits; q++)
			if (fscanf(file, "%d ", &((*terms)[t][q])) != 1) {
				printf("ERROR: hamiltonian file (%s) has a line missing some terms!\n", filename);
				return 1;
			}
			
		// the newline is magically eaten
	}
	
	// indicate success
	return 0;
}


Hamiltonian loadPauliHamilFromFile(char *filename) {

	Hamiltonian hamil;
	hamil.type = PAULI_TERMS;
	getPauliHamilFromFile(filename, &hamil.termCoeffs, &hamil.terms, &hamil.numTerms, &hamil.numQubits);
	hamil.numAmps = pow(2LL, hamil.numQubits);
	
	return hamil;
}



/**
 * Prints a hamiltonian of pauli terms with format (coeff) * X2 Y3 X4
 */
void printPauliHamil(double* coeffs, int** terms, int numTerms, int numQubits) {
	
	for (int t=0; t < numTerms; t++) {
		
		// print term coefficient
		printf("(%lf) * ", coeffs[t]);
		
		// print gates and qubit inds
		for (int q=0; q < numQubits; q++) {
			if (terms[t][q] == 1)
				printf("X");
			if (terms[t][q] == 2)
				printf("Y");
			if (terms[t][q] == 3)
				printf("Z");
			if (terms[t][q] != 0)
				printf("%d ", q);
		}
		
		// print a + except after the last term
		if (t < numTerms - 1)
			printf("+");
		printf("\n");
	}
}



void printHamil(Hamiltonian hamil) {
	
	if (hamil.type == DIAGONAL)
		printDiagHamil(hamil.diagHamil, hamil.numAmps);
		
	if (hamil.type == PAULI_TERMS)
		printPauliHamil(hamil.termCoeffs, hamil.terms, hamil.numTerms, hamil.numQubits);
}



void freeHamil(Hamiltonian hamil) {
	
	if (hamil.type == DIAGONAL) {
		free(hamil.diagHamil);
	} 
	
	if (hamil.type == PAULI_TERMS) {
		free(hamil.termCoeffs);
		for (int t=0; t < hamil.numTerms; t++)
			free(hamil.terms[t]);
		free(hamil.terms);
	}
}



void applyDiagHamil(double complex* hamilState, MultiQubit qubits, double* diagHamil) {
	for (long long int i=0LL; i < qubits.numAmps; i++)
		hamilState[i] = diagHamil[i]*(getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
}

void applyPauliHamil(double complex* hamilState, MultiQubit qubits, double* coeffs, int** terms, int numTerms) {
	
	// clear hamilState
	for (long long int i=0LL; i < qubits.numAmps; i++)
		hamilState[i] = 0;
		
	// for every term in the hamiltonian
	for (int t=0; t < numTerms; t++) {
		
		// apply each gate in the term
		for (int q=0; q < qubits.numQubits; q++) {
			if (terms[t][q] == 1)
				sigmaX(qubits, q);
			if (terms[t][q] == 2)
				sigmaY(qubits, q);
			if (terms[t][q] == 3)
				sigmaZ(qubits, q);
		}
		
		// add this term's contribution to the hamilState
		for (long long int i=0LL; i < qubits.numAmps; i++)
			hamilState[i] += coeffs[t] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
			
		// undo our change to qubits, exploiting XX = YY = ZZ = I
		for (int q=0; q < qubits.numQubits; q++) {
			if (terms[t][q] == 1)
				sigmaX(qubits, q);
			if (terms[t][q] == 2)
				sigmaY(qubits, q);
			if (terms[t][q] == 3)
				sigmaZ(qubits, q);
		}
	}
}

void applyHamil(double complex* hamilState, MultiQubit qubits, Hamiltonian hamil) {
	
	if (hamil.type == DIAGONAL)
		applyDiagHamil(hamilState, qubits, hamil.diagHamil);
	
	if (hamil.type == PAULI_TERMS)
		applyPauliHamil(hamilState, qubits, hamil.termCoeffs, hamil.terms, hamil.numTerms);
}


double getExpectedEnergy(double complex* hamilState, MultiQubit qubits, Hamiltonian hamil) {
	
	if (hamil.type == DIAGONAL) {
		double energy = 0;
		for (long long int i=0LL; i < hamil.numAmps; i++)
			energy += getProbEl(qubits, i)*hamil.diagHamil[i];
		return energy;
	}
	
	if (hamil.type == PAULI_TERMS) {
		// if <E> hasn't changed much from prev, previous hamilState is a good estimate
		applyPauliHamil(hamilState, qubits, hamil.termCoeffs, hamil.terms, hamil.numTerms);
		double innerprod = 0;
		for (long long int i=0LL; i < hamil.numAmps; i++)
			innerprod += creal(
				(getRealAmpEl(qubits, i) - I*getImagAmpEl(qubits, i)) *
				hamilState[i]);
		return innerprod;
	}
	
	printf("If you see this message, something has truly crapped up\n");
	return -666;
}



gsl_matrix_complex* getMatrixFromPauliHamil(Hamiltonian hamil, MultiQubit qubits) {
	
	double complex* hamilState = malloc(hamil.numAmps * sizeof *hamilState);
	
	gsl_matrix_complex* matr = gsl_matrix_complex_calloc(hamil.numAmps, hamil.numAmps);
	if (hamil.type == DIAGONAL)
		return matr;
	
	for (long long i=0; i < hamil.numAmps; i++) {
		initClassicalState(&qubits, i);
		applyHamil(hamilState, qubits, hamil);
		
		for (long long j=0; j < hamil.numAmps; j++)
			gsl_matrix_complex_set(matr, i, j, gsl_complex_rect(creal(hamilState[j]), cimag(hamilState[j])));
	}
	
	free(hamilState);

	return matr;
}


void getPauliHamilEigvals(Hamiltonian hamil, MultiQubit qubits, int numModes, double** eigvals, double complex ***eigvecs) {
	
	// passing 0 means all modes are desired
	if (numModes == 0)
		numModes = hamil.numAmps;
	
	// get matrix Hamiltonian
	gsl_matrix_complex* hamilMatr = getMatrixFromPauliHamil(hamil, qubits);
	
	// make GSL space for eigenvals and eigenvectors
	gsl_eigen_hermv_workspace* space = gsl_eigen_hermv_alloc(hamil.numAmps);
	gsl_vector* eigValsVec = gsl_vector_alloc(hamil.numAmps);
	gsl_matrix_complex *eigVecsMatr = gsl_matrix_complex_alloc(hamil.numAmps, hamil.numAmps); 
	
	// diagaonlise Hamiltonian
	int failed = gsl_eigen_hermv(hamilMatr, eigValsVec, eigVecsMatr, space);
	if (failed)
		printf("WTF diagonalisation failed\n");
	
	// sort spectrum by increasing energy
	gsl_eigen_genhermv_sort(eigValsVec, eigVecsMatr, GSL_EIGEN_SORT_VAL_ASC);
	
	// copy from GSL objects to pointers
	*eigvals = malloc(numModes * sizeof **eigvals);
	for (int i=0; i < numModes; i++)
		(*eigvals)[i] = gsl_vector_get(eigValsVec, i);
		
	*eigvecs = malloc(numModes * sizeof **eigvecs);
	for (int i=0; i < numModes; i++) {
		(*eigvecs)[i] = malloc(hamil.numAmps * sizeof ***eigvecs);
		for (int j=0; j < hamil.numAmps; j++) {
			gsl_complex val = gsl_matrix_complex_get(eigVecsMatr, j, i);
			(*eigvecs)[i][j] = GSL_REAL(val) + I*GSL_IMAG(val);
		}
	}
	
	// free GSL objects
	gsl_eigen_hermv_free(space);
	gsl_matrix_complex_free(hamilMatr);
	gsl_vector_free(eigValsVec);
	gsl_matrix_complex_free(eigVecsMatr);
}





void freePauliHamilEigvals(double *eigvals, double complex **eigvecs, int numModes) {
	
	free(eigvals);
	for (int i=0; i < numModes; i++)
		free(eigvecs[i]);
	free(eigvecs);
}