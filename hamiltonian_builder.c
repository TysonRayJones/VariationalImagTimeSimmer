/** @file 
 * Builds diagonal Hamiltonians from 3SAT problems
 */

#include "sat_generator.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "hamiltonian_builder.h"
 

/**
 * Builds a diagonal hamiltonian from the given 3SAT equation,
 * where the energy of a state (bitstring) is the total number of
 * clauses that it fails. The returned hamiltonian must be freed!
 */
double* getHamil(int* equ, int numBools, int numClauses) {
	
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
 * Prints a diagonal hamiltonian, with format bitstring: energy (as int)
 */
void printHamil(double* hamil, int numBools) {
	
	long long int length = pow(2LL, numBools);
	
	// for every basis vector
	for (long long int i=0LL; i < length; i++) {
		
		// print its bitstring
		for (unsigned int mask = length >> 1; mask; mask >>= 1)
			printf("%d", !!(mask & i));
		
		// print energy
		printf(": %d\n", (int) hamil[i]);
	}
}


/**
 * Given a filename, loads a 3SAT equation (setting equ, numBools and numClauses),
 * finds its solution (setting sol) by brute-force, and builds the corresponding
 * hamiltonian (setting hamil) where the energy of a state (bitstring) is the number
 * of clauses in equ that it fails. equ, sol and hamil must be freed!
 */
void loadEquSolHamil(
	char *filename, 
	int **equ, int **sol, double **hamil, 
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
	*hamil = getHamil(*equ, *numBools, *numClauses);
	
	// equ, sol, hamiltonian must be freed!
}


/**
 * Creates a random 3SAT equation (setting equ and numClauses) of the given size
 * with a unique solution (and other constraints in Simon's paper) (setting sol), 
 * and builds the  corresponding hamiltonian (setting hamil) where the energy of a 
 * state (bitstring) is the number of clauses in equ that it fails. 
 * equ, sol and hamil must be freed!
 */
void makeRandomEquSolHamil(
	int numBools, 
	int **equ, int **sol, double **hamil, int *numClauses
) {
	
	getRandomEquAndSol(numBools, equ, sol, numClauses);
	*hamil = getHamil(*equ, numBools, *numClauses);
	
	// equ, sol hamiltonian must be freed!
}


void printEquSol(int *equ, int *sol, int numBools, int numClauses) {
	
	printf("equ:\n");
	printEquation(equ, numClauses);
	
	printf("\nsol:\n");
	for (int i=0; i<numBools; i++)
		printf("%d", sol[i]);
	printf("\n");
}


void printEquSolHamil(int *equ, int *sol, double *hamil, int numBools, int numClauses) {
	
	printEquSol(equ, sol, numBools, numClauses);
		
	printf("\nhamil:\n");
	printHamil(hamil, numBools);
	printf("\n");
}
