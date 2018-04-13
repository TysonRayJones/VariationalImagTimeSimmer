#ifndef HAMILTONIAN_BUILDER_H_
#define HAMILTONIAN_BUILDER_H_


/**
 * Given a filename, loads a 3SAT equation (setting equ, numBools and numClauses),
 * finds its solution (setting sol) by brute-force, and builds the corresponding
 * hamiltonian (setting hamil) where the energy of a state (bitstring) is the number
 * of clauses in equ that it fails. equ, sol and hamil must be freed!
 */
void loadEquSolHamil(
	char *filename, 
	int **equ, int **sol, double **hamil, int *numBools, int *numClauses
);


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
);

void printHamil(double* hamiltonian, int numBools);

void printEquSol(int *equ, int *sol, int numBools, int numClauses);

void printEquSolHamil(int *equ, int *sol, double *hamil, int numBools, int numClauses);


#endif // HAMILTONIAN_BUILDER_H_
