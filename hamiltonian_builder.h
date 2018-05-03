#ifndef HAMILTONIAN_BUILDER_H_
#define HAMILTONIAN_BUILDER_H_

#include <QuEST.h>
#include <complex.h>


typedef enum {DIAGONAL, PAULI_TERMS} hamilType;


typedef struct {
	
	int numQubits;
	long long int numAmps;
	
	// specifies the format of the stored Hamiltonian
	hamilType type;
	
	
	double* diagHamil;
	
	int numTerms;
	double* termCoeffs;
	int** terms;
	
} Hamiltonian;


/**
 * Given a filename, loads a 3SAT equation (setting equ, numBools and numClauses),
 * finds its solution (setting sol) by brute-force, and builds the corresponding
 * diagonal hamiltonian where the energy of a state (bitstring) is the number
 * of clauses in equ that it fails. The returned Hamiltonian must eventually
 * be freed by freeHamil!
 */
Hamiltonian load3SATHamilFromFile(
	char *filename, 
	int **equ, int **sol,
	int *numBools, int *numClauses
);


/**
 * Creates a random 3SAT equation (setting equ and numClauses) of the given size
 * with a unique solution (and other constraints in Simon's paper) (setting sol), 
 * and builds the  corresponding hamiltonian (setting hamil) where the energy of a 
 * state (bitstring) is the number of clauses in equ that it fails. 
 * The returned Hamiltonian must eventually be freed by freeHamil!
 */
Hamiltonian getRandom3SATHamil(
	int numBools, 
	int **equ, int **sol, int *numClauses
);


/**
 * Given a filename, loads a Hamiltonian specified in a file, with format
 * coeff 0 1 0 2 3 0
 * where each line corresponds to a Hamiltonian sum term, coeff is its coefficient
 * and numbers 0,1,2,3 correspond to gates I,X,Y,Z on that index qubit.
 * The returned Hamiltonian must eventually be freed by freeHamil!
 */
Hamiltonian loadPauliHamilFromFile(char *filename);

/**
 * Prints a Hamiltonian. The format depends on the type of Hamiltonian
 * (diagonal or a sum of pauli terms)
 */
void printHamil(Hamiltonian hamil);

/**
 * Modifies hamilState to be the result of applying hamil to qubits.
 * qubits should remain unchanged, though may slightly vary due to numerical imprecision
 */
void applyHamil(double complex* hamilState, MultiQubit qubits, Hamiltonian hamil);

double getExpectedEnergy(double complex* hamilState, MultiQubit qubits, Hamiltonian hamil);

/**
 * Frees the memory allocated to a hamiltonian
 */
void freeHamil(Hamiltonian hamil);

void print3SATEquSol(int *equ, int *sol, int numBools, int numClauses);


#endif // HAMILTONIAN_BUILDER_H_
