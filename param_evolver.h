#ifndef PARAM_EVOLVER_H_
#define PARAM_EVOLVER_H_

#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "../QuEST/qubits.h"


/**
 * Returns the number of parameters needed to repeat the default ancilla circuit
 * numBlocks times, resulting in an equal number of gates (each with an individual parameter)
 * on each qubit
 */
int chooseDefaultNumParams(MultiQubit qubits, int numBlocks);

/**
 * Applies a default parameterised ancilla circuit to the zero state, modifying
 * the wavefunction in qubits according to the values in params. This can be passed
 * to evolveParams in lieu of a custom ancilla circuit
 */
void defaultAncillaCircuit(MultiQubit qubits, double* params, int numParams);

/**
 * A container for the memory used by evolveParams. Is to be created once
 * via prepareEvolverMemory and passed to calls to evolveParams,
 * and finally freed by freeEvolverMemory
 */
typedef struct {
	long long int stateSize;
	int numParams;
	double complex **derivs;
	double complex *hamilState;
	gsl_matrix *matrA;
	gsl_permutation *permA;
	gsl_vector *vecC, *paramChange;
} evolverMemory;

/**
 * Given a list of parameters, a parameterised wavefunction (through ancillaCircuit) and
 * a diagonal Hamiltonian, modifies the parameters by a single time-step under imaginary time 
 * evolution, using Euler's method. defaultAncillaCircuit may be passed in lieu of a custom one.
 * Param evolution is done by repeatedly simulating a parameterised circuit and using finite-
 * difference approximations of derivatives to populate and here solve a family of
 * linear equations of the parameters. This function should update the parameters so that
 * the parameterised wavefunction moves closer to the ground-state of the given Hamiltonian.
 * The output is 0 if updating the parameters is successful, otherwise a 1 is returned which 
 * indicates timeStepSize was too large and a parameter experienced too great a change.
 * Updates the wavefunction in qubits under the new parameters.
 * mem contains memory for matrices and arrays which are modified
 */
int evolveParams(
	evolverMemory *mem, void (*ancillaCircuit)(MultiQubit, double*, int),
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize);

/**
 * Allocates memory for the data structures needed by the evolveParams function,
 * for a given number of blocks. The memory should be kept for the lifetime of the
 * evolver simulation, and afterward freed with freeEvolverMemory
 */
evolverMemory prepareEvolverMemory(MultiQubit qubits, int numBlocks);

/**
 * Frees the memory allocated with prepareEvolverMemory. 
 * This should be done when there are no more calls to evolveParams
 */
void freeEvolverMemory(evolverMemory *memory);



/*
void evolveParamsByGD(
	evolverMemory *mem,
	MultiQubit qubits, double* params, double* diagHamiltonian, double timeStepSize);
*/

#endif // PARAM_EVOLVER_H_