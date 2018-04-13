/** @file 
 * Solves 3SAT problems by Wick rotation
 */

#include <stdio.h>
#include <limits.h> // for rand param start testing

#include "param_evolver.h"
#include "hamiltonian_builder.h"

#include "../QuEST/qubits.h"


double getExpectedEnergy(MultiQubit qubits, double *hamil, long long int stateSize) {
	
	double energy = 0;
	for (long long int i=0LL; i < stateSize; i++)
		energy += getProbEl(qubits, i)*hamil[i];
		
	return energy;
}


int main(int narg, char *varg[]) {
	
	if (narg != 6) {
		printf("ERROR! Call with arguments: ");
		printf("num_bools num_params rseed threshold timestep\n");
		return 1;
	}
	
	// collect cmd params
	int numBools = atoi(varg[1]);
	int numParams = atoi(varg[2]);
	srand(atoi(varg[3]));
	double threshold; sscanf(varg[4], "%lf", &threshold);
	double timeStep; sscanf(varg[5], "%lf", &timeStep);
	
	if (numBools < 4) {
		printf("ERROR! Minimum num_bools is 4\n");
		return 1;
	}
	if (numParams < 1) {
		printf("ERROR! Minimum num_params is 1\n");
		return 1;
	}
	if (threshold <= 0 || threshold >= 1) {
		printf("ERROR! Threshold must be in (0, 1)\n");
		return 1;
	}
	if (timeStep <= 0) {
		printf("ERROR! time_step must be positive");
	}
	
	// generate a random 3SAT problem
	int *equ, *sol;
	double *hamil;
	int numClauses;
	makeRandomEquSolHamil(numBools, &equ, &sol, &hamil, &numClauses);
	
	// get state index of solution
	int solState = 0;
	for (int i=0; i < numBools; i++)
		solState += sol[numBools - i - 1] << i;
	
	printEquSol(equ, sol, numBools, numClauses);
	printf("\nsol ind:\n%d\n\n", solState);
	//printHamil(hamil, numBools);

	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numBools, env);
	
	// prepare the param evolver
	evolverMemory mem = prepareEvolverMemory(qubits, numParams);
	double params[numParams];
	for (int i=0; i < numParams; i++)
		params[i] = 1.0;
		
	// evolve until a prob(sol) threshold
	int paramUnstable = 0;
	int step=0;
	double prob = 0;
	double energy = 0;
	while (prob < threshold) {
		paramUnstable = evolveParams(&mem, defaultAncillaCircuit, qubits, params, hamil, timeStep);
		prob = getProbEl(qubits, solState);
		energy = getExpectedEnergy(qubits, hamil, mem.stateSize);
		printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step++, prob, energy);
		
		if (paramUnstable) {
			timeStep *= 0.9;
			printf("A param changed unstably! Updating timestep to %lf\n", timeStep);
		}
	}
	
	printf("\nfinal param values:\n");
	for (int i=0; i < numParams; i++)
		printf("p%d = %lf\n", i, params[i]);
	
	// cleanup
	freeEvolverMemory(&mem);
	free(equ);
	free(sol);
	free(hamil);
	
	// unload QuEST
	destroyMultiQubit(qubits, env); 
	closeQuESTEnv(env);
	return 0;
}