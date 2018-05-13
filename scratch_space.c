/** @file 
 * Solves 3SAT and Chemistry problems by variational imaginary time propogation
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>

#include "ansatz_circuits.h"
#include "linear_solvers.h"
#include "param_evolver.h"
#include "mmaformatter.h"

#define OUTPUT_FILE "scratch_data.txt"
 
int main(int narg, char *varg[]) {
	
	int rSeed = 0;
	int numIters = 500;
	int numParams = 70;
	int numQubits = 10;
	char* filename = "chemHamil10qb.txt";
	
	double timeStep = 0.1;
	int wrapParams = 0;
	int derivAccuracy = 1;
	double matrNoise = 0;
	
	srand(rSeed);
	Hamiltonian hamil = loadPauliHamilFromFile(filename);
	printHamil(hamil);
	
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numQubits, env);
	
	EvolverMemory mem = prepareEvolverMemory(qubits, numParams);
	initStateZero(&qubits);
	setAnsatzInitState(&mem, qubits);
	
	double params[numParams];
	for (int p=0; p < numParams; p++)
		params[p] = (rand()/(double) RAND_MAX) * 2 * M_PI;
		
	
	for (int step=0; step < numIters; step++) {
		
		evolveOutcome outcome = evolveParams(
			&mem, hardwareEfficientChemistryAnsatzCircuit, approxParamsByTikhonov,
			qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
		if (outcome == FAILED) {
			printf("FAILED!\n");
			return 1;
		}
			
		double energy = getExpectedEnergy(mem.hamilState, qubits, hamil);
		printf("E: %lf\n", energy);
	}
		
		
	freeHamil(hamil);
	freeEvolverMemory(&mem);
	destroyMultiQubit(qubits, env); 
	closeQuESTEnv(env);
}