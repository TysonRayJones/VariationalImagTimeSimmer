/** @file 
 * Solves 3SAT problems by Wick rotation
 */

#include <stdio.h>

#include "hamiltonian_builder.h"	// for building 3SATs
#include "param_evolver.h"			// for variational imag-time simulation
#include "mmaformatter.h"			// for outputting results to mathematica

#include "QuEST/qubits.h"		// for simulating the ansatz circuit


#define OUTPUT_FILE "wickSATdata.txt"


double getExpectedEnergy(MultiQubit qubits, double *hamil, long long int stateSize) {
	
	double energy = 0;
	for (long long int i=0LL; i < stateSize; i++)
		energy += getProbEl(qubits, i)*hamil[i];
		
	return energy;
}


double getStableTimeStep(double *hamil, long long int stateSize) {
	
	double maxEigVal = 0;
	for (long long int i=0LL; i < stateSize; i++)
		if (hamil[i] > maxEigVal)
			maxEigVal = hamil[i];
			
	return 1/maxEigVal;
}


int main(int narg, char *varg[]) {
	
	
	/*
	 * GET CMD ARGS 
	 */
	 	 
	if (narg != 10) {
		printf("ERROR! Call with arguments: ");
		printf("num_bools num_params rseed threshold[0 to 1] timestep[0 for auto] max_iters wrap_params ");
		printf("deriv_accuracy[1 to 4] matrix_noise[0 to 1]\n");
		return 1;
	}
	
	int numBools = atoi(varg[1]);
	int numParams = atoi(varg[2]);
	srand(atoi(varg[3]));
	double threshold; sscanf(varg[4], "%lf", &threshold);
	double timeStep; sscanf(varg[5], "%lf", &timeStep);
	int maxIterations = atoi(varg[6]);
	int wrapParams = atoi(varg[7]);
	int derivAccuracy = atoi(varg[8]);
	double matrNoise; sscanf(varg[9], "%lf", &matrNoise);
	
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
	if (timeStep < 0) {
		printf("ERROR! time_step must be positive or 0 for auto");
		return 1;
	}
	
	
	/*
	 * PREPARE SIMULATION
	 */
	
	// generate a random 3SAT problem
	int *equ, *sol;
	double *hamil;
	int numClauses;
	makeRandomEquSolHamil(numBools, &equ, &sol, &hamil, &numClauses);
	
	// get state index of solution (for monitoring progress of QuEST state)
	int solState = 0;
	for (int i=0; i < numBools; i++)
		solState += sol[numBools - i - 1] << i;
	
	printEquSol(equ, sol, numBools, numClauses);
	printf("\nsol ind:\n%d\n\n", solState);
	//printHamil(hamil, numBools);
	
	// choose a sufficiently small time-step
	if (timeStep == 0)
		timeStep = getStableTimeStep(hamil, pow(2, numBools));
	printf("Time step: %lf\n", timeStep);
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numBools, env);
	
	// prepare the param evolver
	evolverMemory mem = prepareEvolverMemory(qubits, numParams);
	
	// set initial param values
	double params[numParams];
	for (int i=0; i < numParams; i++)
		params[i] = 1.0;
		
		
	/*
	 * PREPARE DATA RECORDING
	 */
	
	// prepare records of param values
	double paramEvo[numParams][maxIterations];
	for (int i=0; i < numParams; i++)
		for (int j=0; j < maxIterations; j++)
			paramEvo[i][j] = 0;
	
	// prepare records of sol prob and expected energy
	double solProbEvo[maxIterations];
	double expectedEnergyEvo[maxIterations];
	for (int i=0; i < maxIterations; i++) {
		solProbEvo[i] = -1;
		expectedEnergyEvo[i] = -1;
	}
	
	// prepare records of recovered iterations
	int numRecoveries = 0;
	int recoveredIterations[maxIterations];
	
	
	/*
	 * PERFORM SIMULATION
	 */
	
	// evolve until a prob(sol) threshold
	evolveOutcome outcome;
	int step=0;
	double prob = 0;
	double energy = 0;
	
	while (step < maxIterations && prob < threshold ) {
		
		outcome = evolveParams(
			&mem, defaultAnsatzCircuit, approxParamsByTSVD,
			qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			
		if (outcome == RECOVERED) 
			recoveredIterations[numRecoveries++] = step;
		if (outcome == FAILED) {
			printf("Numerical recovery failed! Aborting entire sim!\n");
			return 1;
		}
			
		prob = getProbEl(qubits, solState);
		energy = getExpectedEnergy(qubits, hamil, mem.stateSize);
		printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step, prob, energy);
		
		solProbEvo[step] = prob;
		expectedEnergyEvo[step] = energy;
		for (int i=0; i < numParams; i++)
			paramEvo[i][step] = params[i];
		
		step += 1;
	}
	
	
	/*
	 * SAVE RESULTS
	 */
	
	// record <E>, prob evolutions, recovered iterations
	FILE* file = openAssocWrite(OUTPUT_FILE);
	writeIntToAssoc(file, "derivAccuracy", derivAccuracy);
	writeDoubleToAssoc(file, "matrNoise", matrNoise, 5);
	writeIntToAssoc(file, "wrapParams", wrapParams);
	writeDoubleToAssoc(file, "timeStep", timeStep, 10);
	writeIntToAssoc(file, "numBools", numBools);
	writeIntToAssoc(file, "numParams", numParams);
	writeDoubleToAssoc(file, "threshold", threshold, 5);
	writeIntToAssoc(file, "convergeTime", step);
	writeDoubleArrToAssoc(file, "solProbEvo", solProbEvo, step, 10);
	writeDoubleArrToAssoc(file, "expectedEnergyEvo", expectedEnergyEvo, step, 10);
	writeIntToAssoc(file, "numRecoveries", numRecoveries);
	writeIntArrToAssoc(file, "recoveredIterations", recoveredIterations, numRecoveries);
	
	// record evolution of every param
	char buf[1024];
	for (int p=0; p < numParams; p++) {
		sprintf(buf, "param%dEvo", p);
		writeDoubleArrToAssoc(file, buf, paramEvo[p], step, 10);
	}

	closeAssocWrite(file);
	
	
	/*
	 * TIDY UP
	 */
	
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