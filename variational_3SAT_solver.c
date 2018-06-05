/** @file 
 * Solves 3SAT and Chemistry problems by variational imaginary time propogation
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <QuEST.h>					// for simulating ansatz

#include "hamiltonian_builder.h"	// for building and applying Hamiltonians
#include "ansatz_circuits.h"		// for selecting ansatz circuits
#include "linear_solvers.h" 		// for selecting numerical solving method
#include "param_evolver.h"			// for variational imag-time simulation
#include "true_evolver.h"			// for verifying variational sim
#include "mmaformatter.h"			// for outputting results to mathematica

/** whether to record the projection of the state-vec on the true spectrum every iteration 
 * (expensive) */
#define DIAGONALISE_HAMILTONIAN 1
#define NUM_EIGSTATES_TO_SAVE 20

#define RECORD_SPECTRUM_PROJECTIONS 0
#define NUM_EIGSTATES_TO_RECORD_PROJS 15

/** max energy distance beteween non-degen eigvals for compactifying */
#define DEGEN_EIGVAL_MAX_DIST 0.01

/** whether to re-randomise param values when exciting the Hamiltonian */
#define RANDOMISE_PARAMS_WHEN_RESETTING 1

/** threshold for deciding whether a state is ground */
#define ENERGY_TO_GROUND_MAX_DIST 1E-3

/**
 * Initialises a random wavefunction, uniform on the unit hypersphere
 * http://www.qetlab.com/RandomStateVector (k=0)
 * https://en.wikipedia.org/wiki/Normal_distribution#Generating_values_from_normal_distribution
 */
void initStateRandom(MultiQubit *qubits) {
	
	double norm_sq = 0;
	
	// generate unnormalised, normal real and imag components
	for (long long int i=0LL; i < qubits->numAmps; i++) {
		
		double x_uni = rand()/(double) RAND_MAX;
		double y_uni = rand()/(double) RAND_MAX;
		
		// generate pair of normal numbers via Box-Muller
		double x_norm = sqrt(- 2 * log(x_uni)) * cos(2 * M_PI * y_uni);
		double y_norm = sqrt(- 2 * log(x_uni)) * sin(2 * M_PI * y_uni);
		
		qubits->stateVec.real[i] = x_norm;
		qubits->stateVec.imag[i] = y_norm;
		
		norm_sq += pow(x_norm,2) + pow(y_norm,2);
	}
	
	// normalise components
	double scaling = 1/sqrt(norm_sq);
	for (long long int i=0LL; i < qubits->numAmps; i++) {
		qubits->stateVec.real[i] *= scaling;
		qubits->stateVec.imag[i] *= scaling;
	}
}



void getSpectrum(
	Hamiltonian hamilObj, long long int hamilSize, 
	double** spectrum, int** degeneracy, int** stateToSpecMap, long long int* spectrumSize) {
		
	// monkey-patch: abort for non-diagonal Hamiltonians
	if (hamilObj.type != DIAGONAL) {
		printf("Spectrum analysis skipped for non-diagonal Hamiltonian.\n");
		*spectrumSize=0;
		return;
	}
		
		
	double *hamil = hamilObj.diagHamil;
	
	// collect all energies, pops in roomy arrays
	*spectrumSize = 0;
	double energiesContainer[hamilSize];
	int populationContainer[hamilSize];
	*stateToSpecMap = malloc(hamilSize * sizeof **stateToSpecMap);
	
	for (long long int i=0LL; i < hamilSize; i++) {
		double energy = hamil[i];
		
		// find where this energy lives in our arrays
		int located=0;
		for (long long int j=0LL; j < *spectrumSize; j++) {
			if (energiesContainer[j] == energy) {
				populationContainer[j] += 1;
				(*stateToSpecMap)[i] = j;
				located=1;
				break;
			}
		}
		
		// if it's a new energy, add to the arrays
		if (!located) {
			energiesContainer[*spectrumSize] = energy;
			populationContainer[*spectrumSize] = 1;
			(*stateToSpecMap)[i] = *spectrumSize; 
			*spectrumSize += 1;
		}
	}
	
	// copy over arrays to tighter, dynamic ones
	*spectrum = malloc(*spectrumSize * sizeof **spectrum);
	*degeneracy = malloc(*spectrumSize * sizeof **degeneracy);
	for (long long int i=0LL; i < *spectrumSize; i++) {
		(*spectrum)[i] = energiesContainer[i];
		(*degeneracy)[i] = populationContainer[i];
	}
	
}


/*
 * eigvals must previously be sorted 
 */
 
void countSpectrumDegeneracy(double* eigvals, int numEigvals, double** uniqueEigVals, int** degeneracy, int* numUniqueEigVals) {
	
	double myEigVals[numEigvals];
	int myDegeneracy[numEigvals];
	
	int myInd = 0;
	myEigVals[myInd] = eigvals[0];
	myDegeneracy[myInd] = 1;
	
	// find unique energies
	for (int i=1; i < numEigvals; i++) {
		
		if (fabs(eigvals[i] - myEigVals[myInd]) < DEGEN_EIGVAL_MAX_DIST) {
			myDegeneracy[myInd]++;
		}
		else {
			myInd++;
			myEigVals[myInd] = eigvals[i];
			myDegeneracy[myInd] = 1;
		}
	}
	
	// move arrays to dynamic memory
	*numUniqueEigVals = myInd + 1;
	*uniqueEigVals = malloc(*numUniqueEigVals * sizeof **uniqueEigVals);
	*degeneracy = malloc(*numUniqueEigVals * sizeof **degeneracy);
	for (int i=0; i < *numUniqueEigVals; i++) {
		(*uniqueEigVals)[i] = myEigVals[i];
		(*degeneracy)[i] = myDegeneracy[i];
	}
}



int main(int narg, char *varg[]) {
	
	
	
	/*
	 * GET CMD ARGS 
	 */
	
	// list of args
	int numQubits;
	char* hamilTypeString;
	char* hamilFilename;
	int numParams;
	int randSeed;
	double timeStep;
	int maxIterations;
	int wrapParams;
	int derivAccuracy;
	double matrNoise;
	int useGradDesc;
	char* ansatzCircString;
	double paramChangeThreshold;
	int exciteWhenStuck;
	int simRepetitions;
	int progressPrintFrequency;
	char* outputFilename;
	
	// interactive arg input
	if (narg == 2 && !strcmp(varg[1], "help")) {
		
		printf("\n");
		printf("Enter the number of qubits (or booleans in the 3SAT equation)\n(must be >= 4 for 3SAT)\n");
		printf("num_qubits: ");
		scanf("%d", &numQubits);
		printf("\n");
		
		// TODO: get type
		printf("Enter the type of the Hamiltonian, '3SAT' or 'Chemistry'\n");
		printf("hamil_type: ");
		hamilTypeString = "3SAT";
		printf("\nUSING DEFAULT '3SAT' (until I fix this)\n");
		
		// TODO: get filename
		printf("Enter the filename of the Hamiltonian, or 'auto' to generate a random Hamiltonian\n");
		printf("hamil_filename: ");
		hamilFilename="auto";
		printf("\nUSING DEFAULT 'auto' (until I fix this)\n");
		
		printf("Enter the number of parameters total in the ansatz circuit\n");
		printf("num_params: ");
		scanf("%d", &numParams);
		printf("\n");
		
		printf("Enter the random seed for generating 3SATs, initial param values, matrix noise " 
			   "and GSL numerics\n");
		printf("rseed: ");
		scanf("%d", &randSeed);
		printf("\n");
		
		printf("Enter the time-step (float), or 0 for auto\n");
		printf("timestep: ");
		scanf("%lf", &timeStep);
		printf("\n");
		
		printf("Enter the maximum number of variational iterations in each simulation\n");
		printf("max_iters: ");
		scanf("%d", &maxIterations);
		printf("\n");
		
		printf("Enter whether to wrap-around params to keep them in [0, 2PI)\n"
			   "(1 for yes, 0 for no)\n");
		printf("wrap_params: ");
		scanf("%d", &wrapParams);
		printf("\n");
		
		printf("Enter the accuracy of the derivative estimates, as a finite-difference order\n"
			   "(1 for fastest, 4 for most accurate)\n");
		printf("deriv_accuracy: ");
		scanf("%d", &derivAccuracy);
		printf("\n");
		
		printf("Enter the fraction of random noise in the A and C matrices each iteration\n"
			   "(0 for no noise, e.g. 0.3 for maximum +- 30%% random fluctuation in each element)\n");
		printf("matrix_noise: ");
		scanf("%lf", &matrNoise);
		printf("\n");
		
		printf("Enter whether to use imaginary-time evolution (0) or gradient descent (1)\n");
		printf("use_gd: ");
		scanf("%d", &useGradDesc);
		printf("\n");
		
		// TODO: get ansatz
		printf("Enter ansatz circuit\n");
		printf("ansatz: ");
		ansatzCircString = "lowDepth";
		printf("\nUSING DEFAULT lowDepth (until I fix this)\n\n");
		
		printf("Enter param change threshold, below which, the state is assumed converged\n");
		printf("param_change: ");
		scanf("%lf", &paramChangeThreshold);
		printf("\n");
		
		printf("Enter whether to excite out of stuck states (Suguru's method)\n");
		printf("0 for never, 1 for always, 2 for until ground is reached\n");
		printf("excite_when_stuck: ");
		scanf("%d", &exciteWhenStuck);
		printf("\n");
		
		printf("Enter the number of times to resimulate the given system with different "
			   "random initial parameters\n");
		printf("sim_reps: ");
		scanf("%d", &simRepetitions);
		printf("\n");
		
		printf("Enter how frequently convergence progress should be printed, which slows execution\n"
			   "(0 for never, 1 for every iteration, n for every n iterations)\n");
		printf("print_progress_every: ");
		scanf("%d", &progressPrintFrequency);
		printf("\n");
		
		// TODO: get output filename
		printf("Enter the output filename\n");
		printf("output_fn: ");
		outputFilename = "VarData.txt";
		printf("\nUSING DEFAULT VarData.txt (until I fix this)\n\n");
		
		
	// invalid number of args input
	} else if (narg != 18) {
		printf("\nERROR! Call with arguments:\n");
		printf(
			"num_bools\n"
			"hamil_type['3SAT' or 'Chemistry']\n"
			"hamil_filename['auto' for new, random Hamiltonian]\n"
			"num_params[0 for auto]\n"
			"rseed\n"
			"timestep[0 for auto]\n"
			"max_iters\n"
			"wrap_params\n"
			"deriv_accuracy[1 to 4]\n"
			"matrix_noise[0 to 1]\n"
			"use_gd[0 or 1]\n"
			"ansatz\n"
			"param_change\n"
			"excite_when_stuck[0 never, 1 always, 2 only when in excited]\n"
			"sim_reps\n"
			"print_progress_every\n"
			"output_fn\n"
			"\n"
			"Run './Variational3SATSolver help' to enter arguments interactively\n\n");
			
		return 1;
		
	// cmd arg input
	} else {
	
		numQubits = atoi(varg[1]);
		hamilTypeString = varg[2];
		hamilFilename = varg[3];
		numParams = atoi(varg[4]);
		randSeed = atoi(varg[5]); 
		sscanf(varg[6], "%lf", &timeStep);
		maxIterations = atoi(varg[7]);
		wrapParams = atoi(varg[8]);
		derivAccuracy = atoi(varg[9]);
		sscanf(varg[10], "%lf", &matrNoise);
		useGradDesc = atoi(varg[11]);
		ansatzCircString = varg[12];
		sscanf(varg[13], "%lf", &paramChangeThreshold);
		exciteWhenStuck = atoi(varg[14]);
		simRepetitions = atoi(varg[15]);
		progressPrintFrequency = atoi(varg[16]);
		outputFilename = varg[17];
	}
	
	// seed before 3SAT Hamil generation
	srand(randSeed);
	
	// invalid arg values
	if (numQubits < 4) {
		printf("ERROR! Minimum num_bools is 4\n");
		return 1;
	}
	if (numParams < 0) {
		printf("ERROR! Minimum num_params is 1, or 0 for auto\n");
		return 1;
	}
	if (timeStep < 0) {
		printf("ERROR! time_step must be positive or 0 for auto\n");
		return 1;
	}
	if (simRepetitions < 1) {
		printf("ERROR! sim_reps must be greater than 0\n");
		return 1;
	}
	if (progressPrintFrequency < 0) {
		printf("ERROR! print_progress_every must be positive (1 for every iteration, 0 for never)\n");
		return 1;
	}
	if (exciteWhenStuck == 2 && DIAGONALISE_HAMILTONIAN == 0 && strcmp(hamilTypeString, "Chemistry") == 0) {
		printf("ERROR! exciteWhenStucK=2 (excite until ground is reached) but DIAGONALISE_HAMILTONIAN=0, ");
		printf("so the ground state energy isn't known\n");
		return 1;
	}
	if (DIAGONALISE_HAMILTONIAN == 0 && RECORD_SPECTRUM_PROJECTIONS == 1) {
		printf("ERROR! RECORD_SPECTRUM_PROJECTIONS=1 but DIAGONALISE_HAMILTONIAN=0\n");
		return 1;
	}
	if (timeStep == 0 && DIAGONALISE_HAMILTONIAN == 0 && strcmp(hamilTypeString, "Chemistry")==0) {
		printf("ERROR! Chemistry Hamiltonians can only give auto time-step if diagonalised!\n");
		return 1;
	}
	AnsatzCircuit ansatz = getAnsatzFromString(ansatzCircString);
	if (ansatz == NULL) {
		printf("ERROR! Unrecognised ansatz circuit '%s'\n", ansatzCircString);
		return 1;
	}
	if (strcmp(hamilTypeString, "3SAT") != 0 && strcmp(hamilTypeString, "Chemistry") != 0) {
		printf("ERROR! hamil_type must be '3SAT' or 'Chemistry'\n");
		return 1;
	}
	
	// apply (some of the) auto-parameters
	if (numParams == 0)
		numParams = getIdealNumParamsInAnsatz(ansatz, numQubits);
	
	// confirm args
	printf("\n");
	printf("numQubits: %d\n", numQubits);
	printf("hamilType: %s\n", hamilTypeString);
	printf("hamilFilename: %s\n", hamilFilename);
	printf("numParams: %d\n", numParams);
	printf("randSeed: %d\n", randSeed);
	if (timeStep == 0)
		printf("timeStep: TBC\n");
	else
		printf("timeStep: %lf\n", timeStep);
	printf("iterations: %d\n", maxIterations);
	printf("wrapParams: %d\n", wrapParams);
	printf("derivAccuracy: %d\n", derivAccuracy);
	printf("matrNoise: %lf\n", matrNoise);
	printf("useGradDesc: %d\n", useGradDesc);
	printf("ansatzCircuit: %s\n", ansatzCircString);
	printf("paramChangeThreshold: %lf\n", paramChangeThreshold);
	printf("exciteWhenStuck: %d\n", exciteWhenStuck);
	printf("simRepettions: %d\n", simRepetitions);
	printf("printEvery: %d\n", progressPrintFrequency);
	printf("outputFilename: %s\n", outputFilename);
	printf("\n");
	
	
		
	// record results
	printf("Testing access to output file ..\n");
	FILE* testfile = openAssocWrite(outputFilename);
	writeStringToAssoc(testfile, "STATUS", "RUNNING...");
	closeAssocAppend(testfile);
	
	
	/*
	 * PREPARE HAMILTONIAN
	 */
	 
	// 3SAT structures unused if hamilType = 'Chemistry'
	int *equ = NULL;
	int *sol = NULL;
	int numClauses = -1;
	 
	printf("Preparing Hamiltonian...\n");
	Hamiltonian hamil;
	
	// prepare 3SAT Hamiltonian
	if (strcmp(hamilTypeString, "3SAT") == 0) {
		
		// generate random 3SAT
		if (strcmp(hamilFilename, "auto") == 0)
			hamil = getRandom3SATAndHamil(numQubits, &equ, &sol, &numClauses);
			
		// load 3SAT from file
		else {
			int numBools = -1;
			int failed = 0;
			hamil = load3SATAndHamilFromFile(hamilFilename, &equ, &sol, &numBools, &numClauses, &failed);
			if (failed) {
				printf("ERROR! 3SAT Hamiltonian file (%s) not found!\n", hamilFilename);
				return 1;
			}
			if (numBools != numQubits) {
				printf(
					"ERROR! User specified %d qubits, but supplied 3SAT Hamiltonian (%s) had %d booleans\n",
					numQubits, hamilFilename, numBools);
				free(equ);
				free(sol);
				return 1;
			}
		}
	}
	
	// prepare Chemistry Hamiltonian
	if (strcmp(hamilTypeString, "Chemistry") == 0) {
		
		// generate random Chemistry problem
		if (strcmp(hamilFilename, "auto") == 0) {
			// TODO: implement this
			printf("ERROR! Random chemistry Hamiltonians are not yet supported!\n");
			return 1;
			
		// load Chemistry problem from file
		} else {
			int failed = 0;
			hamil = loadPauliHamilFromFile(hamilFilename, &failed);
			if (failed) {
				printf("ERROR! Chemistry Hamiltonian file (%s) not found!\n", hamilFilename);
				return 1;
			}
		}
	}
	
	
	
	/*
	 * PREPARE SIMULATOR
	 */
	 
	printf("Preparing simulator...\n");
	
	// prepare QuEST
	QuESTEnv env;
	initQuESTEnv(&env);
	MultiQubit qubits; 
	createMultiQubit(&qubits, numQubits, env);
	
	// allocate memory for param evolver
	EvolverMemory mem = prepareEvolverMemory(qubits, numParams);
	
	
	
	/*
	 * STUDY SPECTRUM
	 */
	 
	// Chemistry spectra
	double* eigvals = NULL;
	double complex** eigvecs = NULL;
	int numUniqEigvals = -1;
	double* uniqEigvals = NULL;
	int* eigvalDegen = NULL;
	
	// lowest energy state
	double groundStateEnergy = -666;
	
	// find Chemistry spectrum by diagonalisation
	if (hamil.type == PAULI_TERMS && DIAGONALISE_HAMILTONIAN) {
		
		printf("Diagonalising chemistry Hamiltonian...\n");
		getPauliHamilEigvals(hamil, qubits, hamil.numAmps, &eigvals, &eigvecs);
		printf("Compacting chemistry spectrum...\n");
		countSpectrumDegeneracy(eigvals, hamil.numAmps, &uniqEigvals, &eigvalDegen, &numUniqEigvals);
		
		groundStateEnergy = eigvals[0];
	}
	
	// 3SAT spectrum is already diagonal
	if (hamil.type == DIAGONAL) {
		
		groundStateEnergy = 0;
	}

	// get state index of solution (for monitoring progress of QuEST state)
	/*
	long long int solState = 0;
	for (int i=0; i < numQubits; i++)
		solState += sol[numQubits - i - 1] << i;
	*/
	
	// choose a sufficiently small time-step
	
	if (timeStep == 0) {

		if (hamil.type == DIAGONAL)
			timeStep = getStableImagTimestep(hamil);
		if (hamil.type == PAULI_TERMS)
			timeStep = 1/uniqEigvals[numUniqEigvals-1];
			
		printf("Automatic time-step: %lf\n", timeStep);
	}
	
	
	
	/*
	 * PREPARE DATA RECORDING
	 */
	 
	printf("Initialising data structures...\n");
	 
	//double samsInitParams[42] ={3.243044, 4.81527666, 0.62511592, 0.41921514, 4.24673659, 3.43911744, 0.0877967, 6.06145306, 1.08071278, 6.10536731, 4.16570543, 2.37626469, 5.71440099, 2.77030973, 0.46173575, 3.36354362, 4.70542023, 6.11022712, 2.45688627, -0.11424261, 3.69617933, 1.098452, 4.71331006, 0.99540201, 1.01557048, -0.05416214, 1.13185388, 4.36649297, 0.44485088, 1.89563555, 0.847721051, 2.44248901, 2.76156372, 3.24083925, 0.78477016, 5.50861965, 0.21425296, 2.59175576, 2.44989661, 4.19120921, 1.88714856, 3.23731343};
	//double samsFinalParams[42] = {-10.85422203, 87.96330825, 1.16092286E-05, -0.05925503, -3.08235636, 0.00205020046, -18.84501712, 45.39287869, 3.54246012, 51.8369899, 35.95813837, -32.98477717, -1011.59145386, -3.39717838, -18.78717354, 3.22934895, 3.14376364, 8.69201468, 15.70342282, 74.22959909, 11.0544325, -48.36063456, -221.26090968, 45.59825606, 13.99581773, -0.000161766399, 3.14141416, 0.000886289673, 15.71206725, 3.14159531, 54.64032903, -7.85421824, 4.77108247, -1.7355999, -199.49040412, -12.10106662, 3.14141678, 6.2835601, 6.27902231, 2.39810243, 3.14181229, 0.00264922488};
	
	
	// pre-allocate initial param values (GSL messes with seeding)
	double params[numParams];
	double initParams[simRepetitions][numParams];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < numParams; i++)
			initParams[s][i] = (rand()/(double) RAND_MAX) * 2 * M_PI; // /samsInitParams[i]; // M_PI; //0.005;    ///////////// CHANGE THIS BACCCCCCKKKKK
	
	// prepare records of param values
	double*** paramEvo = malloc(simRepetitions * sizeof *paramEvo);
	for (int i=0; i < simRepetitions; i++) {
		paramEvo[i] = malloc(numParams * sizeof **paramEvo);
		for (int j=0; j < numParams; j++) {
			paramEvo[i][j] = malloc(maxIterations * sizeof ***paramEvo);
			for (int k=0; k < maxIterations; k++)
				paramEvo[i][j][k] = -666;
		}
	}
	
	/*
	double paramEvo[simRepetitions][numParams][maxIterations];
	for (int i=0; i < simRepetitions; i++)
		for (int j=0; j < numParams; j++)
			for (int k=0; k < maxIterations; k++)
				paramEvo[i][j][k] = -666;
	*/
	
	
	// prepare records expected energy...
	double expectedEnergyEvo[simRepetitions][maxIterations+1];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < maxIterations+1; i++)
			expectedEnergyEvo[s][i] = -1;
	
	// and solution prob (only relevant for DIAGONAL Hamiltonians)
	/*
	double solProbEvo[simRepetitions][maxIterations];
	for (int s=0; s < simRepetitions; s++)
		for (int i=0; i < maxIterations; i++)
			solProbEvo[s][i] = -1;
	*/
			
	// prepare records of discovered eigenstates (only relevant if exciteWhenStuck=1)
	int numExcitationsFound[simRepetitions];
	int itersOfExcitationsFound[simRepetitions][hamil.numAmps];
	double energiesOfExcitationsFound[simRepetitions][hamil.numAmps];
	for (int s=0; s < simRepetitions; s++) {
		numExcitationsFound[s] = 0;
		for (int n = 0; n < hamil.numAmps; n++) {
			itersOfExcitationsFound[s][n] = -1;
			energiesOfExcitationsFound[s][n] = -1;
		}
	}
	
	// prepare records of projections of pauli_term eigstates on state-vector
	double*** eigstateProbsEvo = malloc(simRepetitions * sizeof *eigstateProbsEvo);
	for (int s=0; s < simRepetitions; s++) {
		eigstateProbsEvo[s] = malloc(NUM_EIGSTATES_TO_RECORD_PROJS * sizeof **eigstateProbsEvo);
		for (int n=0; n < NUM_EIGSTATES_TO_RECORD_PROJS; n++) {
			eigstateProbsEvo[s][n] = malloc(maxIterations * sizeof ***eigstateProbsEvo);
			for (int i=0; i < maxIterations; i++)
				eigstateProbsEvo[s][n][i] = -1;
		}
	}
	
	// analyse spectrum (only valid for diagonal hamiltonians)
	
	printf("Compacting 3SAT spectrum...\n");
	
	double* spectrum;
	int* degeneracy;
	int* stateToSpecMap;
	long long int spectrumSize;
	long long int solStateSpecInd = 0;
	getSpectrum(
		hamil, qubits.numAmps, 
		&spectrum, &degeneracy, &stateToSpecMap, &spectrumSize
	);
	double specProbEvo[simRepetitions][spectrumSize][maxIterations];
		
	if (hamil.type == DIAGONAL) {
		printf("\nSpectrum size:\t%lld\n", spectrumSize);
		for (long long int i=0LL; i < spectrumSize; i++)
			printf("energy:\t%lf,\tdegeneracy:\t%d\n", spectrum[i], degeneracy[i]);
		for (long long int i=0LL; i < spectrumSize; i++)
			if (spectrum[i] == 0)
				solStateSpecInd = i;
		printf("spectrum[%lld] = 0\n\n", solStateSpecInd);
		
		// prepare records of spectrum evolution
		for (int r=0; r < simRepetitions; r++)
			for (long long int i=0LL; i < spectrumSize; i++)
				for (int j=0; j < maxIterations; j++)
					specProbEvo[r][i][j] = -666;
	}
	
	
	
	/*
	 * PERFORM SIMULATION
	 */
		
	// set the state we'll feed into the ansatz
	initStateZero(&qubits);
	// sigmaX(qubits, 0);    ///////////////////////////////////////////// NOTTED FOR HARTREE FOCK
	setAnsatzInitState(&mem, qubits);
	
	evolveOutcome outcome;
	double energy;
	
	// re-simulate many times
	for (int rep=0; rep < simRepetitions; rep++) {
		
		printf("Beginning simulation %d:\n", rep);
			
		// set random initial param values
		for (int i=0; i < numParams; i++)
			params[i] = initParams[rep][i];
		
		// remove any static states from the Hamiltonian excitations
		clearExcitedStates(&mem);
		
		// get energy of initial param values
		ansatz(&mem, qubits, params, numParams);
		energy = getExpectedEnergy(mem.hamilState, qubits, hamil);
		expectedEnergyEvo[rep][0] = energy;
		printf("initial\t<E> = %f\n", energy);
		
		// keep evolving until we converge or reach max iterations
		for (int step = 0; step < maxIterations; step++) {
			
			// update params under parameterised evolution
			if (useGradDesc)
				outcome = evolveParamsByGradientDescent(
					&mem, ansatz, qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			else
				outcome = evolveParamsByImaginaryTime(
					&mem, ansatz, 
					//approxParamsByTSVD,
					approxParamsByTikhonov,
					//approxParamsByLUDecomp,
					qubits, params, hamil, timeStep, wrapParams, derivAccuracy, matrNoise);
			
			if (outcome == FAILED) {
				printf("Numerical inversion failed! Aborting entire sim!\n");
				return 1;
			}
			
			// evolve wavefunction (non-parameterised) under exact evolution
			//evolveWavefunction(qubits, hamil, timeStep);
			
			// monitor convergence
			
			//prob = 0; // getProbEl(qubits, solState);
			energy = getExpectedEnergy(mem.hamilState, qubits, hamil);
		
			// print progress
			if (progressPrintFrequency != 0 && step % progressPrintFrequency == 0) {
				//if (hamil.type == DIAGONAL)
				//	printf("t%d: \t prob(sol) = %f \t <E> = %f\n", step, prob, energy);
				//else
				printf("t%d: \t<E> = %f\n", step, energy);
			}
			
			
			//solProbEvo[rep][step] = prob;
			
			// record expected energy
			expectedEnergyEvo[rep][step+1] = energy;
			
			// record params
			for (int i=0; i < numParams; i++)
				paramEvo[rep][i][step] = params[i];
				
			// record spectrum evo data
			if (hamil.type == DIAGONAL) {
				for (int i=0; i < spectrumSize; i++)
					specProbEvo[rep][i][step] = 0;
				for (long long int i=0LL; i < qubits.numAmps; i++)
					specProbEvo[rep][stateToSpecMap[i]][step] += getProbEl(qubits, i);
			}
			
			//record spectrum projection data
			if (RECORD_SPECTRUM_PROJECTIONS && hamil.type == PAULI_TERMS) {
				for (int n=0; n < NUM_EIGSTATES_TO_RECORD_PROJS; n++) {
					double complex eigProj = innerProductOnQubits(eigvecs[n], qubits, qubits.numAmps);
					double eigProb = pow(creal(eigProj),2) + pow(cimag(eigProj),2);
					eigstateProbsEvo[rep][n][step] = eigProb;
				}
			}
			
			// randomly wiggle params
			/*
			if (step > 0 && step % 50 == 0) {
				printf("Wiggling params!");
				
				for (int i=0; i < numParams; i++)
					params[i] += 0.01 * 2*M_PI*(rand() / (double) RAND_MAX);
			}
			*/
			
			
			if (
				(exciteWhenStuck==1 || (    // always excite converged state
					exciteWhenStuck==2 &&   // only excite non-ground states
					fabs(energy - groundStateEnergy) > ENERGY_TO_GROUND_MAX_DIST)
				) 
				&& 
				isStuck(paramEvo, rep, numParams, step, paramChangeThreshold)
			) {
				
				// excite the state in the Hamiltonian
				printf("Exciting state:\n");
				exciteStateInHamiltonian(&mem, qubits);
				
				// record the excitation
				itersOfExcitationsFound[rep][numExcitationsFound[rep]] = step;
				energiesOfExcitationsFound[rep][numExcitationsFound[rep]] = energy;
				numExcitationsFound[rep]++;
				
				// reset the parameters
				if (RANDOMISE_PARAMS_WHEN_RESETTING) {
					printf("Assigning new random params!\n");
					for (int n=0; n < numParams; n++)
						params[n] = (rand()/(double) RAND_MAX) * 2 * M_PI;
				} else {
					printf("Restoring initial params!\n");
					for (int n=0; n < numParams; n++)
						params[n] = initParams[rep][n];
				}
			}
		}
	
		/*
		 * SAVE RESULTS TO FILE EVERY RE-SIMULATION
		 */
		
		printf("Writing to file...\n");
		
		// record results
		FILE* file = openAssocWrite(outputFilename);
		
		// meta-data
		writeIntToAssoc(file, "numQubits", numQubits);
		writeStringToAssoc(file, "hamilFilename", hamilFilename);
		writeIntToAssoc(file, "numParams", numParams);
		writeIntToAssoc(file, "randSeed", randSeed);
		writeIntToAssoc(file, "simRepetitions", simRepetitions);
		writeIntToAssoc(file, "maxIterations", maxIterations);
		writeIntToAssoc(file, "derivAccuracy", derivAccuracy);
		writeDoubleToAssoc(file, "matrNoise", matrNoise, 5);
		writeIntToAssoc(file, "wrapParams", wrapParams);
		writeDoubleToAssoc(file, "timeStep", timeStep, 10);
		writeIntToAssoc(file, "useGradDesc", useGradDesc);
		writeStringToAssoc(file, "ansatz", ansatzCircString);
		writeDoubleToAssoc(file, "paramChangeThreshold", paramChangeThreshold, 10);
		writeIntToAssoc(file, "exciteWhenStuck", exciteWhenStuck);
		
		writeIntToAssoc(file, "RANDOMISE_PARAMS_WHEN_RESETTING", RANDOMISE_PARAMS_WHEN_RESETTING);
		writeIntToAssoc(file, "NUM_ITERS_IN_STUCK_CHECK", NUM_ITERS_IN_STUCK_CHECK);
		writeDoubleToAssoc(file, "DERIV_STEP_SIZE", DERIV_STEP_SIZE, 15);
		writeIntToAssoc(file, "MAX_NUM_SAVED_STATES", MAX_NUM_SAVED_STATES);
		writeDoubleToAssoc(file, "EXCITATION_OF_SAVED_STATES", EXCITATION_OF_SAVED_STATES, 10);
		writeIntToAssoc(file, "ENERGY_TO_GROUND_MAX_DIST", ENERGY_TO_GROUND_MAX_DIST);
		writeIntToAssoc(file, "DIAGONALISE_HAMILTONIAN", DIAGONALISE_HAMILTONIAN);
		
		// solution, energy, param evolution
		writeNestedDoubleArrToAssoc(file, "expectedEnergyEvos", expectedEnergyEvo, 2, (int []) {simRepetitions, maxIterations+1}, maxIterations+1, 10);
		writeNestedDoubleArrToAssoc(file, "initParams", initParams, 2, (int []) {simRepetitions, numParams}, numParams, 10);
		//writeNestedDoubleArrToAssoc(file, "paramEvos", paramEvo, 3, (int []) {simRepetitions, numParams, maxIterations}, maxIterations, 10);
		writeNestedDoubleListToAssoc(file, "paramEvos", paramEvo, 3,  (int []) {simRepetitions, numParams, maxIterations}, 10);
		
		// only include spectrum data relevant to the Hamiltonian type
		if (hamil.type == PAULI_TERMS && DIAGONALISE_HAMILTONIAN == 1) {
			
			writeStringToAssoc(file, "hamilType", "PAULI_TERMS");
			
			// write spectrum
			int numToSave = (NUM_EIGSTATES_TO_SAVE < hamil.numAmps)? NUM_EIGSTATES_TO_SAVE : hamil.numAmps;
			writeIntToAssoc(file, "spectrumSize", numToSave);
			writeDoubleArrToAssoc(file, "spectrum", eigvals, numToSave, 10);
			
			numToSave = (NUM_EIGSTATES_TO_SAVE < numUniqEigvals)? NUM_EIGSTATES_TO_SAVE : numUniqEigvals;
			writeIntToAssoc(file, "uniqueSpectrumSize", numToSave);
			writeDoubleArrToAssoc(file, "uniqueSpectrum", uniqEigvals, numToSave, 10);
			writeIntArrToAssoc(file, "uniqueSpectrumDegeneracy", eigvalDegen, numToSave);
		
		}
		if (hamil.type == DIAGONAL) {
			
			// 3SAT equ
			writeStringToAssoc(file, "hamilType", "DIAGONAL");
			writeIntArrToAssoc(file, "3SATEqu", equ, numClauses*3);
			writeIntArrToAssoc(file, "3SATSol", sol, numQubits);
			
			// only 3SATs have explicit single solution
			//writeNestedDoubleArrToAssoc(file, "solProbEvos", solProbEvo, 2, (int []) {simRepetitions, maxIterations}, maxIterations, 10);
			
			// spectrum evolution
			writeIntToAssoc(file, "spectrumSize", spectrumSize);
			writeDoubleArrToAssoc(file, "spectrum", spectrum, spectrumSize, 1);
			writeIntArrToAssoc(file, "spectrumDegeneracy", degeneracy, spectrumSize);
			writeIntToAssoc(file, "solStateSpecInd", solStateSpecInd);
			writeNestedDoubleArrToAssoc(file, "specEvos", specProbEvo, 3, (int []) {simRepetitions, spectrumSize, maxIterations}, maxIterations, 10);
		}
		
		// only write discovered excitation data if populated
		if (exciteWhenStuck) {
			writeUnevenOnceNestedIntArrToAssoc(file, "itersOfExcitationsFound", itersOfExcitationsFound, simRepetitions, numExcitationsFound, hamil.numAmps);
			writeUnevenOnceNestedDoubleArrToAssoc(file, "energiesOfExcitationsFound", energiesOfExcitationsFound, simRepetitions, numExcitationsFound, hamil.numAmps, 10);
			writeIntArrToAssoc(file, "numExcitationsFound", numExcitationsFound, simRepetitions);
		}
		
		if (RECORD_SPECTRUM_PROJECTIONS && hamil.type == PAULI_TERMS) {
			writeIntToAssoc(file, "numSpectrumProbEvos", NUM_EIGSTATES_TO_RECORD_PROJS);
			writeNestedDoubleListToAssoc(file, "spectrumProbEvos", eigstateProbsEvo, 3, 
				(int []) {simRepetitions, NUM_EIGSTATES_TO_RECORD_PROJS, maxIterations}, 10);
		}

		// finish writing data this repetition
		closeAssocWrite(file);
	
	}
	
	
	
	/*
	 * TIDY UP
	 */
	 
	// free param evolution
	for (int i=0; i < simRepetitions; i++) {
		for (int j=0; j < numParams; j++)
			free(paramEvo[i][j]);
		free(paramEvo[i]);
	}
	free(paramEvo);
	
	// free eigstate evolution
	for (int s=0; s < simRepetitions; s++) {
		for (int n=0; n < NUM_EIGSTATES_TO_RECORD_PROJS; n++)
			free(eigstateProbsEvo[s][n]);
		free(eigstateProbsEvo[s]);
	}
	
	if (hamil.type == PAULI_TERMS) {
		freePauliHamilEigvals(eigvals, eigvecs, hamil.numAmps);
	}
	if (hamil.type == DIAGONAL) {
		free(equ);
		free(sol);
		free(spectrum);
		free(degeneracy);
		free(stateToSpecMap);
	}
	
	freeHamil(hamil);
	freeEvolverMemory(&mem);
	destroyMultiQubit(qubits, env); 
	closeQuESTEnv(env);
	
	return 0;
}