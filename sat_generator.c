/** @file 
 * Builds, loads and solves 3SAT problems
 */
 
 #include "sat_generator.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CLAUSE_SIZE 3
#define CLAUSE_BOOL_RATIO 4.267


/*
 * speedups (18 Nov):
 * - reordered while loop conds (iterate lists last)
 * - replace qsort with hard-coded CLAUSE_SIZE=3 sort
 * - reduced number of pow(2, numBools) calls to 1
 * - reduced reserving space for candidate solutions to once
 */


int *loadEquation(char *filename, int *numBools, int *numClauses, int *failed) {
	/*
	 * reads N-SAT equation from file, indexed-from-1
	 * mallocs an array; must be freed by caller.
	 */
		
	*numBools = 0;
	*numClauses = 0;
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		*failed=1;
		return NULL;
	}
	
	
	// count parameters in equation
	int ind = 0, qb = 0;
	do {
		fscanf(file, "%d", &qb);
		if (abs(qb) > *numBools)
			*numBools = abs(qb);
		if (++ind == CLAUSE_SIZE)
			*numClauses += 1;
		ind %= CLAUSE_SIZE;
	} while (!feof(file));
	
	// create equation structure
	rewind(file);
	int *equ = malloc(*numClauses * CLAUSE_SIZE * sizeof(*equ));
	for (ind=0; ind < *numClauses * CLAUSE_SIZE; ind++)
		fscanf(file, "%d", &equ[ind]);
		
	fclose(file);
	*failed=0;
	return equ;
}


void saveEquation(char *filename, int *equ, int numClauses) {
	
	FILE *file = fopen(filename, "w");
	
	for (int clause=0; clause < numClauses; clause++) {
		for (int term=0; term < CLAUSE_SIZE; term++) {
			int ind = clause*CLAUSE_SIZE + term;
			
			// avoid trailing space on each line
			if (term == CLAUSE_SIZE-1)
				fprintf(file, "%d", equ[ind]);
			else
				fprintf(file, "%d ", equ[ind]);
		}
		
		// avoid trailing newline at end of file
		if (clause < numClauses - 1)
			fprintf(file, "\n");
	}
	
	fclose(file);
}


void saveSolution(char *filename, int *sol, int numBools) {
	
	FILE *file = fopen(filename, "w");
	
	for (int i=0; i < numBools; i++)
		fprintf(file, "%d ", sol[i]);
		
	fclose(file);
}


void printEquation(int *equ, int numClauses) {
	
	int ind=0;
	for (int i=0; i < numClauses; i++) {
		for (int j=0; j < CLAUSE_SIZE; j++)
			printf("%d ", equ[ind++]);
		printf("\n");
	}
}


void convertToBinary(int number, int numBits, int *outputBits) {
	// reversed, ofc
	int bit;
	for (int i=0; i < numBits; i++) {
		bit = number % 2;
		number /= 2;
		outputBits[i] = bit;
	}
}


int isValidSolution(int *equ, int *sol, int numClauses) {
	
	int clauseInd; 	// index of the clause being checked
	int localInd;	// iterates 0 to CLAUSE_SIZE
	int termInd;	// index of the considered term in equ
	int termQb;		// index of the considered term's qubit in register
	int termBool;	// 0 if considered term appears negated, else 1
	
	// check each clause
	for (clauseInd=0; clauseInd<numClauses; clauseInd++) {
		
		// whether the current clause is satisfied by sol
		int passedClause = 0;
		
		// check each bool/qubit in clause
		for (localInd=0; localInd<CLAUSE_SIZE; localInd++) {
			
			// locate qubit
			termInd = clauseInd*CLAUSE_SIZE + localInd;
			termQb = abs(equ[termInd]) - 1;
			termBool = (equ[termInd] > 0);
			
			// check if qubit satisfies clause
			if (sol[termQb] == termBool) {
				passedClause = 1;
				break;
			}
		}
		
		// if clause isn't satisfied, solution is invalid
		if (!passedClause)
			return 0;
	}
	
	// if all clauses aere satisfied, solution is valid
	return 1;
}


void nextBinaryNumber(int *bits, int numBits) {
	
	for (int i=numBits-1; i >= 0; i--) {
		if (bits[i] == 0) {
			bits[i] = 1;
			return;
		}
		bits[i] = 0;
	}
}


int findSingleSolution(
	int *equ, int numBools, int numClauses, int *sol, 
	int *candidate, int numCandidates) {
	
	// prepare candidate solutions
	for (int i=0; i < numBools; i++)
		candidate[i] = 0;
	
	// whether an existing solution has been found
	int foundSol = 0;
	
	// check every possible bool sequence
	for (int i=0; i < numCandidates; i++) {
		if (isValidSolution(equ, candidate, numClauses)) {
			
			// if there's an existing solution, there's no "single" solution
			if (foundSol)
				return 0;
			
			// copy candidate to sol
			foundSol = 1;
			for (int j=0; j < numBools; j++)
				sol[j] = candidate[j];
			
		}
		nextBinaryNumber(candidate, numBools);
	}
	
	// return whether a single sol was found
	return foundSol;
}


int getRandomInteger(int min, int max) {
	return min + rand()/(RAND_MAX/(max - min + 1) + 1);
}


int arrayContains(int *array, int length, int element) {
	for (int i=0; i < length; i++) {
		if (array[i] == element)
			return 1;
	}
	return 0;
}


/*
int compareIntegers(const void *a,const void *b) {
	int *x = (int *) a;
	int *y = (int *) b;
	return *x - *y;
}
*/

void swapArrayInts(int *arr, int i, int j) {
	int buffer = arr[i];
	arr[i] = arr[j];
	arr[j] = buffer;
}


void sortTripleIntArray(int *triple) {
	if (triple[1] < triple[0])
		swapArrayInts(triple, 0, 1);
	if (triple[2] < triple[1])
		swapArrayInts(triple, 1, 2);
	if (triple[1] < triple[0])
		swapArrayInts(triple, 0, 1);
}


void getRandomClause(int numBools, int *clause) {
	
	// randomly pick CLAUSE_SIZE distinct indices in [1, numBools]
	for (int i=0; i < CLAUSE_SIZE; i++) {
		do clause[i] = getRandomInteger(1, numBools);
		while (arrayContains(clause, i, clause[i]));
	}
	
	// sort clause indices for easy comparing
	/*
	qsort(clause, CLAUSE_SIZE, sizeof(int), compareIntegers);
	 */
	sortTripleIntArray(clause);
	
	// negate some clauses
	for (int i=0; i < CLAUSE_SIZE; i++) {
		if (getRandomInteger(0, 1))
			clause[i] *= -1;
	}
}


int allTrue(int *array, int length) {
	for (int i=0; i < length; i++) 
		if (!array[i])
			return 0;
	return 1;
}


int equContainsClause(int *equ, int numClauses, int* clause) {
	
	// iterate number of clauses already in equ
	for (int clauseNum=0; clauseNum < numClauses; clauseNum++) {
		
		// starting index of the current clause in equ
		int startInd = clauseNum * CLAUSE_SIZE;
		
		// compare current clause to passed clause
		int clauseFound = 1;
		for (int i=0; i < CLAUSE_SIZE; i++) {
			int equInd = startInd + i;
			if (equ[equInd] != clause[i]) {
				clauseFound = 0;
				break;
			}
		}
		
		// the passed clause didn't disagree with the current clause
		if (clauseFound)
			return 1;
	}
	
	// passed clause was not found in equ
	return 0;
}


void getRandomEquAndSol(int numBools, int **equ, int **sol, int *numClauses) {
	
	*numClauses = getNumClauses(numBools);
	*equ = malloc(*numClauses * CLAUSE_SIZE * sizeof(int));
	*sol = malloc(numBools * sizeof(int));
	
	// prepare state-checking structures
	int numClausesMade = 0;	// the number of clauses constructed so far in equ
	int foundSingleSol = 0;	// whether a single solution was found for the equ
	int *indPresent;			// whether bool[ind] is present anywhere in the equ
	int *indNotPresent;		// whether NOT(bool[ind]) is present anywhere in the equ
	int clause[CLAUSE_SIZE];	// a randomly chosen clause, not necessarily satisfactory
	int *candidate;				// memory for candidate solution iteration
	
	indPresent = malloc(numBools * sizeof(int));
	indNotPresent = malloc(numBools * sizeof(int));
	candidate = malloc(numBools * sizeof(int));
	for (int i=0; i < numBools; i++)
		indPresent[i] = indNotPresent[i] = candidate[i] = 0;
	
	// solution space to explore for each considered equ
	int numCandidates = pow(2, numBools);
	
	// find equ with clauses subject to constraints
	while (
		!foundSingleSol ||
		!allTrue(indPresent, numBools) ||
		!allTrue(indNotPresent, numBools)
	) {
		
		// clear progress
		foundSingleSol = 0;
		numClausesMade = 0;
		for (int i=0; i < numBools; i++)
			indPresent[i] = indNotPresent[i] = 0;
		
		// collect unique clauses
		while (numClausesMade < *numClauses) {
						
			// generate a random (sorted) clause (indices >= 1)
			getRandomClause(numBools, clause);
			
			// discard clauses already present in equ
			if (equContainsClause(*equ, numClausesMade, clause))
				continue;
				
			// record that bools-in-clause have entered equ
			for (int i=0; i < CLAUSE_SIZE; i++) {
				if (clause[i] > 0)
					indPresent[clause[i]-1] = 1; // (indices >= 1)
				else
					indNotPresent[-clause[i]-1] = 1;
			}
				
			// starting index of new clause in equ
			int startInd = numClausesMade * CLAUSE_SIZE;
			numClausesMade++;
			
			// write new clause to equ
			for (int i=0; i < CLAUSE_SIZE; i++)
				(*equ)[startInd + i] = clause[i];
		}
		
		// check if this equ is solvable (with at most 1 sol)
		foundSingleSol = findSingleSolution(
			*equ, numBools, *numClauses, *sol, 
			candidate, numCandidates);
	}
	
	// free memory
	free(indPresent);
	free(indNotPresent);
	free(candidate);
	
	// equ and sol must be freed!
}


int getNumClauses(int numBools) {
	return round(CLAUSE_BOOL_RATIO * numBools);
}


int mainPrivate(int narg, char* varg[]) {
	
	
	/*
	 * FETCH PARAMETERS
	 */
	
	if (narg < 3 || narg > 6) {
		printf("ERROR: call as ./3SATGenerator seed "
			   "numBools[, printEqu][, equFN][, solFN]\n");
		return 1;
	}
	
	// 3SATGenerator rSeed, numBools[, printFlag][, equFN][, solFN]
	int rSeed = atoi(varg[1]);
	int numBools = atoi(varg[2]);
	int printFlag = (narg >= 3+1)? atoi(varg[3]):1;
	int equFNPassed = (narg >= 4+1);
	int solFNPassed = (narg == 5+1);
	
	
	/*
	 * PREPARE 3SAT PROBLEM
	 */
	 
	 // ensure randomness
	 srand(rSeed);
	 
	// generate random problem
	int numClauses;
	int *equ, *sol;
	getRandomEquAndSol(numBools, &equ, &sol, &numClauses);
	
	// output equation and solution
	if (equFNPassed)
		saveEquation(varg[4], equ, numClauses);
	if (printFlag)
		printEquation(equ, numClauses);
	if (solFNPassed)
		saveSolution(varg[5], sol, numBools);
	
	
	/*
	 * FREE MEMORY
	 */
	
	free(equ);
	free(sol);
	return 0;
}