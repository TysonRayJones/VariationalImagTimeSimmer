#ifndef SAT_GEN_H_
#define SAT_GEN_H_


int *loadEquation(char *filename, int *numBools, int *numClauses);

void saveEquation(char *filename, int *equ, int numClauses);

void saveSolution(char *filename, int *sol, int numBools);

void printEquation(int *equ, int numClauses);

int isValidSolution(int *equ, int *sol, int numClauses);

int findSingleSolution(int *equ, int numBools, int numClauses, int *sol, int *candidates, int numCandidates);

void getRandomEquAndSol(int numBools, int **equ, int **sol, int *numClauses);

int getNumClauses(int numBools);

void nextBinaryNumber(int *bits, int numBits);

void convertToBinary(int number, int numBits, int *outputBits);

int getRandomInteger(int min, int max);


#endif // SAT_GEN_H_