#ifndef TRUE_EVOLVER_H_
#define TRUE_EVOLVER_H_

#include <QuEST.h>
#include "hamiltonian_builder.h"

void evolveWavefunction(MultiQubit qubits, Hamiltonian hamil, double timeStep, complex double* hamilState);

#endif // TRUE_EVOLVER_H_