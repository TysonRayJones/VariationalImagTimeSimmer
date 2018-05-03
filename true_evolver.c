
#include "true_evolver.h"

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <QuEST.h>

#include "hamiltonian_builder.h"


/**
 * Evolves the wavefunction under the schrodinger equation by amount timeStep 
 * in imaginary time, using fourth-order runge-kutta
 */
void evolveUnderDiagHamil(MultiQubit qubits, double* hamil, double timeStep) {
	
	// DE: dy/ds = - H y
	
	// get fourth-order runge-kutta imag-time-Schrodinger-equ terms
	double complex k1[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k1[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i));
	
	double complex k2[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k2[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k1[i]/2.0);
		
	double complex k3[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k3[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k2[i]/2.0);
	
	double complex k4[qubits.numAmps];
	for (long long int i=0LL; i < qubits.numAmps; i++) 
		k4[i] = - timeStep * hamil[i] * (getRealAmpEl(qubits, i) + I*getImagAmpEl(qubits, i) + k3[i]);
	
	// update wavefunction
	double norm = 0;
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		double complex dPsi = (k1[i] + k4[i])/6.0 + (k2[i] + k3[i])/3.0;
		qubits.stateVec.real[i] += creal(dPsi);
		qubits.stateVec.imag[i] += cimag(dPsi);
		
		norm += pow(qubits.stateVec.real[i],2) + pow(qubits.stateVec.imag[i], 2);
	}
	
	// normalise wavefunction
	norm = sqrt(norm);
	for (long long int i=0LL; i < qubits.numAmps; i++) {
		qubits.stateVec.real[i] /= norm;
		qubits.stateVec.imag[i] /= norm;
	}
}


void evolveUnderPauliHamil(MultiQubit qubits, Hamiltonian pauliHamil, double timeStep, complex double* hamilState) {
	
	// load H |psi> into hamilState
	applyHamil(hamilState, qubits, pauliHamil);
	
	// TODO
}


void evolveWavefunction(MultiQubit qubits, Hamiltonian hamil, double timeStep, complex double* hamilState) {
	
	if (hamil.type == DIAGONAL)
		evolveUnderDiagHamil(qubits, hamil.diagHamil, timeStep);
		
	if (hamil.type == PAULI_TERMS) {
		
		printf("evolveWavefunction (true imaginary evolution) for Pauli Hamiltonians is not yet supported!!\n");
		evolveUnderPauliHamil(qubits, hamil, timeStep, hamilState);
		
	}
}