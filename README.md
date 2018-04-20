Variational Imag-Time 3SAT Solver
========

# Setup

You'll need to install the GNU Scientific library; instructions [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903).

The executable `Variational3SATSolver` has main in `variational_3SAT_solver.c` and includes
- QuEST
- `sat_generator` which generates the 3SAT problems **:^)**
- `hamiltonian_builder` which builds the Hamiltonians **:^)**
- `param_evolver` which evolves the parameters **:^)** with the variational imagianry time method
- `mmaformatter` which formats the output for importing into Mathematica

The QuEST `makefile` must be adapted for GSL. The changes are outlined [here](https://gist.github.com/TysonRayJones/af7bedcdb8dc59868c7966232b4da903) and my makefile is included in this repo.

# Running

`Variational3SATSolver` accepts 9 command-line arguments:
- `num_bools`: 3SAT problem size (a random instance will be generated)
 - `num_params`: the number of parameterised gates in the ansatz circuit
 -  `rseed`: random seed for 3SAT generation 
 - `threshold`: the probability of the solution state which when achieved, parameter evolution ends
 - `timestep`: multiplies the change the parameters experience at each imag-time evolution step. Passing `0` will use the smallest stable time-step, ~**max(energyEigVal)**
 - `max_iters`: maximum iterations to perform before halting (if not already converged to `threshold`)
 -  `wrap_params`: whether (0 or 1) parameters will be wrapped-around to remain in **[0, 2 PI)** after each evolution step.
 - `deriv_accuracy`: order of finite difference deriv approx to use (1 to 4) in constructing A and C
 - `matrix_noise`: percent (0 to 1) noise to vary A and C elements by before solving for parameter evolution


# Customising

The parameters are updated by the method
```C
evolveOutcome evolveParams(
	evolverMemory *mem, 
	void (*ansatzCircuit)(MultiQubit, double*, int), 
	int (*illPosedRecoveryMethod)(evolverMemory*),
	MultiQubit qubits, 
	double* params, 
	double* diagHamiltonian, 
	double timeStepSize, 
	int wrapParams
);
```
where `ansatzCircuit` and `illPosedRecoveryMethod` are functions you can define and pass, and `mem` is a data structure you must first create with
```C
prepareEvolverMemory(MultiQubit qubits, int numParams);
```
-------
 `ansatzCircuit` must be a function which accepts a `MultiQubit` instance, the parameter list and its length, and applies an ansatz circuit for those parameter values to `MultiQubit`. It's important this ansatz circuit sets the initial state of the qubits. For example
```C
void myAnsatz(MultiQubit qubits, double* params, int numParams) {
	
	initStateZero(&qubits);

	for (int p=0; p < numParams; p++)
		rotateX(qubits, p % qubits.numQubits, params[p]);
}
```
could be passed to ```evolveOutcome(..., myAnsatz, ...) ```
In lieu of your own code, `defaultAnsatzCircuit` can be passed.

-------

`illPosedRecoveryMethod` is called when LU decomposition fails to solve for the change in parameters in the given step of the algorithm, and allows an alternative approximation to be used.
The function must modify `mem->paramChange` (which is a GSL vector) using the `mem->matrA` and `mem->vecC` matrix/vector, to generate a solution to `matrA paramChange = vecC`, and should return `1` if it also numerically fails.

In lieu of your own code, you can pass `approxParamsByLeastSquares`, `approxParamsByRemovingVar` or `approxParamsByTSVD`; see the code-doc in `param_evolver.h` for details.

Underdetermined equations of `paramChange` can be realised by changing `initStateZero` to `initStatePlus` in `defaultAnsatzCircuit`.

# Output

`Variational3SATSolver` (through `variational_3SAT_solver.c`) outputs a Mathematica association (using [this](https://gist.github.com/TysonRayJones/f25cb847aadd70aef3f5e0f8fae04947)) to file, which contains
- expected energy over iterations
- probability of the solution state over iterations
- evolution of every parameter
- which iterations underwent recovery numerical methods

and is simply read via 
```Mathematica 
SetDirectory @ NotebookDirectory[];
data = Get["wickSATdata.txt"];
```

An example `analysis.nb` notebook is included which plots the above data:
![mm](https://qtechtheory.org/wp-content/uploads/2018/04/example.png)
