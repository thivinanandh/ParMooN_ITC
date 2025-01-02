#include "IntelPardisoSolver.h"
#include <cstdlib>
#include <omp.h>
#include <mkl.h>
#include <iostream>

IntelPardisoSolver::IntelPardisoSolver() 
    : is_initialized(false)
{
    // Initialize solver parameters
    mtype = 11;    // Real unsymmetric matrix
    maxfct = 1;    // Maximum number of numerical factorizations
    mnum = 1;      // Which factorization to use
    msglvl = 0;    // Print statistical information
    error = 0;     // Initialize error flag

    // Initialize iparm
    for(int i = 0; i < 64; i++) 
        iparm[i] = 0;

    // Set up parameters for PARDISO
    iparm[0] = 1;    // No solver default
    iparm[1] = 2;    // Fill-in reordering from METIS
    iparm[3] = 0;    // No iterative-direct algorithm
    iparm[4] = 0;    // No user fill-in reducing permutation
    iparm[5] = 0;    // Write solution into x
    iparm[7] = 2;    // Max numbers of iterative refinement steps
    iparm[9] = 13;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;   // Use nonsymmetric permutation and scaling MPS
    iparm[12] = 1;   // Maximum weighted matching algorithm is switched-on
    iparm[34] = 1;   // Zero-based indexing

    // Initialize the internal solver memory pointer
    for(int i = 0; i < 64; i++) 
        pt[i] = 0;
}

IntelPardisoSolver::~IntelPardisoSolver() 
{
    cleanup();
}

void IntelPardisoSolver::initialize(MKL_INT N_DOF, int* rowptr, int* colIndex, double* entries) 
{
    n = N_DOF;
    MKL_INT nrhs = 1;
    double ddum = 0.0;    // Dummy double variable
    MKL_INT idum = 0;     // Dummy integer variable

    // Symbolic factorization
    phase = 11;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, entries, rowptr, colIndex, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    
    if(error != 0) {
        printf("\nERROR during symbolic factorization: %d\n", error);
        exit(1);
    }

    printf("\nReordering completed ...\n");
    printf("Number of nonzeros in factors = %d\n", iparm[17]);
    printf("Number of factorization MFLOPS = %d\n", iparm[18]);
    
    is_initialized = true;
}

void IntelPardisoSolver::solve(double* entries, int* rowptr, int* colIndex, 
                              double* rhs, double* solution) 
{
    if(!is_initialized) {
        printf("\nERROR: Solver not initialized! Call initialize() first.\n");
        return;
    }

    MKL_INT nrhs = 1;
    double ddum = 0.0;    // Dummy double variable
    MKL_INT idum = 0;     // Dummy integer variable

    // Numerical factorization
    phase = 22;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, entries, rowptr, colIndex, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    
    if(error != 0) {
        printf("\nERROR during numerical factorization: %d\n", error);
        exit(2);
    }

    // Back substitution and iterative refinement
    phase = 33;
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, entries, rowptr, colIndex, &idum, &nrhs,
            iparm, &msglvl, rhs, solution, &error);
    
    if(error != 0) {
        printf("\nERROR during solution: %d\n", error);
        exit(3);
    }
}

void IntelPardisoSolver::cleanup() 
{
    if(!is_initialized) 
        return;

    MKL_INT nrhs = 1;
    double ddum = 0.0;    // Dummy double variable
    MKL_INT idum = 0;     // Dummy integer variable

    phase = -1;  // Release internal memory
    PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &n, &ddum, NULL, NULL, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    
    is_initialized = false;
}