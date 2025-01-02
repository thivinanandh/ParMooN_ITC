#ifndef __INTELPARDISOSOLVER__
#define __INTELPARDISOSOLVER__

#include <mkl.h>
// #include <cstdio>

class IntelPardisoSolver 
{
private:
    void *pt[64];          // Internal solver memory pointer
    MKL_INT iparm[64];     // Pardiso control parameters
    MKL_INT maxfct;        // Maximum number of numerical factorizations
    MKL_INT mnum;          // Which factorization to use
    MKL_INT mtype;         // Matrix type
    MKL_INT phase;         // Current solution phase
    MKL_INT error;         // Error flag
    MKL_INT msglvl;        // Message level
    MKL_INT n;             // Matrix dimension
    bool is_initialized;    // Initialization flag

public:
    // Constructor
    IntelPardisoSolver();

    // Destructor
    ~IntelPardisoSolver();

    // Initialize solver and perform symbolic factorization
    void initialize(MKL_INT N_DOF, int* rowptr, int* colIndex, double* entries);

    // Solve the system (numerical factorization + solution)
    void solve(double* entries, int* rowptr, int* colIndex, 
               double* rhs, double* solution);

    // Clean up when completely done
    void cleanup();

    // Utility functions
    bool isInitialized() const { return is_initialized; }
    // MKL_INT getError() const { return error; }
    
    // Disable copy constructor and assignment operator
    IntelPardisoSolver(const IntelPardisoSolver&) = delete;
    IntelPardisoSolver& operator=(const IntelPardisoSolver&) = delete;
};

#endif // __INTELPARDISOSOLVER__