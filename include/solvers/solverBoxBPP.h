#ifndef RIGIDBODYTUTORIAL_SOLVERBOXBPP_H
#define RIGIDBODYTUTORIAL_SOLVERBOXBPP_H
#pragma once

#include "solvers/Solver.h"

// Block Principal Pivoting (BPP) Boxed LCP solver.
//
class SolverBoxBPP : public Solver
{
public:
    explicit SolverBoxBPP(RigidBodySystem* _rigidBodySystem);

    // Implement Block Principal Pivoting method.
    //
    virtual void solve(float h);

};
#endif //RIGIDBODYTUTORIAL_SOLVERBOXBPP_H
