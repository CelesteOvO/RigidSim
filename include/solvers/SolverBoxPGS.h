#pragma once

#include "solvers/Solver.h"

// Projected Gauss-Seidel (PGS) Boxed LCP solver.
// This implementation uses a matrix-free approach, where only
// the non-zero blocks of the lead matrix are assembled and stored.
//
class SolverBoxPGS : public Solver
{
public:

    explicit SolverBoxPGS(RigidBodySystem* _rigidBodySystem);

    // Implement PGS method that solves for the constraint forces in @a m_rigidBodySystem.
    //
    void solve(float h) override;

};
