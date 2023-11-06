//
// Created by LiYifan on 2023/11/6.
//
#include "solverBoxBPP.h"

#include "rigidbody/RigidBodySystem.h"
#include "rigidbody/RigidBody.h"
#include "contact/Contact.h"

#include <Eigen/Dense>

SolverBoxBPP::SolverBoxBPP(RigidBodySystem *_rigidBodySystem) : Solver(_rigidBodySystem) {

}

void SolverBoxBPP::solve(float h) {
    auto contacts = m_rigidBodySystem->getContacts();
    const unsigned int numContacts = contacts.size();
    if( numContacts > 0 )
    {

    }
}

