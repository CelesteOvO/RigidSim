#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"

Contact::Contact() : p(), n(), t1(), t2(), mu(0.4f),
    body0(nullptr), body1(nullptr), k(1e6f), b(1e5f), index(-1), pene(0.0f)
{

}

Contact::Contact(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _p, const Eigen::Vector3f& _n, float _pene) :
    p(_p), n(_n), t1(), t2(), mu(0.4f), pene(_pene), body0(_body0), body1(_body1), k(1e6f), b(1e5f), index(-1)
{
    J0.setZero(3, 6);
    J1.setZero(3, 6);
    J0Minv.setZero(3, 6);
    J1Minv.setZero(3, 6);
    lambda.setZero(3);
    phi.setZero(3);
    phi(0) = _pene;

    body0->contacts.push_back(this);
    body1->contacts.push_back(this);
}

Contact::~Contact()
{

}

void Contact::computeContactFrame()
{
    Eigen::Vector3f dir = Eigen::Vector3f(1,0,0);
    t1 = dir - (dir.dot(n))*n;
    if( t1.norm() < 1e-5f )
    {
        // Fail-safe: use axis-aligned direction (0,0,-1)
        t1 = -n.cross( Eigen::Vector3f(0,0,-1) );
    }
    t1.normalize();

    t2 = n.cross(t1);
    t2.normalize();

    /*Eigen::Vector3f tangent1 = Eigen::Vector3f::Zero();
    Eigen::Vector3f tangent2 = Eigen::Vector3f::Zero();
    Eigen::Vector3f n = this->n;

    if (n.x() == 0 && n.y() == 0) {
        tangent1 = Eigen::Vector3f(1, 0, 0);
    } else {
        tangent1 = Eigen::Vector3f(n.y(), -n.x(), 0);
    }
    tangent1.normalize();
    tangent2 = n.cross(tangent1);
    tangent2.normalize();

    this->t1 = tangent1;
    this->t2 = tangent2;*/
}

void Contact::computeJacobian()
{
    const Eigen::Vector3f r0 = p - body0->x;
    const Eigen::Vector3f r1 = p - body1->x;

    J0.setZero(3, 6);
    J1.setZero(3, 6);
    J0Minv.setZero(3, 6);
    J1Minv.setZero(3, 6);
    lambda.setZero(3);
    phi.setZero(3);
    phi(0) = pene;

    // normal row (non-interpenetration)
    J0.block(0, 0, 1, 3) = n.transpose();
    J0.block(0, 3, 1, 3) = r0.cross(n).transpose();
    J1.block(0, 0, 1, 3) = -n.transpose();
    J1.block(0, 3, 1, 3) = -r1.cross(n).transpose();
    // tangent 1 (friction)
    J0.block(1, 0, 1, 3) = t1.transpose();
    J0.block(1, 3, 1, 3) = r0.cross(t1).transpose();
    J1.block(1, 0, 1, 3) = -t1.transpose();
    J1.block(1, 3, 1, 3) = -r1.cross(t1).transpose();
    // tangent 2 (friction)
    J0.block(2, 0, 1, 3) = t2.transpose();
    J0.block(2, 3, 1, 3) = r0.cross(t2).transpose();
    J1.block(2, 0, 1, 3) = -t2.transpose();
    J1.block(2, 3, 1, 3) = -r1.cross(t2).transpose();

    J0Minv.block(0,0,3,3) = (1.0f/body0->mass) * J0.block(0, 0, 3, 3);
    J0Minv.block(0,3,3,3) = J0.block(0, 3, 3, 3) * body0->Iinv;
    J1Minv.block(0,0,3,3) = (1.0f/body1->mass) * J1.block(0, 0, 3, 3);
    J1Minv.block(0,3,3,3) = J1.block(0, 3, 3, 3) * body1->Iinv;
}
