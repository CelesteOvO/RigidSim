#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

namespace
{
    // Plane-point collision test.
    //
    // Inputs:
    //   p - The point to test.
    //   plane_p - The plane origin.
    //   plane_n  - Direction perpendicular to the plane (the normal).
    // Outputs:
    //   phi - The penetration depth.
    // Returns:
    //   True if the point intersects the plane.
    //
    inline bool collisionDetectPointPlane(const Eigen::Vector3f& p, const Eigen::Vector3f& plane_p, const Eigen::Vector3f& plane_n, float& phi)
    {
        const float dp = (p - plane_p).dot(plane_n);
        if (dp < 0.0f)
        {
            phi = std::min(0.0f, dp);
            return true;
        }
        return false;
    }
}

CollisionDetect::CollisionDetect(RigidBodySystem* rigidBodySystem) : m_rigidBodySystem(rigidBodySystem)
{

}

void CollisionDetect::detectCollisions()
{
    // First, clear any existing contacts.
    //
    clear();

    // Next, loop over all pairs of bodies and test for contacts.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for(unsigned int i = 0; i < bodies.size(); ++i)
    {
        for(unsigned int j = i+1; j < bodies.size(); ++j)
        {
            RigidBody* body0 = bodies[i];
            RigidBody* body1 = bodies[j];

            // Special case: skip tests for pairs of static bodies.
            //
            if (body0->fixed && body1->fixed) 
                continue;

            // Test for sphere-sphere collision.
            if( body0->geometry->getType() == kSphere &&
                body1->geometry->getType() == kSphere )
            {
                collisionDetectSphereSphere(body0, body1);
            }
            // Test for sphere-box collision
            else if( body0->geometry->getType() == kSphere &&
                     body1->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body0, body1);
            }
            // Test for box-sphere collision (order swap)
            else if( body1->geometry->getType() == kSphere &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body1, body0);
            }
            // Test for plane-box collision
            else if( body1->geometry->getType() == kPlane &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectBoxPlane(body0, body1);
            }
            // Test for plane-box collision (order swap)
            else if ( body0->geometry->getType() == kPlane &&
                      body1->geometry->getType() == kBox )
            {
                collisionDetectBoxPlane(body1, body0);
            }
            // Test for plane-sphere collision
            else if( body1->geometry->getType() == kPlane &&
                     body0->geometry->getType() == kSphere )
            {
                collisionDetectSpherePlane(body0, body1);
            }
            // Test for plane-sphere collision (order swap)
            else if ( body0->geometry->getType() == kPlane &&
                      body1->geometry->getType() == kSphere )
            {
                collisionDetectSpherePlane(body1, body0);
            }
        }
    }
}

void CollisionDetect::computeContactJacobians()
{
    // Build constraint Jacobians for all contacts
    //
    for(auto c : m_contacts)
    {
        c->computeContactFrame();
        c->computeJacobian();
    }
}

void CollisionDetect::clear()
{
    // First, remove all contacts from rigid bodies.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for (auto b : bodies)
    {
        b->contacts.clear();
    }

    // Then, cleanup the local contact array.
    //
    for(auto c : m_contacts)
    {
        delete c;
    }
    m_contacts.clear();

}

void CollisionDetect::collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1)
{
    auto* sphere0 = dynamic_cast<Sphere*>(body0->geometry.get());
    auto* sphere1 = dynamic_cast<Sphere*>(body1->geometry.get());

    // Implement sphere-sphere collision detection.
    // The function should check if a collision exists, and if it does
    // compute the contact normal, contact point, and penetration depth.
    //
    Eigen::Vector3f vec = body0->x - body1->x;

    const float rsum = (sphere0->radius + sphere1->radius);
    const float dist = vec.norm();
    if( dist < rsum )
    {
        const Eigen::Vector3f n = vec / dist;
        const Eigen::Vector3f p = 0.5f * ((body0->x - sphere0->radius*n) + (body1->x + sphere1->radius*n));
        const float phi = dist-rsum;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectSphereBox(RigidBody* body0, RigidBody* body1)
{
    auto* sphere = dynamic_cast<Sphere*>(body0->geometry.get());
    Box* box = dynamic_cast<Box*>(body1->geometry.get());

    Eigen::Vector3f clocal = body1->q.toRotationMatrix().transpose() * (body0->x - body1->x);

    Eigen::Vector3f q(0,0,0);
    for(unsigned int i = 0; i < 3; ++i)
    {
        q[i] = std::max(-box->dim[i]/2.0f, std::min(box->dim[i]/2.0f, clocal[i]));
    }

    const Eigen::Vector3f dx = clocal - q;
    const float dist = dx.norm();
    if( dist < sphere->radius )
    {
        const Eigen::Vector3f n = body1->q.toRotationMatrix() * (dx/dist);
        const Eigen::Vector3f p = body1->q.toRotationMatrix() * q + body1->x;
        const float phi = dist - sphere->radius;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectBoxPlane(RigidBody *body0, RigidBody *body1) {
    Box* box = dynamic_cast<Box*>(body0->geometry.get());
    auto* plane = dynamic_cast<Plane*>(body1->geometry.get());
    const Eigen::Vector3f pplane = body1->x;
    const Eigen::Vector3f nplane = body1->q.toRotationMatrix() * plane->normal;
    const Eigen::Vector3f plocal[8] = {
            0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1), -box->dim(2)),
            0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1),  box->dim(2)),
            0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1), -box->dim(2)),
            0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1),  box->dim(2)),
            0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1), -box->dim(2)),
            0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1),  box->dim(2)),
            0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1), -box->dim(2)),
            0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1),  box->dim(2))
    };

    for (const auto & i : plocal)
    {
        const Eigen::Vector3f pbox = body0->q.toRotationMatrix() * i + body0->x;
        float phi;
        if  ( collisionDetectPointPlane(pbox, pplane, nplane, phi) )
        {
            m_contacts.push_back( new Contact(body0, body1, pbox, nplane, phi) );
        }
    }
}

void CollisionDetect::collisionDetectSpherePlane(RigidBody *body0, RigidBody *body1) {
    auto* sphere = dynamic_cast<Sphere*>(body0->geometry.get());
    auto* plane = dynamic_cast<Plane*>(body1->geometry.get());

    Eigen::Vector3f clocal = body1->q.toRotationMatrix().transpose() * (body0->x - body1->x);

    Eigen::Vector3f q(0,0,0);
    for(unsigned int i = 0; i < 3; ++i)
    {
        q[i] = std::max(-plane->dim[i]/2.0f, std::min(plane->dim[i]/2.0f, clocal[i]));
    }

    const Eigen::Vector3f dx = clocal - q;
    const float dist = dx.norm();
    if( dist < sphere->radius )
    {
        const Eigen::Vector3f n = body1->q.toRotationMatrix() * (dx/dist);
        const Eigen::Vector3f p = body1->q.toRotationMatrix() * q + body1->x;
        const float phi = dist - sphere->radius;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}
