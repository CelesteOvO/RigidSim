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
            // Test for box-box collision
            else if( body0->geometry->getType() == kBox &&
                     body1->geometry->getType() == kBox )
            {
                collisionDetectBoxBox(body0, body1);
            }
            // Test for box-box collision (order swap)
            else if( body1->geometry->getType() == kBox &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectBoxBox(body1, body0);
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

///TODO: 有问题
void CollisionDetect::collisionDetectBoxBox(RigidBody *body0, RigidBody *body1) {
    std::vector<Eigen::Vector3f> axis;

    for (int i = 0; i < 3; i++) {
        axis.emplace_back(body0->q.toRotationMatrix().col(i));
        axis.emplace_back(body1->q.toRotationMatrix().matrix().col(i));
    }

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            Eigen::Vector3f temp = axis[i].cross(axis[3+j]);
            if(temp.norm() > 1e-6f)
                axis.emplace_back(temp.normalized());
        }
    }

    float minOverlap = 1e10f;
    unsigned int bestAxis = 15;
    for (unsigned int i = 0; i < 15; i++)
    {
        Eigen::Vector3f temp = axis[i];
        if(temp.norm() < 1e-6f)
            continue;
        temp.normalize();
        float overlap = penetrationOnAxis(body0, body1, temp);
        if(overlap < 0)
            return;
        if(overlap < minOverlap)
        {
            minOverlap = overlap;
            bestAxis = i;
        }
    }

    if(bestAxis <= 2)
    {
        collisionDetectFaceVertex(body0, body1, axis[bestAxis], minOverlap);
    }else if(bestAxis <= 5)
    {
        collisionDetectFaceVertex(body1, body0, axis[bestAxis], minOverlap);
    }else{
        int oneAxisIndex = (bestAxis - 6) / 3;
        int twoAxisIndex = bestAxis % 3;
        collisionDetectEdgeEdge(body0, body1, axis[bestAxis], minOverlap, oneAxisIndex, twoAxisIndex);
    }
}

float CollisionDetect::penetrationOnAxis(RigidBody *pBody, RigidBody *pBody1, const Eigen::Vector3f& axis) {
    float one = transformToAxis(pBody, axis);
    float two = transformToAxis(pBody1, axis);
    float center = abs((pBody->x - pBody1->x).dot(axis));
    return one + two - center;
}

float CollisionDetect::transformToAxis(RigidBody *pBody, const Eigen::Vector3f &axis) {
    Box* box0 = dynamic_cast<Box*>(pBody->geometry.get());
    float projectionMax = -1e10f;
    float projectionMin = 1e10f;
    for (float cx = -1.0f * box0->dim(0) / 2;cx < box0->dim(0);cx += box0->dim(0))
    {
        for (float cy = -1.0f * box0->dim(1) / 2;cy < box0->dim(1);cy += box0->dim(1))
        {
            for (float cz = -1.0f * box0->dim(2) / 2;cz < box0->dim(2);cz += box0->dim(2))
            {
                Eigen::Vector3f vertex = Eigen::Vector3f(cx,cy,cz);
                float projection = vertex.dot(axis);
                projectionMax = std::max(projectionMax,projection);
                projectionMin = std::min(projectionMin,projection);
            }
        }
    }
    return (projectionMax - projectionMin) / 2.0f;
}

void CollisionDetect::collisionDetectFaceVertex(RigidBody *body0, RigidBody *body1, Eigen::Vector3f axis,
                                                float penetration) {
    Eigen::Vector3f toCenter = body1->x - body0->x;
    if(axis.dot(toCenter) > 0)
        axis *= -1;
    Box* box1 = dynamic_cast<Box*>(body1->geometry.get());
    Eigen::Vector3f vertex = body1->x + box1->dim(0) * body1->q.toRotationMatrix().col(0) + box1->dim(1) * body1->q.toRotationMatrix().col(1) + box1->dim(2) * body1->q.toRotationMatrix().col(2);
    if(body1->q.toRotationMatrix().col(0).dot(axis) < 0)
        vertex -= body1->q.toRotationMatrix().col(0);
    if(body1->q.toRotationMatrix().col(1).dot(axis) < 0)
        vertex -= body1->q.toRotationMatrix().col(1);
    if(body1->q.toRotationMatrix().col(2).dot(axis) < 0)
        vertex -= body1->q.toRotationMatrix().col(2);

    m_contacts.push_back(new Contact(body0, body1, vertex, axis, penetration));
}

void
CollisionDetect::collisionDetectEdgeEdge(RigidBody *body0, RigidBody *body1, Eigen::Vector3f axis, float penetration, int oneAxisIndex, int twoAxisIndex) {
    Eigen::Vector3f toCenter = body1->x - body0->x;
    if(axis.dot(toCenter) > 0)
        axis *= -1;
    Eigen::Vector3f ptOnEdgeOne = body0->x;
    Eigen::Vector3f ptOnEdgeTwo = body1->x;
    Box* box0 = dynamic_cast<Box*>(body0->geometry.get());
    Box* box1 = dynamic_cast<Box*>(body1->geometry.get());
    for(int i = 0; i < 3; i++)
    {
        if(i == oneAxisIndex){}
        else if(body0->q.toRotationMatrix().col(i).dot(axis) > 0)
            ptOnEdgeOne -= box0->dim(i) / 2.0f * body0->q.toRotationMatrix().col(i);
        else
            ptOnEdgeOne += box0->dim(i) / 2.0f * body0->q.toRotationMatrix().col(i);
        if(i == twoAxisIndex){}
        else if(body1->q.toRotationMatrix().col(i).dot(axis) < 0)
            ptOnEdgeTwo -= box1->dim(i) / 2.0f * body1->q.toRotationMatrix().col(i);
        else
            ptOnEdgeTwo += box1->dim(i) / 2.0f * body1->q.toRotationMatrix().col(i);
    }
    Eigen::Vector3f axisOne = body0->q.toRotationMatrix().col(oneAxisIndex);
    Eigen::Vector3f axisTwo = body1->q.toRotationMatrix().col(twoAxisIndex);
    Eigen::Vector3f toSt = ptOnEdgeOne - ptOnEdgeTwo;

    float dpStaOne = axisOne.dot(toSt);
    float dpStaTwo = axisTwo.dot(toSt);

    float smOne = axisOne.squaredNorm();
    float smTwo = axisTwo.squaredNorm();

    double dotProductEdges = axisTwo.dot(axisOne);
    double denom = smOne * smTwo - dotProductEdges * dotProductEdges;
    double ta = (dotProductEdges * dpStaTwo - smTwo * dpStaOne) / denom;
    double tb = (smOne * dpStaTwo - dotProductEdges * dpStaOne) / denom;

    Eigen::Vector3f ptOne = ptOnEdgeOne + ta * axisOne;
    Eigen::Vector3f ptTwo = ptOnEdgeTwo + tb * axisTwo;

    m_contacts.push_back(new Contact(body0, body1, 0.5f * (ptOne + ptTwo), axis, penetration));
}