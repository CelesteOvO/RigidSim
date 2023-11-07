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
    // 求出两个box的八个顶点
    std::vector<Eigen::Vector3f> corners1, corners2;
    corners1 = getCorners(body0);
    corners2 = getCorners(body1);

    std::vector<Eigen::Vector3f> faceNormals1, faceNormals2;
    faceNormals1 = getFaceNormal(body0, corners1);
    faceNormals2 = getFaceNormal(body1, corners2);

    std::vector<Eigen::Vector3f> edgeNormals;
    edgeNormals = getEdgeNormal(body0, body1, corners1, corners2);

    // 把所有的normal放到一个vector中
    std::vector<Eigen::Vector3f> normals;
    normals.insert(normals.end(), faceNormals1.begin(), faceNormals1.end());
    normals.insert(normals.end(), faceNormals2.begin(), faceNormals2.end());
    normals.insert(normals.end(), edgeNormals.begin(), edgeNormals.end());

    Eigen::Vector3f d = body0->x - body1->x;
    // 测试所有的normal
    float minOverlap = 1000000.0f;
    int axis = testSAT(minOverlap, normals, body0, body1, corners1, corners2);

    if (axis == -1)
        return;
    else if(axis <= 2) /// 发生face-vertex collision
    {
        collisionDetectFaceVertex(body0, body1, normals[axis], minOverlap);
    }
    else if(axis <= 5)
    {
        collisionDetectFaceVertex(body1, body0, normals[axis], minOverlap);
    }
    else
    {
        int oneAxisIndex = (axis - 6) / 3;
        int twoAxisIndex = (axis - 6) % 3;
        collisionDetectEdgeEdge(body0, body1, normals[axis], minOverlap, oneAxisIndex, twoAxisIndex);
    }
}

std::vector<Eigen::Vector3f> CollisionDetect::getCorners(RigidBody *body) {
    std::vector<Eigen::Vector3f> corners;

    Box* box = dynamic_cast<Box*>(body->geometry.get());

    Eigen::Vector3f center = body->x;
    Eigen::Matrix3f orientation = body->q.toRotationMatrix();

    Eigen::Vector3f halfExtents = box->dim / 2.0f;

    // 盒子的局部坐标系的轴向量
    Eigen::Vector3f localXAxis = orientation.col(0);
    Eigen::Vector3f localYAxis = orientation.col(1);
    Eigen::Vector3f localZAxis = orientation.col(2);

    // 通过半边长和轴向量计算顶点位置
    corners.emplace_back(center + halfExtents[0] * localXAxis + halfExtents[1] * localYAxis + halfExtents[2] * localZAxis);
    corners.emplace_back(center - halfExtents[0] * localXAxis + halfExtents[1] * localYAxis + halfExtents[2] * localZAxis);
    corners.emplace_back(center + halfExtents[0] * localXAxis - halfExtents[1] * localYAxis + halfExtents[2] * localZAxis);
    corners.emplace_back(center - halfExtents[0] * localXAxis - halfExtents[1] * localYAxis + halfExtents[2] * localZAxis);
    corners.emplace_back(center + halfExtents[0] * localXAxis + halfExtents[1] * localYAxis - halfExtents[2] * localZAxis);
    corners.emplace_back(center - halfExtents[0] * localXAxis + halfExtents[1] * localYAxis - halfExtents[2] * localZAxis);
    corners.emplace_back(center + halfExtents[0] * localXAxis - halfExtents[1] * localYAxis - halfExtents[2] * localZAxis);
    corners.emplace_back(center - halfExtents[0] * localXAxis - halfExtents[1] * localYAxis - halfExtents[2] * localZAxis);

    return corners;
}

std::vector<Eigen::Vector3f> CollisionDetect::getFaceNormal(RigidBody *body, std::vector<Eigen::Vector3f> corners) {
    std::vector<Eigen::Vector3f> normals;
    Eigen::Vector3f U = corners[1] - corners[0];
    Eigen::Vector3f V = corners[2] - corners[0];
    Eigen::Vector3f N = U.cross(V).normalized();
    normals.push_back(N);

    U = corners[6] - corners[2];
    V = corners[3] - corners[2];
    N = U.cross(V).normalized();
    normals.push_back(N);

    U = corners[4] - corners[6];
    V = corners[2] - corners[6];
    N = U.cross(V).normalized();
    normals.push_back(N);

    return normals;
}

std::vector<Eigen::Vector3f>
CollisionDetect::getEdgeNormal(RigidBody *body1, RigidBody *body2, std::vector<Eigen::Vector3f> corners1,
                               std::vector<Eigen::Vector3f> corners2) {
    std::vector<Eigen::Vector3f> edgeNormals;
    Eigen::Vector3f body1Edge[3];
    Eigen::Vector3f body2Edge[3];

    body1Edge[0] = corners1[1] - corners1[0];
    body1Edge[1] = corners1[2] - corners1[0];
    body1Edge[2] = corners1[4] - corners1[0];

    body2Edge[0] = corners2[1] - corners2[0];
    body2Edge[1] = corners2[2] - corners2[0];
    body2Edge[2] = corners2[4] - corners2[0];

    for (auto & i : body1Edge) {
        for (const auto & j : body2Edge) {
            Eigen::Vector3f N = i.cross(j).normalized();
            edgeNormals.push_back(N);
        }
    }
    return edgeNormals;
}

int CollisionDetect::testSAT(float minOverlap, std::vector<Eigen::Vector3f> normals, RigidBody* body0, RigidBody* body1, const std::vector<Eigen::Vector3f>& corners1, const std::vector<Eigen::Vector3f>& corners2) {
    int axisNum = 0;

    for (int i = 0; i < 15; i++)
    {
        Eigen::Vector3f axis = normals[i];
        if(axis.norm() < 1e-6f)
            continue;
        axis.normalize();

        float hc1 = transformToAxis(body0, axis, corners1);
        float hc2 = transformToAxis(body1, axis, corners2);
        float center = abs((body0->x - body1->x).dot(axis));
        float overlap = hc1 + hc2 - center;

        if(overlap < 0)
            return -1;
        if(overlap < minOverlap)
        {
            minOverlap = overlap;
            axisNum = i;
        }
    }
    return axisNum;
}

float CollisionDetect::transformToAxis(RigidBody *pBody, const Eigen::Vector3f &axis, const std::vector<Eigen::Vector3f>& corners) {
    float projectionMax = -1e10f;
    float projectionMin = 1e10f;
    for (auto & corner : corners) {
        float projection = corner.dot(axis);
        if (projection > projectionMax)
            projectionMax = projection;
        if (projection < projectionMin)
            projectionMin = projection;
    }
    float hc = 0.5f * (projectionMax - projectionMin);
    return hc;
}

void CollisionDetect::collisionDetectFaceVertex(RigidBody *body0, RigidBody *body1, Eigen::Vector3f axis,
                                                float penetration) {
    Eigen::Vector3f toCenter = body0->x - body1->x;
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

void CollisionDetect::collisionDetectEdgeEdge(RigidBody *body0, RigidBody *body1, Eigen::Vector3f axis, float penetration, int oneAxisIndex, int twoAxisIndex) {
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






