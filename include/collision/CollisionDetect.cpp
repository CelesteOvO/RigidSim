#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

float CollisionDetect::s_contactTime;
float CollisionDetect::s_minAxisPenetrationDepth;
Eigen::Vector3f CollisionDetect::s_contactNormal;

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
    float stepsize = 0.0167f; // m_dt / m_numSteps

    if (DynamicCheck(body0, body1, stepsize))
    {
        DeriveContacts(body0, body1);
    }
}


bool CollisionDetect::DeriveContacts(RigidBody *body0, RigidBody *body1) {

    Box* box0 = dynamic_cast<Box*>(body0->geometry.get());
    Box* box1 = dynamic_cast<Box*>(body1->geometry.get());

    s_contactTime = 0.0f;
    float tlast = std::numeric_limits<float>::max();

    Eigen::Vector3f relVelocity = (body0->xdot - body1->xdot) * 0.0167f;

    int side = intersectConfig::NONE;
    intersectConfig box0Cfg, box1Cfg;

    // 求出两个box的八个顶点
    std::vector<Eigen::Vector3f> corners1, corners2;
    corners1 = getCorners(body0);
    corners2 = getCorners(body1);

    std::vector<Eigen::Vector3f> faceNormals1, faceNormals2;
    faceNormals1 = getFaceNormal(body0, corners1);
    faceNormals2 = getFaceNormal(body1, corners2);

    for(int i = 0; i < 3; i++)
    {
        Eigen::Vector3f axis = faceNormals1[i];
        if(!FindintersectionOnAxis(body0, body1, axis, relVelocity,  s_contactTime, tlast,side, box0Cfg, box1Cfg))
            return false;
    }

    for(int i = 0; i < 3; i++)
    {
        Eigen::Vector3f axis = faceNormals2[i];
        if(!FindintersectionOnAxis(body0, body1, axis, relVelocity, s_contactTime, tlast,side, box0Cfg, box1Cfg))
            return false;
    }

    /*std::vector<Eigen::Vector3f> edgeNormals;
    edgeNormals = getEdgeNormal(body0, body1, corners1, corners2);*/

    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; ++j)
        {
            Eigen::Vector3f axis = faceNormals1[i].cross(faceNormals2[j]);
            if (abs(axis.dot(axis)) <= std::numeric_limits<float>::epsilon())// skip almost parallel edges
            {
                FindContactSet(body0, body1, side, box0Cfg, box1Cfg, s_contactTime);
                return true;
            }
            if(!FindintersectionOnAxis(body0, body1, axis, relVelocity, s_contactTime, tlast,side, box0Cfg, box1Cfg))
                return false;
        }
    }

    // velocity cross box 0 edges
    for(int i = 0; i < 3; i++)
    {
        Eigen::Vector3f axis = faceNormals1[i].cross(relVelocity);
        if(!FindintersectionOnAxis(body0, body1, axis, relVelocity, s_contactTime, tlast,side, box0Cfg, box1Cfg))
            return false;
    }

    // velocity cross box 1 edges
    for(int i = 0; i < 3; i++)
    {
        Eigen::Vector3f axis = faceNormals2[i].cross(relVelocity);
        if(!FindintersectionOnAxis(body0, body1, axis, relVelocity, s_contactTime, tlast,side, box0Cfg, box1Cfg))
            return false;
    }

    if (s_contactTime <= 0.0f || side == intersectConfig::NONE)
        return false;

    FindContactSet(body0, body1, side, box0Cfg, box1Cfg, s_contactTime);
    return true;
}

std::vector<Eigen::Vector3f> CollisionDetect::getCorners(RigidBody *body) {
    std::vector<Eigen::Vector3f> corners;

    Box* box = dynamic_cast<Box*>(body->geometry.get());

    Eigen::Vector3f center = body->x;
    Eigen::Matrix3f orientation = body->q.toRotationMatrix();

    // 计算顶点的局部坐标
    corners.emplace_back(0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1), -box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1),  box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1), -box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1),  box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1), -box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1),  box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1), -box->dim(2)));
    corners.emplace_back(0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1),  box->dim(2)));

    // 计算顶点的世界坐标
    for (auto & corner : corners) {
        corner = orientation * corner + center;
    }

    return corners;
}

std::vector<Eigen::Vector3f> CollisionDetect::getFaceNormal(RigidBody *body, std::vector<Eigen::Vector3f> corners) {
    /*std::vector<Eigen::Vector3f> normals;
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
    normals.push_back(N);*/

    std::vector<Eigen::Vector3f> normals;
    Eigen::Vector3f N = body->q.toRotationMatrix() * Eigen::Vector3f(1, 0, 0);
    normals.push_back(N.normalized());

    N = body->q.toRotationMatrix() * Eigen::Vector3f(0, 1, 0);
    normals.push_back(N.normalized());

    N = body->q.toRotationMatrix() * Eigen::Vector3f(0, 0, 1);
    normals.push_back(N.normalized());

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

bool CollisionDetect::FindintersectionOnAxis(RigidBody *obb0, RigidBody *obb1, Eigen::Vector3f &axis,
                                             Eigen::Vector3f &relVelocity, float &tfirst, float &tlast, int &side, intersectConfig &box0Cfg,
                                             intersectConfig &box1Cfg) {
    intersectConfig box0CfgStart{};
    box0CfgStart.setConfiguration(axis, obb0);

    intersectConfig box1CfgStart{};
    box1CfgStart.setConfiguration(axis, obb1);

    float t;
    float speed = relVelocity.dot(axis);

    if (box1CfgStart.m_max < box0CfgStart.m_min) // object1 left of object0
    {
        if(speed <= 0.0f)
            return false;
        t = (box0CfgStart.m_min - box1CfgStart.m_max) / speed;
        if (t > tfirst)
        {
            tfirst = t;
            side = intersectConfig::LEFT;
            box0Cfg = box0CfgStart;
            box1Cfg = box1CfgStart;
            s_contactNormal = axis;
        }

        t = (box0CfgStart.m_max - box1CfgStart.m_min) / speed;
        if (t < tlast)
            tlast = t;

        if (tfirst > tlast)
            return false;
    }
    else if (box0CfgStart.m_max < box1CfgStart.m_min)
    {
        if (speed >= 0.0f) // object1 moving away from object0
            return false;

        // find first time of contact on this axis
        t = (box0CfgStart.m_max - box1CfgStart.m_min) / speed;

        // If this is the new maximum first time of contact,  set side and
        // configuration.
        if (t > tfirst) {
            tfirst          = t;
            side            = intersectConfig::RIGHT;
            box0Cfg    = box0CfgStart;
            box1Cfg    = box1CfgStart;
            s_contactNormal = axis;
        }

        // find last time of contact on this axis
        t = (box0CfgStart.m_min - box1CfgStart.m_max) / speed;
        if (t < tlast)
            tlast = t;

        // quick out: intersection before desired interval
        if (tfirst > tlast)
            return false;
    }
    else
    {
        if (speed > 0.0f) {
            // find last time of contact on this axis
            t = (box0CfgStart.m_max - box1CfgStart.m_min) / speed;
            if (t < tlast) {
                tlast           = t;
                s_contactNormal = axis;
            }

            // quick out: intersection before desired interval
            if (tfirst > tlast)
                return false;
        } else if (speed < 0.0f) {
            // find last time of contact on this axis
            t = (box0CfgStart.m_min - box1CfgStart.m_max) / speed;
            if (t < tlast) {
                tlast           = t;
                s_contactNormal = axis;
            }

            // quick out: intersection before desired interval
            if (tfirst > tlast)
                return false;
        }
    }
    return true;
}

void CollisionDetect::FindContactSet(RigidBody *obb0, RigidBody *obb1, int side, intersectConfig &box0Cfg,
                                     intersectConfig &box1Cfg, float tfirst) {
    Eigen::Vector3f pts[8];
    int numPts = 0;
    
    RigidBody *box0Final, *box1Final;
    box0Final          = obb0;
    box1Final          = obb1;
    box0Final->x = obb0->x + tfirst * obb0->xdot;
    box1Final->x = obb1->x + tfirst * obb1->xdot;

    const int *b0Index = box0Cfg.m_index;
    const int *b1Index = box1Cfg.m_index;

    std::cout << "b0: ";
    for (int i = 0; i < 8; i++) {
        std::cout << b0Index[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "b1: ";
    for (int i = 0; i < 8; i++) {
        std::cout << b1Index[i] << " ";
    }
    std::cout << std::endl;
    
    if (side == intersectConfig::LEFT)
    {
        if (box0Cfg.m_map == intersectConfig::m1_1)
        {
            // box0point-box1face
            numPts = 1;
            pts[0] = GetPointFromIndex(box0Final, b0Index[0]);
            std::cout << "0 - point-face collision!" << std::endl;
        }else if (box1Cfg.m_map == intersectConfig::m1_1) {
            // box0face-box1point
            numPts = 1;
            pts[0] = GetPointFromIndex(box1Final, b1Index[7]);
            std::cout << "1 - face-point collision!" << std::endl;
        } else if (box0Cfg.m_map == intersectConfig::m2_2) {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0edge-box1edge intersection
                std::vector<Eigen::Vector3f> edge0(2), edge1(2);
                edge0[0] = GetPointFromIndex(box0Final, b0Index[0]);
                edge0[1] = GetPointFromIndex(box0Final, b0Index[1]);
                edge1[0] = GetPointFromIndex(box1Final, b1Index[6]);
                edge1[1] = GetPointFromIndex(box1Final, b1Index[7]);
                segmentSegment(edge0, edge1, numPts, pts);
                std::cout << "2 - edge-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0edge-box1face intersection
                std::vector<Eigen::Vector3f> edge0(2), face1(4);
                edge0[0] = GetPointFromIndex(box0Final, b0Index[0]);
                edge0[1] = GetPointFromIndex(box0Final, b0Index[1]);
                face1[0] = GetPointFromIndex(box1Final, b1Index[4]);
                face1[1] = GetPointFromIndex(box1Final, b1Index[5]);
                face1[2] = GetPointFromIndex(box1Final, b1Index[6]);
                face1[3] = GetPointFromIndex(box1Final, b1Index[7]);
                coplanarSegmentRectangle(edge0, face1, numPts, pts);
                std::cout << "3 - edge-face collision!" << std::endl;
            }
        } else
        {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0face-box1edge intersection
                std::vector<Eigen::Vector3f> face0(4), edge1(2);
                face0[0] = GetPointFromIndex(box0Final, b0Index[0]);
                face0[1] = GetPointFromIndex(box0Final, b0Index[1]);
                face0[2] = GetPointFromIndex(box0Final, b0Index[2]);
                face0[3] = GetPointFromIndex(box0Final, b0Index[3]);
                edge1[0] = GetPointFromIndex(box1Final, b1Index[6]);
                edge1[1] = GetPointFromIndex(box1Final, b1Index[7]);
                coplanarSegmentRectangle(edge1, face0, numPts, pts);
                std::cout << "4 - face-edge collision!" << std::endl;
            } else {
                // box0face-box1face intersection
                std::vector<Eigen::Vector3f> face0(4), face1(4);
                face0[0] = GetPointFromIndex(box0Final, b0Index[0]);
                face0[1] = GetPointFromIndex(box0Final, b0Index[1]);
                face0[2] = GetPointFromIndex(box0Final, b0Index[2]);
                face0[3] = GetPointFromIndex(box0Final, b0Index[3]);
                face1[0] = GetPointFromIndex(box1Final, b1Index[4]);
                face1[1] = GetPointFromIndex(box1Final, b1Index[5]);
                face1[2] = GetPointFromIndex(box1Final, b1Index[6]);
                face1[3] = GetPointFromIndex(box1Final, b1Index[7]);
                coplanarRectangleRectangle(face0, face1, numPts, pts);
                std::cout << "5 - face-face collision!" << std::endl;
            }
        }
    }else // side == RIGHT
    {
        // box1 on right of box0
        if (box0Cfg.m_map == intersectConfig::m1_1) {
            // box0point-box1face
            numPts = 1;
            pts[0] = GetPointFromIndex(box0Final, b0Index[7]);
            std::cout << "6 - point-face collision!" << std::endl;
        } else if (box1Cfg.m_map == intersectConfig::m1_1) {
            // box0face-box1point
            numPts = 1;
            pts[0] = GetPointFromIndex(box1Final, b1Index[0]);
            std::cout << "7 - face-point collision!" << std::endl;
        } else if (box0Cfg.m_map == intersectConfig::m2_2) {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0edge-box1edge intersection
                std::vector<Eigen::Vector3f> edge0(2), edge1(2);
                edge0[0] = GetPointFromIndex(box0Final, b0Index[6]);
                edge0[1] = GetPointFromIndex(box0Final, b0Index[7]);
                edge1[0] = GetPointFromIndex(box1Final, b1Index[0]);
                edge1[1] = GetPointFromIndex(box1Final, b1Index[1]);
                segmentSegment(edge0, edge1, numPts, pts);
                std::cout << "8 - edge-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0edge-box1face intersection
                std::vector<Eigen::Vector3f> edge0(2), face1(4);
                edge0[0] = GetPointFromIndex(box0Final, b0Index[6]);
                edge0[1] = GetPointFromIndex(box0Final, b0Index[7]);
                face1[0] = GetPointFromIndex(box1Final, b1Index[0]);
                face1[1] = GetPointFromIndex(box1Final, b1Index[1]);
                face1[2] = GetPointFromIndex(box1Final, b1Index[2]);
                face1[3] = GetPointFromIndex(box1Final, b1Index[3]);
                coplanarSegmentRectangle(edge0, face1, numPts, pts);
                std::cout << "9 - edge-face collision!" << std::endl;
            }
        } else // box0Cfg.m_map == m44
        {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0face-box1edge intersection
                std::vector<Eigen::Vector3f> face0(4), edge1(2);
                face0[0] = GetPointFromIndex(box0Final, b0Index[4]);
                face0[1] = GetPointFromIndex(box0Final, b0Index[5]);
                face0[2] = GetPointFromIndex(box0Final, b0Index[6]);
                face0[3] = GetPointFromIndex(box0Final, b0Index[7]);
                edge1[0] = GetPointFromIndex(box1Final, b1Index[0]);
                edge1[1] = GetPointFromIndex(box1Final, b1Index[1]);
                coplanarSegmentRectangle(edge1, face0, numPts, pts);
                std::cout << "10 - face-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0face-box1face intersection
                std::vector<Eigen::Vector3f> face0(4), face1(4);
                face0[0] = GetPointFromIndex(box0Final, b0Index[4]);
                face0[1] = GetPointFromIndex(box0Final, b0Index[5]);
                face0[2] = GetPointFromIndex(box0Final, b0Index[6]);
                face0[3] = GetPointFromIndex(box0Final, b0Index[7]);
                face1[0] = GetPointFromIndex(box1Final, b1Index[0]);
                face1[1] = GetPointFromIndex(box1Final, b1Index[1]);
                face1[2] = GetPointFromIndex(box1Final, b1Index[2]);
                face1[3] = GetPointFromIndex(box1Final, b1Index[3]);
                coplanarRectangleRectangle(face0, face1, numPts, pts);
                std::cout << "11 - face-face collision!" << std::endl;
            }
        }
    }

    s_contactNormal.normalized();
    Eigen::Vector3f direction= obb1->x - obb0->x;
    float diff = direction.dot(s_contactNormal);
    if (diff > 0.0f)
        s_contactNormal = -s_contactNormal;

    /*std::cout << "normal: " << s_contactNormal << std::endl
              << "Points: " << numPts << std::endl
              << "End" << std::endl
              << std::endl;*/

    for(int i = 0; i < numPts; i++)
    {
        Eigen::Vector3f p = pts[i];
        m_contacts.push_back( new Contact(obb0, obb1, p, s_contactNormal, s_minAxisPenetrationDepth) );
    }
}

Eigen::Vector3f CollisionDetect::GetPointFromIndex(RigidBody *obb, int index) {
    Box* box = dynamic_cast<Box*>(obb->geometry.get());
    Eigen::Vector3f point = obb->x;
    if (index & 4)
        point += box->halfDim(2) * obb->q.toRotationMatrix().col(2);
    else
        point -= box->halfDim(2) * obb->q.toRotationMatrix().col(2);
    if (index & 2)
        point += box->halfDim(1) * obb->q.toRotationMatrix().col(1);
    else
        point -= box->halfDim(1) * obb->q.toRotationMatrix().col(1);
    if (index & 1)
        point += box->halfDim(0) * obb->q.toRotationMatrix().col(0);
    else
        point -= box->halfDim(0) * obb->q.toRotationMatrix().col(0);
    return point;
}

void CollisionDetect::segmentSegment(const std::vector<Eigen::Vector3f>& segment0,
                                     const std::vector<Eigen::Vector3f>& segment1, int &numPts,
                                     Eigen::Vector3f *pts) {
    Eigen::Vector3f dir0 = segment0[1] - segment0[0];
    Eigen::Vector3f dir1 = segment1[1] - segment1[0];
    Eigen::Vector3f normal = dir0.cross(dir1);

    float sqrLenU = dir0.dot(dir0);
    float sqrLenV = dir1.dot(dir1);
    float sqrLenN = normal.dot(normal);

    if (sqrLenN < std::numeric_limits<float>::epsilon())
    {
        numPts = 2;
        pts[0] = segment0[0];
        pts[1] = segment0[1];

        // point 0
        Eigen::Vector3f V = segment1[0] - segment0[0];
        float c = V.dot(segment1[0]);
        clipConvexPolygonAgainstPlane(V, c, numPts, pts);

        // point 1
        V = -V;
        c = V.dot(segment1[1]);
        clipConvexPolygonAgainstPlane(V, c, numPts, pts);
    }else
        segmentThroughPlane(segment1, segment0[0], normal.cross(segment0[1] - segment0[0]), numPts, pts);
}

void CollisionDetect::coplanarRectangleRectangle(std::vector<Eigen::Vector3f> rectangle0,
                                                 std::vector<Eigen::Vector3f> rectangle1, int &numPts,
                                                 Eigen::Vector3f *pts) {
    numPts = 4;
    for (int i = 0; i < 4; ++i)
        pts[i] = rectangle0[i];
    for (int i0 = 3, i1 = 0; i1 < 4; i0 = i1++) {
        Eigen::Vector3f normal = rectangle1[i1] - rectangle1[i0];
        float constant = normal.dot(rectangle1[i0]);
        clipConvexPolygonAgainstPlane(normal, constant, numPts, pts);
    }
}

void CollisionDetect::coplanarSegmentRectangle(const std::vector<Eigen::Vector3f> &segment,
                                               const std::vector<Eigen::Vector3f> &rectangle, int &numPts,
                                               Eigen::Vector3f *pts) {
    numPts = 2;
    for (int i = 0; i < 2; ++i)
        pts[i] = segment[i];

    for (int i0 = 3, i1 = 0; i1 < 4; i0 = i1++) {
        Eigen::Vector3f normal = rectangle[i1] - rectangle[i0];
        float constant = normal.dot(rectangle[i0]);
        clipConvexPolygonAgainstPlane(normal, constant, numPts, pts);
    }
}

void CollisionDetect::clipConvexPolygonAgainstPlane(const Eigen::Vector3f &normal, float constant, int &numPts,
                                                    Eigen::Vector3f *pts) {
    int positive = 0, negative = 0, pIndex = -1;
    int currCount = numPts;

    float test[8];
    int i;
    for (i = 0; i < numPts; ++i)
    {
        test[i] = normal.dot(pts[i]) - constant +
                  abs(constant) * std::numeric_limits<float>::epsilon();

        if (test[i] >= 0.0) {
            ++positive;
            if (pIndex < 0) {
                pIndex = i;
            }
        } else {
            ++negative;
        }
    }
    if (numPts == 2) {
        // Lines are a little different, in that clipping the segment
        // cannot create a new segment, as clipping a polygon would
        if (positive > 0) {
            if (negative > 0) {
                int clip;

                if (pIndex == 0) {
                    // vertex0 positive, vertex1 is clipped
                    clip = 1;
                } else // pIndex == 1
                {
                    // vertex1 positive, vertex0 clipped
                    clip = 0;
                }

                float t  = test[pIndex] / (test[pIndex] - test[clip]);
                pts[clip].x() = pts[pIndex].x() + t * (pts[clip] - pts[pIndex]).x();
            }
            // otherwise both positive, no clipping
        } else {
            // Assert:  The entire line is clipped, but we should not
            // get here.
            numPts = 0;
        }
    }else {
        if (positive > 0) {
            if (negative > 0) {
                // plane transversely intersects polygon
                Eigen::Vector3f CV[8];
                int cCount = 0, cur, prv;
                float t;

                if (pIndex > 0) {
                    // first clip vertex on line
                    cur          = pIndex;
                    prv          = cur - 1;
                    t            = test[cur] / (test[cur] - test[prv]);
                    CV[cCount++] = pts[cur] + t * (pts[prv] - pts[cur]);

                    // vertices on positive side of line
                    while (cur < currCount && test[cur] >= 0.0) {
                        CV[cCount++] = pts[cur++];
                    }

                    // last clip vertex on line
                    if (cur < currCount) {
                        prv = cur - 1;
                    } else {
                        cur = 0;
                        prv = currCount - 1;
                    }
                    t            = test[cur] / (test[cur] - test[prv]);
                    CV[cCount++] = pts[cur] + t * (pts[prv] - pts[cur]);
                } else // pIndex is 0
                {
                    // vertices on positive side of line
                    cur = 0;
                    while (cur < currCount && test[cur] >= 0.0) {
                        CV[cCount++] = pts[cur++];
                    }

                    // last clip vertex on line
                    prv          = cur - 1;
                    t            = test[cur] / (test[cur] - test[prv]);
                    CV[cCount++] = pts[cur] + t * (pts[prv] - pts[cur]);

                    // skip vertices on negative side
                    while (cur < currCount && test[cur] < 0.0) {
                        cur++;
                    }

                    // first clip vertex on line
                    if (cur < currCount) {
                        prv          = cur - 1;
                        t            = test[cur] / (test[cur] - test[prv]);
                        CV[cCount++] = pts[cur] + t * (pts[prv] - pts[cur]);

                        // vertices on positive side of line
                        while (cur < currCount && test[cur] >= 0.0) {
                            CV[cCount++] = pts[cur++];
                        }
                    } else {
                        // cur = 0
                        prv          = currCount - 1;
                        t            = test[0] / (test[0] - test[prv]);
                        CV[cCount++] = pts[0] + t * (pts[prv] - pts[0]);
                    }
                }

                currCount = cCount;
                std::memcpy(pts, CV, cCount * sizeof(Eigen::Vector3f));
            }
            // else polygon fully on positive side of plane, nothing to do

            numPts = currCount;
        } else {
            // Polygon does not intersect positive side of plane, clip all.
            // This should not ever happen if called by the findintersect
            // routines after an intersection has been determined.
            numPts = 0;
        }
    }
}

void
CollisionDetect::segmentThroughPlane(const std::vector<Eigen::Vector3f> &segment, const Eigen::Vector3f &planeOrigin,
                                     const Eigen::Vector3f &planeNormal, int &numPts, Eigen::Vector3f *pts) {
    numPts     = 1;

    float u     = planeNormal.dot(planeOrigin);
    float v0    = planeNormal.dot(segment[0]);
    float v1    = planeNormal.dot(segment[1]);

    // Now that there it has been reduced to a 1-dimensional problem via
    // projection, it becomes easy to find the ratio along V that V
    // intersects with U.
    float ratio = (u - v0) / (v1 - v0);
    pts[0]     = segment[0] + ratio * (segment[1] - segment[0]);
}

bool
CollisionDetect::IsSeparated(float min0, float max0, float min1, float max1, float speed, float tmax, float &tlast) {
    s_minAxisPenetrationDepth = std::numeric_limits<float>::max();
    float invSpeed, t;
    if (max1 < min0) // box1 initially on left of box0
    {
        if (speed <= 0.0f)
            // The projection intervals are moving apart.
            return true;

        invSpeed = 1.0f / speed;

        t        = (min0 - max1) * invSpeed;

        if (t > s_contactTime)
            s_contactTime = t;

        if (s_contactTime > tmax)
            // intervals do not intersect during the specified time.
            return true;

        t = (max0 - min1) * invSpeed;

        if (t < tlast)
            tlast = t;

        if (s_contactTime > tlast)
            // Physically inconsistent times--the objects cannot intersect.
            return true;

        float currPenetration = max1 - min0;
        if (currPenetration < s_minAxisPenetrationDepth)
            s_minAxisPenetrationDepth = currPenetration;
    } else if (max0 < min1) // box1 initially on right of box0
    {
        if (speed >= 0.0f)
            // The projection intervals are moving apart.
            return true;

        invSpeed = 1.0f / speed;

        t        = (max0 - min1) * invSpeed;

        if (t > s_contactTime)
            s_contactTime = t;

        if (s_contactTime > tmax)
            // intervals do not intersect during the specified time.
            return true;

        t = (min0 - max1) * invSpeed;

        if (t < tlast)
            tlast = t;

        if (s_contactTime > tlast)
            // Physically inconsistent times--the objects cannot intersect.
            return true;

        float currPenetration = max0 - min1;
        if (currPenetration < s_minAxisPenetrationDepth)
            s_minAxisPenetrationDepth = currPenetration;
    } else // box0 and box1 initially overlap
    {
        if (speed > 0.0f) {
            t = (max0 - min1) / speed;

            if (t < tlast)
                tlast = t;

            if (s_contactTime > tlast)
                // Physically inconsistent times--the objects cannot intersect.
                return true;

            float currPenetration = max1 - min0;
            if (currPenetration < s_minAxisPenetrationDepth)
                s_minAxisPenetrationDepth = currPenetration;
        } else if (speed < 0.0f) {
            t = (min0 - max1) / speed;

            if (t < tlast)
                tlast = t;

            if (s_contactTime > tlast)
                // Physically inconsistent times--the objects cannot intersect.
                return true;

            float currPenetration = max0 - min1;
            if (currPenetration < s_minAxisPenetrationDepth)
                s_minAxisPenetrationDepth = currPenetration;
        }
    }
    return false;
}

bool CollisionDetect::DynamicCheck(RigidBody *obb0, RigidBody *obb1, float tmax) {

    Box* box0 = dynamic_cast<Box*>(obb0->geometry.get());
    Box* box1 = dynamic_cast<Box*>(obb1->geometry.get());

    Eigen::Vector3f velocity0 = obb0->xdot * 0.0167f;
    Eigen::Vector3f velocity1 = obb1->xdot * 0.0167f;
    
    if (velocity0 == velocity1) {
        s_contactTime = 0.0f;
        //return StaticCheck(obb0, obb1);
        return false;
    }

    const float cutoff                  = 1.0f - std::numeric_limits<float>::epsilon();
    bool existsParallelPair           = false;

    Eigen::Vector3f m_axes1[3],m_axes2[3];
    m_axes1[0] = obb0->q.toRotationMatrix() * Eigen::Vector3f(1,0,0);
    m_axes1[1] = obb0->q.toRotationMatrix() * Eigen::Vector3f(0,1,0);
    m_axes1[2] = obb0->q.toRotationMatrix() * Eigen::Vector3f(0,0,1);

    m_axes2[0] = obb1->q.toRotationMatrix() * Eigen::Vector3f(1,0,0);
    m_axes2[1] = obb1->q.toRotationMatrix() * Eigen::Vector3f(0,1,0);
    m_axes2[2] = obb1->q.toRotationMatrix() * Eigen::Vector3f(0,0,1);

    // convenience variables
    const Eigen::Vector3f *A = m_axes1;
    const Eigen::Vector3f *B = m_axes2;
    const Eigen::Vector3f EA = box0->halfDim;
    const Eigen::Vector3f EB = box1->halfDim;
    Eigen::Vector3f D        = obb1->x - obb0->x;
    Eigen::Vector3f W        = velocity1 - velocity0;
    float C[3][3];    // matrix C = A^T B, c_{ij} = Dot(A_i,B_j)
    float AbsC[3][3]; // |c_{ij}|
    float AD[3];      // Dot(A_i,D)
    float AW[3];      // Dot(A_i,W)
    float min0, max0, min1, max1, center, radius, speed;
    unsigned i, j;

    float tlast = std::numeric_limits<float>::max();

    // axes C0+t*A[i]
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            C[i][j]    = A[i].dot(B[j]);
            AbsC[i][j] = abs(C[i][j]);
            if (AbsC[i][j] > cutoff)
                existsParallelPair = true;
        }
        AD[i]  = A[i].dot(D);
        AW[i]  = A[i].dot(W);
        min0   = -EA[i];
        max0   = +EA[i];
        radius = EB[0] * AbsC[i][0] + EB[1] * AbsC[i][1] + EB[2] * AbsC[i][2];
        min1   = AD[i] - radius;
        max1   = AD[i] + radius;
        speed  = AW[i];
        if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
            return false;
    }

    // axes C0+t*B[i]
    for (i = 0; i < 3; ++i) {
        radius = EA[0] * AbsC[0][i] + EA[1] * AbsC[1][i] + EA[2] * AbsC[2][i];
        min0   = -radius;
        max0   = +radius;
        center = B[i].dot(D);
        min1   = center - EB[i];
        max1   = center + EB[i];
        speed  = W.dot(B[i]);
        if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
            return false;
    }

    // At least one pair of box axes was parallel, so the separation is
    // effectively in 2D where checking the "edge" normals is sufficient for
    // the separation of the boxes.
    if (existsParallelPair)
        return true;

    // axis C0+t*A0xB0
    radius = EA[1] * AbsC[2][0] + EA[2] * AbsC[1][0];
    min0   = -radius;
    max0   = +radius;
    center = AD[2] * C[1][0] - AD[1] * C[2][0];
    radius = EB[1] * AbsC[0][2] + EB[2] * AbsC[0][1];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[2] * C[1][0] - AW[1] * C[2][0];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A0xB1
    radius = EA[1] * AbsC[2][1] + EA[2] * AbsC[1][1];
    min0   = -radius;
    max0   = +radius;
    center = AD[2] * C[1][1] - AD[1] * C[2][1];
    radius = EB[0] * AbsC[0][2] + EB[2] * AbsC[0][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[2] * C[1][1] - AW[1] * C[2][1];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A0xB2
    radius = EA[1] * AbsC[2][2] + EA[2] * AbsC[1][2];
    min0   = -radius;
    max0   = +radius;
    center = AD[2] * C[1][2] - AD[1] * C[2][2];
    radius = EB[0] * AbsC[0][1] + EB[1] * AbsC[0][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[2] * C[1][2] - AW[1] * C[2][2];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A1xB0
    radius = EA[0] * AbsC[2][0] + EA[2] * AbsC[0][0];
    min0   = -radius;
    max0   = +radius;
    center = AD[0] * C[2][0] - AD[2] * C[0][0];
    radius = EB[1] * AbsC[1][2] + EB[2] * AbsC[1][1];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[0] * C[2][0] - AW[2] * C[0][0];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A1xB1
    radius = EA[0] * AbsC[2][1] + EA[2] * AbsC[0][1];
    min0   = -radius;
    max0   = +radius;
    center = AD[0] * C[2][1] - AD[2] * C[0][1];
    radius = EB[0] * AbsC[1][2] + EB[2] * AbsC[1][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[0] * C[2][1] - AW[2] * C[0][1];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A1xB2
    radius = EA[0] * AbsC[2][2] + EA[2] * AbsC[0][2];
    min0   = -radius;
    max0   = +radius;
    center = AD[0] * C[2][2] - AD[2] * C[0][2];
    radius = EB[0] * AbsC[1][1] + EB[1] * AbsC[1][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[0] * C[2][2] - AW[2] * C[0][2];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A2xB0
    radius = EA[0] * AbsC[1][0] + EA[1] * AbsC[0][0];
    min0   = -radius;
    max0   = +radius;
    center = AD[1] * C[0][0] - AD[0] * C[1][0];
    radius = EB[1] * AbsC[2][2] + EB[2] * AbsC[2][1];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[1] * C[0][0] - AW[0] * C[1][0];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A2xB1
    radius = EA[0] * AbsC[1][1] + EA[1] * AbsC[0][1];
    min0   = -radius;
    max0   = +radius;
    center = AD[1] * C[0][1] - AD[0] * C[1][1];
    radius = EB[0] * AbsC[2][2] + EB[2] * AbsC[2][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[1] * C[0][1] - AW[0] * C[1][1];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    // axis C0+t*A2xB2
    radius = EA[0] * AbsC[1][2] + EA[1] * AbsC[0][2];
    min0   = -radius;
    max0   = +radius;
    center = AD[1] * C[0][2] - AD[0] * C[1][2];
    radius = EB[0] * AbsC[2][1] + EB[1] * AbsC[2][0];
    min1   = center - radius;
    max1   = center + radius;
    speed  = AW[1] * C[0][2] - AW[0] * C[1][2];
    if (IsSeparated(min0, max0, min1, max1, speed, tmax, tlast))
        return false;

    return true;
}









