#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

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
    if(!TestSAT(body0, body1)) {
        //std::cout << "There is a collision between box and box" << std::endl;

        ///Generate contacts
        DeriveContacts(body0, body1);
    }
}

bool CollisionDetect::TestSAT(RigidBody *body0, RigidBody *body1) {
    std::vector<Eigen::Vector3f> corners0 = getCorners(body0);
    std::vector<Eigen::Vector3f> corners1 = getCorners(body1);

    std::vector<Eigen::Vector3f> aNormals = getNormal(body0, corners0);
    std::vector<Eigen::Vector3f> bNormals = getNormal(body1, corners1);
    std::vector<Eigen::Vector3f> edgeNormals = getEdgeNormal(aNormals, bNormals);

    bool isSeparated = false;

    float minAxisPenetrationDepth = std::numeric_limits<float>::max();

    // For each normal
    for (const auto & aNormal : aNormals)
    {
        // Get the Min and Max projections for each object along the normal.
        float aMin, aMax;
        GetMinMax(body0, aNormal, aMin, aMax);

        float bMin, bMax;
        GetMinMax(body1, aNormal, bMin, bMax);

        isSeparated = aMax < bMin || bMax < aMin;
        float penetrationDepth = std::min(aMax - bMin, bMax - aMin);
        if(isSeparated)
            break;
        else if (penetrationDepth < minAxisPenetrationDepth)
        {
            minAxisPenetrationDepth = penetrationDepth;
        }
    }

    if (!isSeparated)
    {
        for (const auto & bNormal : bNormals)
        {
            // Get the Min and Max projections for each object along the normal.
            float aMin, aMax;
            GetMinMax(body0, bNormal, aMin, aMax);

            float bMin, bMax;
            GetMinMax(body1, bNormal, bMin, bMax);

            isSeparated = aMax < bMin || bMax < aMin;
            float penetrationDepth = std::min(aMax - bMin, bMax - aMin);
            if (isSeparated)
                break;
            else if (penetrationDepth < minAxisPenetrationDepth)
            {
                minAxisPenetrationDepth = penetrationDepth;
            }
        }
    }

if (!isSeparated)
    {
        for (const auto & edgeNormal : edgeNormals)
        {
            // Get the Min and Max projections for each object along the normal.
            float aMin, aMax;
            GetMinMax(body0, edgeNormal, aMin, aMax);

            float bMin, bMax;
            GetMinMax(body1, edgeNormal, bMin, bMax);

            isSeparated = aMax < bMin || bMax < aMin;
            float penetrationDepth = std::min(aMax - bMin, bMax - aMin);
            if (isSeparated)
                break;
            else if (penetrationDepth < minAxisPenetrationDepth)
            {
                minAxisPenetrationDepth = penetrationDepth;
            }
        }
    }
    s_minAxisPenetrationDepth = minAxisPenetrationDepth;
    return isSeparated;
}

std::vector<Eigen::Vector3f> CollisionDetect::getNormal(RigidBody *body, std::vector<Eigen::Vector3f> corners) {
    std::vector<Eigen::Vector3f> normals;

    Eigen::Vector3f U = corners[1] - corners[0];
    Eigen::Vector3f V = corners[2] - corners[0];
    Eigen::Vector3f N = U.cross(V).normalized();
    normals.push_back(N);

    U = corners[3] - corners[2];
    V = corners[6] - corners[2];
    N = U.cross(V).normalized();
    normals.push_back(N);

    U = corners[2] - corners[0];
    V = corners[4] - corners[0];
    N = U.cross(V).normalized();
    normals.push_back(N);

    return normals;
}

std::vector<Eigen::Vector3f> CollisionDetect::getCorners(RigidBody *body) {
    Box* box = dynamic_cast<Box*>(body->geometry.get());
    std::vector<Eigen::Vector3f> corners;
    for(const auto & corner : box->corners)
    {
        corners.emplace_back(body->q.toRotationMatrix() * corner + body->x );
    }
    return corners;
}

std::vector<Eigen::Vector3f>
CollisionDetect::getEdgeNormal(const std::vector<Eigen::Vector3f>& aNormals, const std::vector<Eigen::Vector3f>& bNormals) {
    std::vector<Eigen::Vector3f> edgeNormals;
    for (const auto & aNormal : aNormals)
    {
        for (const auto & bNormal : bNormals)
        {
            Eigen::Vector3f edgeNormal = aNormal.cross(bNormal).normalized();
            edgeNormals.push_back(edgeNormal);
        }
    }
    return edgeNormals;
}

void CollisionDetect::GetMinMax(RigidBody *obj, const Eigen::Vector3f& axis, float &min, float &max) {
    std::vector<Eigen::Vector3f> corners = getCorners(obj);
    min = axis.dot(corners[0]);
    max = min;
    for (const auto & corner : corners)
    {
        float projection = axis.dot(corner);
        if (projection < min)
            min = projection;
        if (projection > max)
            max = projection;
    }
}

void CollisionDetect::DeriveContacts(RigidBody *body0, RigidBody *body1) {
    int side = intersectConfig::NONE;
    intersectConfig box0Cfg{}, box1Cfg{};

    std::vector<Eigen::Vector3f> corners0 = getCorners(body0);
    std::vector<Eigen::Vector3f> corners1 = getCorners(body1);

    std::vector<Eigen::Vector3f> aNormals = getNormal(body0, corners0);
    std::vector<Eigen::Vector3f> bNormals = getNormal(body1, corners1);
    std::vector<Eigen::Vector3f> edgeNormals = getEdgeNormal(aNormals, bNormals);

    for (const auto & aNormal : aNormals)
    {
        if(!FindintersectionOnAxis(body0, body1, aNormal,side, box0Cfg, box1Cfg))
            return;
    }

    for (const auto & bNormal : bNormals)
    {
        if(!FindintersectionOnAxis(body0, body1, bNormal,side, box0Cfg, box1Cfg))
            return;
    }

    for (const auto & edgeNormal : edgeNormals)
    {
        if (abs(edgeNormal.dot(edgeNormal)) <= std::numeric_limits<float>::epsilon())// skip almost parallel edges
        {
            FindContactSet(body0, body1, side, box0Cfg, box1Cfg);
        }
        if(!FindintersectionOnAxis(body0, body1, edgeNormal,side, box0Cfg, box1Cfg))
            return;
    }

    if (side == intersectConfig::NONE)
        return;

    FindContactSet(body0, body1, side, box0Cfg, box1Cfg);
}

bool CollisionDetect::FindintersectionOnAxis(RigidBody *body0, RigidBody *body1, const Eigen::Vector3f &axis, int &side,
                                             intersectConfig &box0Cfg, intersectConfig &box1Cfg) {
    intersectConfig box0CfgStart{};
    box0CfgStart.setConfiguration(axis, body0);

    intersectConfig box1CfgStart{};
    box1CfgStart.setConfiguration(axis, body1);

    if (box1CfgStart.m_max < box0CfgStart.m_min) // 1在0的左边
    {
        if (box0CfgStart.m_min - box1CfgStart.m_max > s_minAxisPenetrationDepth)
            return false;
        else
        {
            s_contactNormal = axis;
            side = intersectConfig::LEFT;
            box0Cfg = box0CfgStart;
            box1Cfg = box1CfgStart;
        }
    }
    else if (box0CfgStart.m_max < box1CfgStart.m_min) // 0在1的左边
    {
        if (box1CfgStart.m_min - box0CfgStart.m_max > s_minAxisPenetrationDepth)
            return false;
        else
        {
            s_contactNormal = -axis;
            side = intersectConfig::RIGHT;
            box0Cfg = box0CfgStart;
            box1Cfg = box1CfgStart;
        }
    }
    else // 0和1相交
    {
        if (box0CfgStart.m_min < box1CfgStart.m_min)
        {
            s_contactNormal = -axis;
            side = intersectConfig::RIGHT;
            box0Cfg = box0CfgStart;
            box1Cfg = box1CfgStart;
        }
        else
        {
            s_contactNormal = axis;
            side = intersectConfig::LEFT;
            box0Cfg = box0CfgStart;
            box1Cfg = box1CfgStart;
        }
    }

    return true;
}

void CollisionDetect::FindContactSet(RigidBody *body0, RigidBody *body1, int side, const intersectConfig &box0Cfg,
                                     const intersectConfig &box1Cfg) {
    Eigen::Vector3f pts[8];
    int numPts;
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

    std::cout << "Contact: ";

    if (side == intersectConfig::LEFT) {
        // box1 on left of box0
        if (box0Cfg.m_map == intersectConfig::m1_1) {
            // box0point-box1face
            numPts = 1;
            pts[0] = GetPointFromIndex(body0, b0Index[0]);
            std::cout << "0 - point-face collision!" << std::endl;
        } else if (box1Cfg.m_map == intersectConfig::m1_1) {
            // box0face-box1point
            numPts = 1;
            pts[0] = GetPointFromIndex(body1, b1Index[7]);
            std::cout << "1 - face-point collision!" << std::endl;
        } else if (box0Cfg.m_map == intersectConfig::m2_2) {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0edge-box1edge intersection
                Eigen::Vector3f edge0[2], edge1[2];
                edge0[0] = GetPointFromIndex(body0, b0Index[0]);
                edge0[1] = GetPointFromIndex(body0, b0Index[1]);
                edge1[0] = GetPointFromIndex(body1, b1Index[6]);
                edge1[1] = GetPointFromIndex(body1, b1Index[7]);
                segmentSegment(edge0, edge1, numPts, pts);
                std::cout << "2 - edge-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0edge-box1face intersection
                Eigen::Vector3f edge0[2], face1[4];
                edge0[0] = GetPointFromIndex(body0, b0Index[0]);
                edge0[1] = GetPointFromIndex(body0, b0Index[1]);
                face1[0] = GetPointFromIndex(body1, b1Index[4]);
                face1[1] = GetPointFromIndex(body1, b1Index[5]);
                face1[2] = GetPointFromIndex(body1, b1Index[6]);
                face1[3] = GetPointFromIndex(body1, b1Index[7]);
                coplanarSegmentRectangle(edge0, face1, numPts, pts);
                std::cout << "3 - edge-face collision!" << std::endl;
            }
        } else // box0Cfg.m_map == m44
        {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0face-box1edge intersection
                Eigen::Vector3f face0[4], edge1[2];
                face0[0] = GetPointFromIndex(body0, b0Index[0]);
                face0[1] = GetPointFromIndex(body0, b0Index[1]);
                face0[2] = GetPointFromIndex(body0, b0Index[2]);
                face0[3] = GetPointFromIndex(body0, b0Index[3]);
                edge1[0] = GetPointFromIndex(body1, b1Index[6]);
                edge1[1] = GetPointFromIndex(body1, b1Index[7]);
                coplanarSegmentRectangle(edge1, face0, numPts, pts);
                std::cout << "4 - face-edge collision!" << std::endl;
            } else {
                // box0face-box1face intersection
                Eigen::Vector3f face0[4], face1[4];
                face0[0] = GetPointFromIndex(body0, b0Index[0]);
                face0[1] = GetPointFromIndex(body0, b0Index[1]);
                face0[2] = GetPointFromIndex(body0, b0Index[2]);
                face0[3] = GetPointFromIndex(body0, b0Index[3]);
                face1[0] = GetPointFromIndex(body1, b1Index[4]);
                face1[1] = GetPointFromIndex(body1, b1Index[5]);
                face1[2] = GetPointFromIndex(body1, b1Index[6]);
                face1[3] = GetPointFromIndex(body1, b1Index[7]);
                coplanarRectangleRectangle(face0, face1, numPts, pts);
                std::cout << "5 - face-face collision!" << std::endl;
            }
        }
    } else // side == RIGHT
    {
        // box1 on right of box0
        if (box0Cfg.m_map == intersectConfig::m1_1) {
            // box0point-box1face
            numPts = 1;
            pts[0] = GetPointFromIndex(body0, b0Index[7]);
            std::cout << "6 - point-face collision!" << std::endl;
        } else if (box1Cfg.m_map == intersectConfig::m1_1) {
            // box0face-box1point
            numPts = 1;
            pts[0] = GetPointFromIndex(body1, b1Index[0]);
            std::cout << "7 - face-point collision!" << std::endl;
        } else if (box0Cfg.m_map == intersectConfig::m2_2) {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0edge-box1edge intersection
                Eigen::Vector3f edge0[2], edge1[2];
                edge0[0] = GetPointFromIndex(body0, b0Index[6]);
                edge0[1] = GetPointFromIndex(body0, b0Index[7]);
                edge1[0] = GetPointFromIndex(body1, b1Index[0]);
                edge1[1] = GetPointFromIndex(body1, b1Index[1]);
                segmentSegment(edge0, edge1, numPts, pts);
                std::cout << "8 - edge-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0edge-box1face intersection
                Eigen::Vector3f edge0[2], face1[4];
                edge0[0] = GetPointFromIndex(body0, b0Index[6]);
                edge0[1] = GetPointFromIndex(body0, b0Index[7]);
                face1[0] = GetPointFromIndex(body1, b1Index[0]);
                face1[1] = GetPointFromIndex(body1, b1Index[1]);
                face1[2] = GetPointFromIndex(body1, b1Index[2]);
                face1[3] = GetPointFromIndex(body1, b1Index[3]);
                coplanarSegmentRectangle(edge0, face1, numPts, pts);
                std::cout << "9 - edge-face collision!" << std::endl;
            }
        } else // box0Cfg.m_map == m44
        {
            if (box1Cfg.m_map == intersectConfig::m2_2) {
                // box0face-box1edge intersection
                Eigen::Vector3f face0[4], edge1[2];
                face0[0] = GetPointFromIndex(body0, b0Index[4]);
                face0[1] = GetPointFromIndex(body0, b0Index[5]);
                face0[2] = GetPointFromIndex(body0, b0Index[6]);
                face0[3] = GetPointFromIndex(body0, b0Index[7]);
                edge1[0] = GetPointFromIndex(body1, b1Index[0]);
                edge1[1] = GetPointFromIndex(body1, b1Index[1]);
                coplanarSegmentRectangle(edge1, face0, numPts, pts);
                std::cout << "10 - face-edge collision!" << std::endl;
            } else // box1Cfg.m_map == m44
            {
                // box0face-box1face intersection
                Eigen::Vector3f face0[4], face1[4];
                face0[0] = GetPointFromIndex(body0, b0Index[4]);
                face0[1] = GetPointFromIndex(body0, b0Index[5]);
                face0[2] = GetPointFromIndex(body0, b0Index[6]);
                face0[3] = GetPointFromIndex(body0, b0Index[7]);
                face1[0] = GetPointFromIndex(body1, b1Index[0]);
                face1[1] = GetPointFromIndex(body1, b1Index[1]);
                face1[2] = GetPointFromIndex(body1, b1Index[2]);
                face1[3] = GetPointFromIndex(body1, b1Index[3]);
                coplanarRectangleRectangle(face0, face1, numPts, pts);
                std::cout << "11 - face-face collision!" << std::endl;
            }
        }
    }

    s_contactNormal.normalized();
    Eigen::Vector3f dir = body0->x - body1->x;
    if (dir.dot(s_contactNormal) < 0.0f)
        s_contactNormal = -s_contactNormal;
    for (int i = 0; i < numPts; ++i)
    {
        m_contacts.push_back( new Contact(body0, body1, pts[i], s_contactNormal, s_minAxisPenetrationDepth) );
    }
}

Eigen::Vector3f CollisionDetect::GetPointFromIndex(RigidBody *body, int index) {
    Box* box = dynamic_cast<Box*>(body->geometry.get());
    return body->q.toRotationMatrix() * box->corners[index] + body->x;
}

void CollisionDetect::segmentSegment(const Eigen::Vector3f *segment0, const Eigen::Vector3f *segment1, int &numPts,
                                     Eigen::Vector3f *pts) {
    Eigen::Vector3f dir0 = segment0[1] - segment0[0];
    Eigen::Vector3f dir1 = segment1[1] - segment1[0];

    Eigen::Vector3f normal = dir0.cross(dir1);
    
    float sqrLen0 = dir0.dot(dir0);
    float sqrLen1 = dir1.dot(dir1);
    float sqrLenN = normal.dot(normal);

    if (sqrLenN < std::numeric_limits<float>::epsilon() * sqrLen0 * sqrLen1)
        colinearSegments(segment0, segment1, numPts, pts);
    else
        segmentThroughPlane(segment1, segment0[0], normal.cross(segment0[1] - segment0[0]), numPts, pts);
}

void
CollisionDetect::coplanarSegmentRectangle(const Eigen::Vector3f *segment, const Eigen::Vector3f *rectangle, int &numPts,
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

void CollisionDetect::coplanarRectangleRectangle(const Eigen::Vector3f *rectangle0, const Eigen::Vector3f *rectangle1,
                                                 int &numPts, Eigen::Vector3f *pts) {
    numPts = 4;
    for (int i = 0; i < 4; ++i)
        pts[i] = rectangle0[i];
    for (int i0 = 3, i1 = 0; i1 < 4; i0 = i1++) {
        Eigen::Vector3f normal = rectangle1[i1] - rectangle1[i0];
        float constant = normal.dot(rectangle1[i0]);
        clipConvexPolygonAgainstPlane(normal, constant, numPts, pts);
    }
}

void CollisionDetect::colinearSegments(const Eigen::Vector3f *segment0, const Eigen::Vector3f *segment1, int &numPts,
                                       Eigen::Vector3f *pts) {
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
}

void CollisionDetect::segmentThroughPlane(const Eigen::Vector3f *segment, const Eigen::Vector3f &planeOrigin,
                                          const Eigen::Vector3f &planeNormal, int &numPts, Eigen::Vector3f *pts) {
    numPts     = 1;

    float u     = planeNormal.dot(planeOrigin);
    float v0    = planeNormal.dot(segment[0]);
    float v1    = planeNormal.dot(segment[1]);

    float ratio = (u - v0) / (v1 - v0);
    pts[0]     = segment[0] + ratio * (segment[1] - segment[0]);
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



















