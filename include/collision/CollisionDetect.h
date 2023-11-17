#pragma once

#include <Eigen/Dense>
#include <vector>
#include "IntersectionConfig.h"

class Contact;
class RigidBody;
class RigidBodySystem;

//
// The main collision detection class,
// it will test the geometries of each pair of rigid bodies
// and populate an array with contacts.
//
class CollisionDetect
{
public:

    // Constructor.
    //
    explicit CollisionDetect(RigidBodySystem* rigidBodySystem);

    // Tests for collisions between all pairs of bodies in the rigid body system
    // and generates contacts for any intersecting geometries.
    // The array of generated contacts can be retrieved by calling getContacts().
    //
    void detectCollisions();

    void clear();

    // Compute the Jacobians for contacts
    void computeContactJacobians();

    const std::vector<Contact*>& getContacts() const { return m_contacts; }
    std::vector<Contact*>& getContacts() { return m_contacts; }

private:

    // Sphere-sphere collision test.
    // Assumes that both @a body0 and @a body1 have a sphere collision geometry.
    //
    void collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1);

    // Sphere-box collision test.
    // Assumes that @a body0 has a sphere geometry and @a body1 has a box geometry.
    //
    void collisionDetectSphereBox(RigidBody* body0, RigidBody* body1);

    void collisionDetectBoxPlane(RigidBody* body0, RigidBody* body1);

    void collisionDetectSpherePlane(RigidBody* body0, RigidBody* body1);

    void collisionDetectBoxBox(RigidBody* body0, RigidBody* body1);

public:
    static bool TestSAT(RigidBody* body0, RigidBody* body1);
    static std::vector<Eigen::Vector3f> getNormal(RigidBody* body, std::vector<Eigen::Vector3f> corners);
    static std::vector<Eigen::Vector3f> getEdgeNormal(const std::vector<Eigen::Vector3f>& aNormals, const std::vector<Eigen::Vector3f>& bNormals);
    static std::vector<Eigen::Vector3f> getCorners(RigidBody* body);
    static void GetMinMax(RigidBody* obj, const Eigen::Vector3f& axis, float& min, float& max);
    void DeriveContacts(RigidBody *body0, RigidBody *body1);
    static bool FindintersectionOnAxis(RigidBody* body0, RigidBody* body1, const Eigen::Vector3f &axis, int &side, intersectConfig &box0Cfg, intersectConfig &box1Cfg);
    void FindContactSet(RigidBody* body0, RigidBody* body1, int side, const intersectConfig &box0Cfg, const intersectConfig &box1Cfg);

    static Eigen::Vector3f GetPointFromIndex(RigidBody* body, int index);
    void segmentSegment(const Eigen::Vector3f segment0[2], const Eigen::Vector3f segment1[2], int &numPts, Eigen::Vector3f *pts);
    static void coplanarSegmentRectangle(const Eigen::Vector3f segment[2], const Eigen::Vector3f rectangle[4], int &numPts, Eigen::Vector3f *pts);
    static void coplanarRectangleRectangle(const Eigen::Vector3f rectangle0[4], const Eigen::Vector3f rectangle1[4], int &numPts, Eigen::Vector3f *pts);

    static void colinearSegments(const Eigen::Vector3f segment0[2], const Eigen::Vector3f segment1[2], int &numPts, Eigen::Vector3f *pts);
    static void segmentThroughPlane(const Eigen::Vector3f segment[2], const Eigen::Vector3f &planeOrigin, const Eigen::Vector3f &planeNormal, int &numPts, Eigen::Vector3f *pts);
    static void clipConvexPolygonAgainstPlane(const Eigen::Vector3f &normal, float constant, int &numPts, Eigen::Vector3f *pts);

private:

    RigidBodySystem* m_rigidBodySystem;     // Rigid body system where collision detection is performed.
    std::vector<Contact*> m_contacts;       // Contact array.
    static Eigen::Vector3f s_contactNormal;
    static float s_minAxisPenetrationDepth;
};
