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

    static std::vector<Eigen::Vector3f> getFaceNormal(RigidBody* body, std::vector<Eigen::Vector3f> corners);
    static std::vector<Eigen::Vector3f> getEdgeNormal(RigidBody* body1, RigidBody* body2, std::vector<Eigen::Vector3f> corners1, std::vector<Eigen::Vector3f> corners2);
    static std::vector<Eigen::Vector3f> getCorners(RigidBody* body);

    bool CollisionDetect::DeriveContacts(RigidBody *body0, RigidBody *body1);
    static bool FindintersectionOnAxis(RigidBody* obb0, RigidBody* obb1, Eigen::Vector3f &axis, Eigen::Vector3f &relVelocity, float &tfirst, float &tlast, int &side, intersectConfig &box0Cfg,  intersectConfig &box1Cfg);
    void FindContactSet(RigidBody* obb0, RigidBody* obb1, int side, intersectConfig &box0Cfg, intersectConfig &box1Cfg, float tfirst);
    static Eigen::Vector3f GetPointFromIndex(RigidBody* body, int index);
    static bool IsSeparated(float min0, float max0, float min1, float max1, float speed, float tmax, float &tlast);
    static bool DynamicCheck(RigidBody* body0, RigidBody* body1, float tmax);
    static void segmentSegment(const std::vector<Eigen::Vector3f>& segment0, const std::vector<Eigen::Vector3f>& segment1, int &numPts, Eigen::Vector3f *pts);
    static void coplanarSegmentRectangle(const std::vector<Eigen::Vector3f>& segment, const std::vector<Eigen::Vector3f>& rectangle, int &numPts, Eigen::Vector3f *pts);
    static void coplanarRectangleRectangle(std::vector<Eigen::Vector3f> rectangle0, std::vector<Eigen::Vector3f> rectangle1, int &numPts, Eigen::Vector3f *pts);
    static void clipConvexPolygonAgainstPlane(const Eigen::Vector3f& normal, float constant, int &numPts, Eigen::Vector3f *pts);
    static void segmentThroughPlane(const std::vector<Eigen::Vector3f>& segment, const Eigen::Vector3f& planeOrigin, const Eigen::Vector3f& planeNormal, int &numPts, Eigen::Vector3f *pts);
private:

    RigidBodySystem* m_rigidBodySystem;     // Rigid body system where collision detection is performed.
    std::vector<Contact*> m_contacts;       // Contact array.
    static float s_contactTime;
    static Eigen::Vector3f s_contactNormal;
    static float s_minAxisPenetrationDepth;
};
