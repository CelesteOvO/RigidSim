#pragma once

#include <Eigen/Dense>
#include <vector>

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

    static int testSAT(float minOverlap, std::vector<Eigen::Vector3f> normals, RigidBody* body1, RigidBody* body2, const std::vector<Eigen::Vector3f>& corners1, const std::vector<Eigen::Vector3f>& corners2);
    static float transformToAxis(RigidBody *pBody, const Eigen::Vector3f& axis, const std::vector<Eigen::Vector3f>& corners);

    void collisionDetectFaceVertex(RigidBody* body0, RigidBody* body1, Eigen::Vector3f axis, float penetration);
    void collisionDetectEdgeEdge(RigidBody* body0, RigidBody* body1, Eigen::Vector3f axis, float penetration, int oneAxisIndex, int twoAxisIndex);

private:

    RigidBodySystem* m_rigidBodySystem;     // Rigid body system where collision detection is performed.
    std::vector<Contact*> m_contacts;       // Contact array.

};
