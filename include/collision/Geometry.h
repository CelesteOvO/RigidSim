#pragma once

#include <Eigen/Dense>
#include <Discregrid/All>
#include <utility>

// List of geometry type ids.
enum eGeometryType { kSphere, kBox , kPlane, kSDF };

// Generic geometry interface.
//
class Geometry
{
public:

    virtual ~Geometry() = default;

    virtual Eigen::Matrix3f computeInertia(float _mass) = 0;

    virtual eGeometryType getType() const  = 0;

protected:
    Eigen::Matrix3f m_I;          // Inertia 3x3 matrix for this. Only used for local computations. (internal)
};


// Sphere geometry.
//
class Sphere : public Geometry
{
public:
    float radius;           // Sphere radius.

    explicit Sphere(float _radius) : radius(_radius) {}
    ~Sphere() override = default;

    Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setZero();
        m_I(0,0) = m_I(1,1) = m_I(2,2) = (2.0f/5.0f) * _mass * radius * radius;
        return m_I;
    }

    eGeometryType getType() const override { return kSphere; }

};

// Box geometry.
//
class Box : public Geometry
{
public:
    Eigen::Vector3f dim;        // Box dimensions.

    explicit Box(Eigen::Vector3f  _dim) : dim(std::move(_dim)) {

    }
    ~Box() override = default;

    Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setZero();
        m_I(0,0) = (1.0f/12.0f)*_mass*(dim[1]*dim[1] + dim[2]*dim[2]);
        m_I(1,1) = (1.0f/12.0f)*_mass*(dim[0]*dim[0] + dim[2]*dim[2]);
        m_I(2,2) = (1.0f/12.0f)*_mass*(dim[0]*dim[0] + dim[1]*dim[1]);
        return m_I;
    }

    eGeometryType getType() const override { return kBox; }

};

// Plane geometry.
class Plane : public Geometry
{
public:
    Eigen::Vector3f normal;         // The plane normal.

    explicit Plane(Eigen::Vector3f  _normal)
            : normal(std::move(_normal)) {}
    ~Plane() override = default;

    Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setIdentity();
        return _mass*m_I;
    }

    eGeometryType getType() const override { return kPlane; }
};

// SDF geometry.
class SDFGeometry : public Geometry
{
protected:
    Eigen::Matrix3f I0;
    // Compute the inertia of the mesh.
    // Assumes a manifold mesh and center of mass is the geometry origin.
    // Based on code from: http://melax.github.io/volint.html (Stan Melax)
    void computeInertia();
public:
    std::unique_ptr<Discregrid::TriangleMesh> mesh;
    std::unique_ptr<Discregrid::CubicLagrangeDiscreteGrid> sdf;

    SDFGeometry(const std::string& filename, const std::array<unsigned int, 3>& resolution);
    ~SDFGeometry() override = default;

    virtual Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I = _mass * I0;
        return m_I;
    }

    virtual eGeometryType getType() const override { return kSDF; }
};

