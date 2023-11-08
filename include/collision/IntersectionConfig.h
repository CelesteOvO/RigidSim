#ifndef RIGIDBODYTUTORIAL_INTERSECTIONCONFIG_H
#define RIGIDBODYTUTORIAL_INTERSECTIONCONFIG_H

#include <Eigen/Dense>
#include <vector>
#include "rigidbody/RigidBody.h"

class intersectConfig {
public:
    // ContactSide (order of the intervals of projection).
    enum { LEFT, RIGHT, NONE };

    // VertexProjectionMap (how the vertices are projected to the minimum
    // and maximum points of the interval).
    enum {
        m2,
        m11, // segments
        m3,
        m21,
        m12,
        m111, // triangles
        m44,
        m2_2,
        m1_1 // boxes
    };

    // The VertexProjectionMap value for the configuration.
    int m_map;

    // The order of the vertices.
    int m_index[8];

    // Projection interval.
    float m_min, m_max;

    /**
     *  \brief		sets up the configuration based on passed in axis and obb
     *	\Param		const Math::Vector3D<Real> & axis
     *	\Param		const Obb * obb
     *	\return		void
     */
    void setConfiguration(const Eigen::Vector3f &axis, const RigidBody *obb);
};
#endif //RIGIDBODYTUTORIAL_INTERSECTIONCONFIG_H
