//
// Created by LiYifan on 2023/11/8.
//
#include "IntersectionConfig.h"

void intersectConfig::setConfiguration(const Eigen::Vector3f &axis, RigidBody *obb) {
    Box* box = dynamic_cast<Box*>(obb->geometry.get());
    float axes[3] = {axis.dot(obb->q.toRotationMatrix().col(0)), axis.dot(obb->q.toRotationMatrix().col(1)), axis.dot(obb->q.toRotationMatrix().col(2))};
    float absAxes[3] = {std::abs(axes[0]), std::abs(axes[1]), std::abs(axes[2])};
    
    float maxProjectedExtent;

    if (absAxes[0] < std::numeric_limits<float>::epsilon())
    {
        if(absAxes[1] < std::numeric_limits<float>::epsilon())
        {
            m_map              = m44;
            maxProjectedExtent = absAxes[2] * box->dim(2) * 0.5f;
            
            // faces have normals along axis[2]
            if (axes[2] > 0.0f)
            {
                m_index[0] = 0;
                m_index[1] = 1;
                m_index[2] = 3;
                m_index[3] = 2;

                m_index[4] = 6;
                m_index[5] = 7;
                m_index[6] = 5;
                m_index[7] = 4;
            }
            else
            {
                m_index[0] = 6;
                m_index[1] = 7;
                m_index[2] = 5;
                m_index[3] = 4;

                m_index[4] = 0;
                m_index[5] = 1;
                m_index[6] = 3;
                m_index[7] = 2;
            }
        }else if(absAxes[2] < std::numeric_limits<float>::epsilon()) {
            // face-face
            m_map = m44;
            maxProjectedExtent = absAxes[1] * box->dim(1) * 0.5f;

            // faces have normals along axis[1]
            if (axes[1] > 0.0f) {
                m_index[0] = 4;
                m_index[1] = 5;
                m_index[2] = 1;
                m_index[3] = 0;

                m_index[4] = 2;
                m_index[5] = 3;
                m_index[6] = 7;
                m_index[7] = 6;
            } else {
                m_index[0] = 2;
                m_index[1] = 3;
                m_index[2] = 7;
                m_index[3] = 6;

                m_index[4] = 4;
                m_index[5] = 5;
                m_index[6] = 1;
                m_index[7] = 0;
            }
        } else // only axes[0] is equal to 0
        {
            // seg-seg
            m_map              = m2_2;

            maxProjectedExtent = absAxes[1] * box->dim(1) * 0.5f + absAxes[2] * box->dim(2) * 0.5f;

            // axis 0 is perpendicular to axis
            if (axes[1] > 0.0f) {
                if (axes[2] > 0.0f) {
                    m_index[0] = 0;
                    m_index[1] = 1;

                    m_index[6] = 6;
                    m_index[7] = 7;
                } else {
                    m_index[0] = 4;
                    m_index[1] = 5;

                    m_index[6] = 2;
                    m_index[7] = 3;
                }
            } else // axes[1] < 0
            {
                if (axes[2] > 0.0f) {
                    m_index[0] = 2;
                    m_index[1] = 3;

                    m_index[6] = 4;
                    m_index[7] = 5;
                } else {
                    m_index[0] = 6;
                    m_index[1] = 7;

                    m_index[6] = 0;
                    m_index[7] = 1;
                }
            }
        }
    }else if (absAxes[1] < std::numeric_limits<float>::epsilon()) {
        if (absAxes[2] < std::numeric_limits<float>::epsilon()) {
            // face-face
            m_map              = m44;

            maxProjectedExtent = absAxes[0] * box->dim(0) * 0.5f;

            // faces have normals along axis[0]
            if (axes[0] > 0.0f) {
                m_index[0] = 0;
                m_index[1] = 2;
                m_index[2] = 6;
                m_index[3] = 4;

                m_index[4] = 5;
                m_index[5] = 7;
                m_index[6] = 3;
                m_index[7] = 1;
            } else {
                m_index[4] = 0;
                m_index[5] = 2;
                m_index[6] = 6;
                m_index[7] = 4;

                m_index[0] = 5;
                m_index[1] = 7;
                m_index[2] = 3;
                m_index[3] = 1;
            }

        } else // only axes[1] is equal to 0
        {
            // seg-seg
            m_map              = m2_2;

            maxProjectedExtent = absAxes[0] * box->dim(0) * 0.5f + absAxes[2] * box->dim(2) * 0.5f;

            // axis 1 is perpendicular to axis
            if (axes[0] > 0.0f) {
                if (axes[2] > 0.0f) {
                    m_index[0] = 0;
                    m_index[1] = 2;

                    m_index[6] = 5;
                    m_index[7] = 7;
                } else {
                    m_index[0] = 4;
                    m_index[1] = 6;

                    m_index[6] = 1;
                    m_index[7] = 3;
                }
            } else // axes[0] < 0
            {
                if (axes[2] > 0.0f) {
                    m_index[0] = 1;
                    m_index[1] = 3;

                    m_index[6] = 4;
                    m_index[7] = 6;
                } else {
                    m_index[0] = 5;
                    m_index[1] = 7;

                    m_index[6] = 0;
                    m_index[7] = 2;
                }
            }
        }
    } else if (absAxes[2] < std::numeric_limits<float>::epsilon()) {
        // only axis2 less than zero
        // seg-seg
        m_map              = m2_2;

        maxProjectedExtent = absAxes[0] * box->dim(0) * 0.5f + absAxes[1] * box->dim(1) * 0.5f;

        // axis 2 is perpendicular to axis
        if (axes[0] > 0.0f) {
            if (axes[1] > 0.0f) {
                m_index[0] = 0;
                m_index[1] = 4;

                m_index[6] = 3;
                m_index[7] = 7;
            } else {
                m_index[0] = 2;
                m_index[1] = 6;

                m_index[6] = 1;
                m_index[7] = 5;
            }
        } else // axes[0] < 0
        {
            if (axes[1] > 0.0f) {
                m_index[0] = 1;
                m_index[1] = 5;

                m_index[6] = 2;
                m_index[7] = 6;
            } else {
                m_index[0] = 3;
                m_index[1] = 7;

                m_index[6] = 0;
                m_index[7] = 4;
            }
        }
    }
    else // no axis is equal to zero
    {
        // point-point (unique maximal and minimal vertex)
        m_map              = m1_1;

        maxProjectedExtent = absAxes[0] * box->dim(0) * 0.5f + absAxes[1] * box->dim(1) * 0.5f +
                             absAxes[2] * box->dim(2) * 0.5f;

        // only these two vertices matter, the rest are irrelevant
        m_index[0] = (axes[0] > 0.0f ? 0 : 1) +
                     (axes[1] > 0.0f ? 0 : 2) +
                     (axes[2] > 0.0f ? 0 : 4);
        // by ordering the vertices this way, opposite corners add up to 7
        m_index[7] = 7 - m_index[0];
    }

    // Find projections onto line
    float origin = axis.dot(obb->x);
    m_min      = origin - maxProjectedExtent;
    m_max      = origin + maxProjectedExtent;
}

