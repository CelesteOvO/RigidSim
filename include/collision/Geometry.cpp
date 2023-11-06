
#include "Geometry.h"

void SDFGeometry::computeInertia() {
    if( mesh )
    {
        Eigen::Matrix3f A;
        float vol = 0.0f;
        Eigen::Vector3f diag = Eigen::Vector3f::Zero();
        Eigen::Vector3f offd = Eigen::Vector3f::Zero();
        for(unsigned int i = 0; i < mesh->nFaces(); ++i)
        {
            const std::array<unsigned int, 3>& face = mesh->face(i);
            A.col(0) = mesh->vertex(face[0]).cast<float>();
            A.col(1) = mesh->vertex(face[1]).cast<float>();
            A.col(2) = mesh->vertex(face[2]).cast<float>();
            const float det = A.determinant();
            vol += det;

            for(int j = 0; j < 3; ++j)
            {
                const int j1 = (j+1) % 3;
                const int j2 = (j+2) % 3;
                diag(j) += (A(j,0)*A(j,1) + A(j,1)*A(j,2) + A(j,2)*A(j,0) +
                            A(j,0)*A(j,0) + A(j,1)*A(j,1) + A(j,2)*A(j,2)  ) * det; // divide by 60.0f later
                offd(j) += (A(j1,0)*A(j2,1) + A(j1,1)*A(j2,2) + A(j1,2)*A(j2,0)  +
                            A(j1,0)*A(j2,2) + A(j1,1)*A(j2,0) + A(j1,2)*A(j2,1)  +
                            2.0f*A(j1,0)*A(j2,0) + 2.0f*A(j1,1)*A(j2,1) + 2.0f*A(j1,2)*A(j2,2) ) * det; // divide by 120.0f later
            }
        }

        diag /= vol * (60.0f / 6.0f);  // divide by total volume (vol/6) since density= 1/volume
        offd /= vol * (120.0f / 6.0f);
        I0(0,0) = diag.y()+diag.z(); I0(0,1) = -offd.z(); I0(0,2) = -offd.y();
        I0(1,0) = -offd.z(); I0(1,1) = diag.x()+diag.z(); I0(1,2) = -offd.x();
        I0(2,0) = -offd.y(); I0(2,1) = -offd.x(); I0(2,2) = diag.x()+diag.y();
    }
    else
    {
        I0.setZero();
    }
}

SDFGeometry::SDFGeometry(const std::string &filename, const std::array<unsigned int, 3> &resolution): mesh(), sdf() {
    Eigen::AlignedBox3d domain;
    domain.setEmpty();

    mesh = std::make_unique<Discregrid::TriangleMesh>(filename);

    Discregrid::MeshDistance md(*mesh);

    // Compute the bounding box of the mesh
    //
    for (auto const& x : mesh->vertices())
    {
        domain.extend(x);
    }

    // Enlarge the bounding box slightly
    //
    domain.max() += 1e-3f * domain.diagonal().norm() * Eigen::Vector3d::Ones();
    domain.min() -= 1e-3f * domain.diagonal().norm() * Eigen::Vector3d::Ones();

    sdf = std::make_unique<Discregrid::CubicLagrangeDiscreteGrid>(domain, resolution);

    // Create function that returns the SD
    auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
    func = [&md](Eigen::Vector3d const& xi) { return md.signedDistanceCached(xi); };
    sdf->addFunction(func, true);

    computeInertia();
}


