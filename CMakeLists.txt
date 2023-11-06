cmake_minimum_required(VERSION 3.12)
project(rigidBodyTutorial)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)

set(discregrid_SOURCE 3rdParty/Discregrid/src/acceleration/bounding_sphere_hierarchy.cpp
		3rdParty/Discregrid/src/cubic_lagrange_discrete_grid.cpp
		3rdParty/Discregrid/src/discrete_grid.cpp
		3rdParty/Discregrid/src/geometry/mesh_distance.cpp
		3rdParty/Discregrid/src/geometry/point_triangle_distance.cpp
		3rdParty/Discregrid/src/mesh/entity_containers.cpp
		3rdParty/Discregrid/src/mesh/entity_iterators.cpp
		3rdParty/Discregrid/src/mesh/triangle_mesh.cpp
		3rdParty/Discregrid/src/utility/timing.cpp )
set(discregrid_HEADERS 3rdParty/Discregrid/include/Discregrid/acceleration/bounding_sphere.hpp
		3rdParty/Discregrid/include/Discregrid/acceleration/bounding_sphere_hierarchy.hpp
		3rdParty/Discregrid/include/Discregrid/acceleration/kd_tree.hpp
		3rdParty/Discregrid/include/Discregrid/acceleration/kd_tree.inl
		3rdParty/Discregrid/include/Discregrid/cubic_lagrange_discrete_grid.hpp
		3rdParty/Discregrid/include/Discregrid/discrete_grid.hpp
		3rdParty/Discregrid/include/Discregrid/geometry/mesh_distance.hpp
		3rdParty/Discregrid/include/Discregrid/mesh/entity_containers.hpp
		3rdParty/Discregrid/include/Discregrid/mesh/entity_iterators.hpp
		3rdParty/Discregrid/include/Discregrid/mesh/halfedge.hpp
		3rdParty/Discregrid/include/Discregrid/mesh/triangle_mesh.hpp
		3rdParty/Discregrid/include/Discregrid/utility/lru_cache.hpp
		3rdParty/Discregrid/include/Discregrid/utility/serialize.hpp
		3rdParty/Discregrid/src/data/z_sort_table.hpp
		3rdParty/Discregrid/src/geometry/point_triangle_distance.hpp
		3rdParty/Discregrid/src/utility/spinlock.hpp
		3rdParty/Discregrid/src/utility/timing.hpp )

set(collision_HEADERS include/collision/CollisionDetect.h 
                      include/collision/Geometry.h )
set(collision_SOURCE include/collision/CollisionDetect.cpp)

set(rigidbody_HEADERS include/rigidbody/RigidBody.h
                      include/rigidbody/RigidBodySystem.h 
                      include/rigidbody/Scenarios.h )
set(rigidbody_SOURCE include/rigidbody/RigidBody.cpp
		include/rigidbody/RigidBodySystem.cpp)

set(contact_HEADERS include/contact/Contact.h )
set(contact_SOURCE include/contact/Contact.cpp)

set(solver_HEADERS include/solvers/Solver.h 
                   include/solvers/SolverBoxPGS.h )
set(solver_SOURCE src/solvers/SolverBoxPGS.cpp )

set(util_HEADERS include/util/OBJLoader.h
                 include/util/MeshAssets.h )
set(util_SOURCE include/util/OBJLoader.cpp
		include/util/MeshAssets.cpp)

set(viewer_HEADERS include/viewer/SimViewer.h)
set(viewer_SOURCE include/viewer/SimViewer.cpp)
        
add_executable (rigidBodyTutorial main.cpp
          ${contact_HEADERS} ${contact_SOURCE} 
          ${collision_HEADERS} ${collision_SOURCE}
          ${rigidbody_HEADERS} ${rigidbody_SOURCE}
          ${solver_HEADERS} ${solver_SOURCE}
		  ${util_HEADERS} ${util_SOURCE}
          ${viewer_HEADERS} ${viewer_SOURCE}
		  ${discregrid_HEADERS} ${discregrid_SOURCE})

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/3rdParty/polyscope)
target_link_libraries(rigidBodyTutorial OpenGL::GL polyscope)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 3rdParty/Eigen3/include/eigen3 3rdParty/Discregrid/include ${COMMON_INCLUDES})

# Setup source code groups
#
source_group(src FILES main.cpp MainWindow.cpp)
source_group(include FILES MainWindow.h)
source_group(ui FILES main.ui)
source_group(3rdParty/discregrid/src FILES ${discregrid_SOURCE})
source_group(3rdParty/discregrid/include FILES ${discregrid_HEADERS})
source_group(src/collision FILES ${collision_SOURCE})
source_group(include/collision FILES ${collision_HEADERS})
source_group(src/contact FILES ${contact_SOURCE})
source_group(include/contact FILES ${contact_HEADERS})
source_group(src/rigidbody FILES ${rigidbody_SOURCE})
source_group(include/rigidbody FILES ${rigidbody_HEADERS})
source_group(src/solvers FILES ${solver_SOURCE})
source_group(include/solvers FILES ${solver_HEADERS})
source_group(src/util FILES ${util_SOURCE})
source_group(include/util FILES ${util_HEADERS})
source_group(src/viewer FILES ${viewer_SOURCE})
source_group(include/viewer FILES ${viewer_HEADERS})
source_group(src/util FILES ${util_SOURCE})

if(MSVC) 
	set_property(TARGET rigidBodyTutorial PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()