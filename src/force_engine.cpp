
#include "force_engine.hpp"
#include <iostream>
#include "tree.hpp"

double ForceEngine::openingAngle = 0.5;
double ForceEngine::softening = 0;


std::vector<Eigen::Vector3d> ForceEngine::computeForce(Method method) const {
    switch (method) {
        case Method::direct:
            return directForce();
        
        case Method::tree_mono:
            return treeForce(false); // TODO: change the way that this is implemented => if no quadrupole, dont compute the force using  0 filled matrix multiplication like we do

        case Method::tree_quad:
            return treeForce(true);

        default:
            throw std::runtime_error("Method not implemented.");
    }
}

std::vector<Eigen::Vector3d> ForceEngine::directForce() const {
    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));

    for (int i = 0; i < particles.size(); i++) {
        for (int j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue;
            }
            Eigen::Vector3d f = Particle::computeForce(particles.get(i), particles.get(j), softening);            
            forces[i] += f;
        }
    }
    return forces;
}


std::vector<Eigen::Vector3d> ForceEngine::treeForce(bool use_quad) const {

    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));
    Octree tree(particles.radius());
    tree.insert(particles);
    if (use_quad) {
        tree.computeQuadrupoles();
    }

    for (int i = 0; i < particles.size(); i++) {
        forces[i] = tree.getForce(particles.get(i), openingAngle);
    }
    
    return forces;
}