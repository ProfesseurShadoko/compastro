
#include "force_engine.hpp"
#include <iostream>
#include "tree.hpp"

std::vector<Eigen::Vector3d> ForceEngine::computeForce(Method method) const {
    switch (method) {
        case Method::direct:
            return directForce();
        
        case Method::tree:
            return treeForce();

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
            forces[j] -= f;
        }
    }
    return forces;
}


std::vector<Eigen::Vector3d> ForceEngine::treeForce() const {
    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));
    Octree tree(1);
    tree.insert(particles);

    for (int i = 0; i < particles.size(); i++) {
        forces[i] = tree.getForce(particles.get(i));
    }
    
    return forces;
}