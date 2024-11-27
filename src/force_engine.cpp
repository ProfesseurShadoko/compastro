

#include "force_engine.hpp"
#include <iostream>

std::vector<Eigen::Vector3d> ForceEngine::computeForce(Method method) const {
    switch (method) {
        case Method::direct:
            return directForce();
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
            Eigen::Vector3d r = particles.get(j).position - particles.get(i).position;
            double dist = r.norm();
            double f = particles.get(i).mass * particles.get(j).mass / (dist * dist + softening * softening);
            forces[j] -= f * r.normalized();
        }
    }
    return forces;
}