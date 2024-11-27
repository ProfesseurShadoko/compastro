
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "particle.hpp"
#include <string>

enum class Method {
    direct,
    tree,
    pm
};



class ForceEngine {
public:
    double softening = 0;
    ParticleSet& particles; // avoid copy here! allows in place modification!
    

private:
    /**
     * Direct computation of force between particles in O(n^2)
     */
    std::vector<Eigen::Vector3d> directForce() const;

public:
    ForceEngine(ParticleSet& particles) : particles(particles) {};

    /**
     * Compute the force between particles given parameters and chosen method
     */
    std::vector<Eigen::Vector3d> computeForce(Method method) const;
};