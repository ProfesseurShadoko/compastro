
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "particle.hpp"
#include <string>

enum class Method {
    direct,
    tree_mono,
    tree_quad,
    pm
};



class ForceEngine {
public:
    static double softening;
    static double openingAngle;
    
    ParticleSet& particles; // avoid copy here! allows in place modification!
    

private:
    /**
     * Direct computation of force between particles in O(n^2)
     */
    std::vector<Eigen::Vector3d> directForce() const;

    /**
     * Tree code computation of force between particles in O(n log n) (n particles that act on a tree of depth log n).
     * If not use_quadrupole, the quadrupole moment will be left to 0 and thus the force computation will only contain the monopole term.
     */
    std::vector<Eigen::Vector3d> treeForce(bool use_quadrupole) const;

public:
    ForceEngine(ParticleSet& particles) : particles(particles) {};

    /**
     * Compute the force between particles given parameters and chosen method
     */
    std::vector<Eigen::Vector3d> computeForce(Method method) const;

    double crossingTime() const;
};