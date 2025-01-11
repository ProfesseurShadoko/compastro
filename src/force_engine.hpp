
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "particle.hpp"
#include <string>

enum class Method {
    direct,
    direct_opt,
    tree_mono,
    tree_quad,
    pm
};

enum class IntegrationMethod {
    euler,
    leapfrog,
    rk2,
    rk4,
    symplectic, // semi implicit euler => conservs energy
    verlet,
};



class ForceEngine {
public:
    static double softening;
    static double openingAngle;
    static bool compute_potential;
    static std::vector<double> totalEnergies; // when 'evolve' is called this gets filled with the total energy of the system at each saved step
    ParticleSet& particles; // avoid copy here! allows in place modification!
    

private:
    double fixed_radius = -1;
    /**
     * Direct computation of force between particles in O(n^2)
     */
    std::vector<Eigen::Vector3d> directForce();

    /**
     * Direct computation of potential energy of the system in O(n^2)
     */
    std::vector<double> directPotential();

    /**
     * Direct computation of force between particles in O(n(n+1)/2) with optimized version
     */
    std::vector<Eigen::Vector3d> directForceOpt();

    /**
     * Direct computation of potential energy of the system in O(n(n+1)/2) with optimized version
     */
    std::vector<double> directPotentialOpt();

    /**
     * Tree code computation of force between particles in O(n log n) (n particles that act on a tree of depth log n).
     * If not use_quadrupole, the quadrupole moment will be left to 0 and thus the force computation will only contain the monopole term.
     */
    std::vector<Eigen::Vector3d> treeForce(bool use_quadrupole);

    std::vector<double> treePotential(bool use_quadrupole);


    void euler(double dt, Method method);

    void leapfrog(double dt, Method method);

    void rk2(double dt, Method method);

    void rk4(double dt, Method method);

    void symplectic(double dt, Method method);


public:
    ForceEngine(ParticleSet& particles) : particles(particles) {};

    
    /**
     * Compute the force between particles given chosen method
     */
    std::vector<Eigen::Vector3d> computeForce(Method method);

    /**
     * Compute the potential energy of the system given chosen method
     */
    std::vector<double> computePotential(Method method);

    /**
     * Apply forces to particles, given chosen method. Does not update particle current time though!
     */
    void updateParticles(double dt, std::vector<Eigen::Vector3d> forces, IntegrationMethod i_method);

    /**
     * 1. Compute the force between particles, given chosen method
     * 2. Updates particle velocities and positions inplace according to dt, given chosen i_method
     * 
     * Basically calls successively euler (or an other method depending on i_method) who calls computeForce(method). And updates particle current time.
     * 
     */
    void evolve(double dt, Method method, IntegrationMethod i_method);

    /**
     * 1. Compute the force between particles, given chosen method
     * 2. Updates particle velocities and positions inplace according to dt, given chosen i_method
     * 3. Copies the current set of particles (or rather the first N_save particles) in the ParticleSet that will be returned in the end
     * 4. Loops N_iter times
     * 
     * Does not return the current state of particles, but a copy of N_save particles at each of the N_iter step. Time information is stored in the current_time of stored particles.
     * Steps 1. and 2. are done in the evolve method (overloaded).
     * If N_skip is 4, then only the 4th, 8th, 12th, ... particle set will be stored. Defaults to 1.
     */
    ParticleSet evolve(double dt, Method method, IntegrationMethod i_method, int N_iter, int N_save = -1, int N_skip = 1, bool cap_radius = false);

    double crossingTime() const;
};