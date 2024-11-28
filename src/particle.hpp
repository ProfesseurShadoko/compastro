
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <initializer_list>

enum class IntegrationMethod {
    euler,
    leapfrog,
};


class Particle {

    public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acceleration;
    double mass;

    Particle();
    Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass);
    Particle(const Particle& p);

    /**
     * Initialize particle by giving position only. Other quantities are set to zero (1 for mass).
     */
    Particle(std::initializer_list<double> init);

    /**
     * @brief Apply a force to the particle. Updates the acceleration, with the formula F = ma.
     */
    void applyForce(Eigen::Vector3d force);

    /**
     * @brief Reset the force acting on the particle. This function needs to be called at each time step.
     */
    void resetForce();

    /**
     * @brief Update the particle's position and velocity based on the current acceleration.
     */
    void update(double dt, IntegrationMethod method = IntegrationMethod::euler);
    void display();
};


/**
 * @brief A set of particles. This class is a wrapper around a vector of particles.
 * The idea is to always use the method get(int i) to access a particle without copy.ADJ_OFFSET_SINGLESHOT
 * The set is initialized in the beginning, and after that no copy shall be made!
 * 
 * With this set, each particle has an index, which is good for tree code and all.
 */
class ParticleSet {
    public:
    std::vector<Particle> particles; // this will be the basic element, so we instantiate it here once, and it never gets copied after

    ParticleSet();

    /**
     * Here a copy is made, but its the last time! The particles vector should never be copied again.
     */
    ParticleSet(std::vector<Particle> particles);

    ParticleSet(std::initializer_list<Particle> init);

    /**
     * Calls the applyForce method for each particle in the set.
     */
    void applyForces(std::vector<Eigen::Vector3d> forces);

    /**
     * Calls the resetForce method for each particle in the set.
     */
    void resetForces();

    /**
     * Calls the update method for each particle in the set.
     */
    void update(double dt, IntegrationMethod method = IntegrationMethod::euler);

    /**
     * @brief Add a particle to the set. Particle is copied! This is important when you store a trajectory.
     */
    void add(Particle p);



    /**
     * @brief Add a particle set to the set. Particles are copied! This is important when you store a trajectory.
     */
    void add(ParticleSet ps);

    /**
     * @brief Add a vector of particles to the set. Particles are copied! This is important when you store a trajectory.
     */
    void add(std::vector<Particle> ps);

    /**
     * @brief Add particles to the set. Particles are copied! This is important when you store a trajectory.
     */
    void add(std::initializer_list<Particle> init);

    /**
     * Access a particle by index. Allows in place modification.
     */
    Particle& get(int i);

    /**
     * @brief Get the number of particles in the set.
     */
    int size();

    /**
     * Initialize the set with a number of particles, randomly distributed in a cube.
     */
    static ParticleSet load(std::string filename);

    /**
     * @brief Export the particle set to a file. Follows the same format as the loadFile method.
     */
    static void save(std::string filename, ParticleSet ps);

    /**
     * Exports the forces acting on the particles to a file. The indeces will match the ones of the particles!
     */
    static void saveForces(std::string filename, std::vector<Eigen::Vector3d> forces);
};