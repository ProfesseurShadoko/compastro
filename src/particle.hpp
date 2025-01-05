
#pragma once
#include <Eigen/Dense>
#include <vector>
#include <initializer_list>




class Particle {

    public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double mass;
    double current_time = 0;
    double potentialEnergy = 0;
    

private:
    /**
     * @brief Counts how many times the function computeForce has been called.
     * 
     * This will provide insights into how much force calculation gains from the tree code.
     * Indeed, the monopole is computed with computeForce(particle, part_com).
     */
    static long long forceCallCounter;
    static long long id_counter;
    int id;
    

public:
    Particle();
    Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass);

    /**
     * Copy constructor. The id of the particle also gets copied!
     */
    Particle(const Particle& p);

    /**
     * Initialize particle by giving position only. Other quantities are set to zero (1 for mass).
     */
    Particle(std::initializer_list<double> init);

    /**
     * Compares id of particles. Usefull in the tree code. Careful! Does not compare positions. And equality survives copying:
     * particle == Particle(particle) is true.
     */
    bool operator==(const Particle& p) const {
        return getId() == p.getId();
    }

    int getId() const;


    /**
     * Increments current time by dt.
     */
    void updateCurrentTime(double dt);
    void display();

    /**
     * @brief Compute the force between two particles. This will return the force applied by the attractor on the particle!
     */
    static Eigen::Vector3d computeForce(Particle& particle, Particle& p_attractor, double eps = 0);

    /**
     * @brief Compute the potential applied by the attracotr on the particle. Returns normalized potential (for actual potential, multiply by particle mass!).
     */
    static double computePotential(Particle& particle, Particle& p_attractor, double eps = 0);

    /**
     * @brief Counts how many times the function computeForce has been called. Reset the counter.
     * 
     * This will provide insights into how much force calculation gains from the tree code.
     * Indeed, the monopole is computed with computeForce(particle, part_com).
     */
    static long long popForceCallCounter();

    /**
     * @brief Counts how many times the function computeForce has been called. In order to reset counter, call popForceCallCounter.
     * 
     * This will provide insights into how much force calculation gains from the tree code.
     * Indeed, the monopole is computed with computeForce(particle, part_com).
     */
    static long long getForceCallCounter();
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
     * FYI: this is a shallow copy, so the particles are not copied!
     */
    ParticleSet(std::vector<Particle> particles);

    ParticleSet(std::initializer_list<Particle> init);

    /**
     * Copy constructor. Actually is a deep copy, by magic of std::vector.
     */
    ParticleSet(const ParticleSet& ps);


    /**
     * Increments current time of each particle by dt.
     */
    void updateCurrentTime(double dt);

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
     * Split the set into all distinct time steps. Assumes particles are sorted in time from the beginning.
     */
    std::vector<ParticleSet> split();


    /**
     * Slice. If start greater than end (eg if end = -1), all particles are added.
     */
    ParticleSet slice(int start, int end);

    /**
     * Take first N elements.
     */
    ParticleSet slice(int end) {
        return slice(0, end);
    }

    /**
     * @brief Get the number of particles in the set.
     */
    int size();

    /**
     * @brief assuming that the particles are centered around the origin. This function returns the radius of the cube => max_{part}(|x|, |y|, |z|)
     */
    double radius();

    /**
     * @brief changes velocities and positions to put everyone into center of mass frame, centered on 0
     */
    void com();

    /**
     * @brief returns sum of kinetic energy and potential energy of the system
     */
    double totalEnergy();

    /**
     * Initialize the set with a number of particles, randomly distributed in a cube.
     */
    static ParticleSet load(std::string filename);

    /**
     * @brief Export the particle set to a file. The resulting file is a .csv file.
     */
    static void save(std::string filename, ParticleSet ps);

    /**
     * Exports the forces acting on the particles to a file. The indeces will match the ones of the particles! Resulting file is a .csv file.
     */
    static void saveForces(std::string filename, std::vector<Eigen::Vector3d> forces);


};