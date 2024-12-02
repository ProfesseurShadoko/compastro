

#include <iostream>
#include <Eigen/Dense>
#include <string>

#include "particle.hpp"
#include "force_engine.hpp"
#include "tree.hpp"
#include "message.hpp"


void testDirect3BodyStable() {
    ParticleSet particles = {
        Particle(Eigen::Vector3d(0.97, -0.243, 0), Eigen::Vector3d(0.4662, 0.4324, 0), 1), // these are values for stable 3 body problem
        Particle(Eigen::Vector3d(-0.97, 0.243, 0), Eigen::Vector3d(0.4662, 0.4324, 0), 1),
        Particle(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(-0.9324, -0.8647, 0), 1)
    };
    ParticleSet particles_over_time = ParticleSet(particles);

    // let's compute the forces and evolve our particles
    ForceEngine engine(particles);
    double dt = 0.01;

    for (int i = 0; i < 5000; i++) {
        std::vector<Eigen::Vector3d> forces = engine.computeForce(Method::direct);
        particles.applyForces(forces);
        particles.update(dt, IntegrationMethod::euler);
        particles_over_time.add(particles); // this makes a copy => takes a snapshot of the particles
    }

    //save the particles over time stored in the particle_set into a .csv file
    ParticleSet::save("files/tests/DIRECT_EULER_N=3_stable.csv", particles_over_time);
}



int main() {

    /**
     * ########################
     * ### DIRECT SUMMATION ###
     * ########################
     */
    /*
    Message("Evolving 3 body problem with stable initial conditions using direct summation.");
    // let's start with the stable 3 body problem. Let's see if it works!
    testDirect3BodyStable();
    Message::print("Done.\n");*/

    Message("Computing the exact force from `data.txt` using direct summation.");
    ParticleSet particles = ParticleSet::load("files/data.txt");
    Message::print("Loaded " + std::to_string(particles.size()) + " particles.");
    Message::print("Radius of the particle set: " + std::to_string(particles.radius())); // this is quite big compared to the distribution saddely. // the tree will get bigger.
    Message::print("First Particle:");
    particles.get(0).display();
    
    ForceEngine engine(particles);
    ForceEngine::softening = 0;
    ForceEngine::openingAngle = 0.5;
    Timer timer;

    /*
    timer = Timer("Direct Solver");
    timer.start();
    std::vector<Eigen::Vector3d> forces = engine.computeForce(Method::direct);
    timer.stop();
    timer.display();

    ParticleSet::saveForces("files/output/DIRECT_FORCES_eps=" + std::to_string(engine.softening) + "_data.csv", forces);
    */

    /**
     * #####################
     * ### TREE MONOPOLE ###
     * #####################
     */
    Message("Computing the force from `data.txt` using the tree code with monopole approximation.");
    timer = Timer("Tree Solver (Monopole)");
    timer.start();
    std::vector<Eigen::Vector3d> forces_tree_mono = engine.computeForce(Method::tree_mono);
    timer.stop();
    timer.display();
    ParticleSet::saveForces("files/output/TREE-MONO_FORCES_data.csv", forces_tree_mono);


    /**
     * #######################
     * ### TREE QUADRUPOLE ###
     * #######################
     */
    Message("Computing the force from `data.txt` using the tree code with quadrupole approximation.");
    timer = Timer("Tree Solver (Quadrupole)");
    timer.start();
    std::vector<Eigen::Vector3d> forces_tree_quad = engine.computeForce(Method::tree_quad);
    timer.stop();
    timer.display();
    ParticleSet::saveForces("files/output/TREE-QUAD_FORCES_data.csv", forces_tree_quad);

    Message("Execution complete!", "#");
    return 0;
}
