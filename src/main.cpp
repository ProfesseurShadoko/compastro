

#include <iostream>
#include <Eigen/Dense>
#include <string>

#include "particle.hpp"
#include "force_engine.hpp"
#include "tree.hpp"
#include "message.hpp"


// ### LOAD PARTICLES ###
ParticleSet& loadData() {
    Message("Computing the exact force from `data.txt` using direct summation.");
    static ParticleSet particles = ParticleSet::load("files/data.txt");
    Message::print("Loaded " + std::to_string(particles.size()) + " particles.");
    Message::print("Radius of the particle set: " + std::to_string(particles.radius())); // this is quite big compared to the distribution saddely. // the tree will get bigger.
    Message::print("First Particle:");
    particles.get(0).display();
    return particles;
}




/**
 * ########################
 * ### DIRECT SUMMATION ###
 * ########################
 */
void testDirect3BodyStable() {
    Message("Evolving 3 body problem with stable initial conditions using direct summation.");
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

    Message::print("Done.\n");
}


void computeExactForcesOnData(ParticleSet& particles) {
   
    ForceEngine engine(particles);
    Timer timer("Direct Solver");
    timer.start();
    std::vector<Eigen::Vector3d> forces = engine.computeForce(Method::direct);
    timer.stop();
    timer.display();

    ParticleSet::saveForces("files/output/DIRECT_FORCES_eps=" + std::to_string(engine.softening) + "_data.csv", forces);
}



/**
 * #####################
 * ### TREE MONOPOLE ###
 * #####################
 */

void treeMonopoleOnData(ParticleSet& particles) {
    Message("Computing the force from `data.txt` using the tree code with monopole approximation.");
    ForceEngine engine(particles);
    Timer timer("Tree Solver (Monopole)");
    timer.start();
    std::vector<Eigen::Vector3d> forces_tree_mono = engine.computeForce(Method::tree_mono);
    timer.stop();
    timer.display();
    ParticleSet::saveForces("files/output/TREE-MONO_FORCES_data.csv", forces_tree_mono);
}

void testDMQ() {
    ParticleSet particles = ParticleSet::load("files/tests/TEST_PARTICLES.txt");
    std::cout << "First Particle: " << std::endl;
    particles.get(0).display();
    ForceEngine engine(particles);

    std::vector<Eigen::Vector3d> forces_direct = engine.computeForce(Method::direct);
    std::vector<Eigen::Vector3d> forces_mono = engine.computeForce(Method::tree_mono);
    std::vector<Eigen::Vector3d> forces_quad = engine.computeForce(Method::tree_quad);

    ParticleSet::saveForces("files/tests/test_particles_direct.csv", forces_direct);
    ParticleSet::saveForces("files/tests/test_particles_mono.csv", forces_mono);
    ParticleSet::saveForces("files/tests/test_particles_quad.csv", forces_quad);
}


void testDMQ2() {
    Particle particle = {0.9, 0.1, 0.1};
    Particle pole1 = {-0.1, 0.12000000000000001, 0.1};
    Particle pole2 = {-0.1, 0.08, 0.1};

    ParticleSet particles = {particle, pole1, pole2};
    ForceEngine engine(particles);

    std::vector<Eigen::Vector3d> forces_direct = engine.computeForce(Method::direct);
    std::vector<Eigen::Vector3d> forces_mono = engine.computeForce(Method::tree_mono);
    std::vector<Eigen::Vector3d> forces_quad = engine.computeForce(Method::tree_quad);

    Eigen::Vector3d force_direct = forces_direct[0];
    Eigen::Vector3d force_mono = forces_mono[0];
    Eigen::Vector3d force_quad = forces_quad[0];

    std::cout << "Direct: " << force_direct.transpose() << std::endl;
    std::cout << "Mono: " << force_mono.transpose() << std::endl;
    std::cout << "Quad: " << force_quad.transpose() << std::endl;
}




/**
 * #######################
 * ### TREE QUADRUPOLE ###
 * #######################
 */

void treeQuadOnData(ParticleSet& particles) {
    Message("Computing the force from `data.txt` using the tree code with quadrupole approximation.");
    ForceEngine engine(particles);
    Timer timer("Tree Solver (Quadrupole)");
    timer.start();
    std::vector<Eigen::Vector3d> forces_tree_quad = engine.computeForce(Method::tree_quad);
    timer.stop();
    timer.display();
    ParticleSet::saveForces("files/output/TREE-QUAD_FORCES_data.csv", forces_tree_quad);
}





int main() {
    ParticleSet particles = loadData();
    //computeExactForcesOnData(particles);
    treeMonopoleOnData(particles);
    treeQuadOnData(particles);

    Message("Execution complete!", "#");
    return 0;
}
