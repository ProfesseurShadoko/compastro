

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
    std::cout << "Force call count (direct): " << Particle::popForceCallCounter() << std::endl;

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

    std::cout << "Force call count (monopole):" << Particle::popForceCallCounter() << std::endl;
}

void evolveMonopoleOnData(ParticleSet& particles) {
    // let's evolve evryone a thousand times

    Message("Evolving with Monopole, hang on tight!", "?");
    ParticleSet particles_over_time = ParticleSet(particles); // snapshot of the particles
    ForceEngine engine(particles);
    double dt = engine.crossingTime() / 100;
    int N_iter = 100;

    for (int i=0; i<N_iter; i++) {
        std::vector<Eigen::Vector3d> forces = engine.computeForce(Method::tree_mono);
        particles.applyForces(forces);
        particles.update(dt, IntegrationMethod::euler);
        particles_over_time.add(particles);
    }

    ParticleSet::save("files/output/evolve_monopole/mono_" + std::to_string(N_iter) + "-iter.csv", particles_over_time);
    
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
    std::cout << "Force call count (quadrupole):" << Particle::popForceCallCounter() << std::endl;

}





int main() {
    ParticleSet particles = loadData();
    //computeExactForcesOnData(particles);
    //treeMonopoleOnData(particles);
    //treeQuadOnData(particles);


    evolveMonopoleOnData(particles);

    Message("Execution complete!", "#");
    return 0;
}
