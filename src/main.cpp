

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


void testIntegrationMethods() {
    Particle sun = Particle(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    double dt = 1e-3;

    double e = 0.5; // eccentricity // in plane x-y
    Particle earth = Particle(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, sqrt(1 + e), 0), 1e-30); // we want the sun to stay still
    // initial velocity is sqrt(1 + e)

    int N_iter = 53315;
    ParticleSet system = {earth, sun};

    ProgressBar::mute();

    // --- Euler ---
    ParticleSet euler_system = system; // deep copy
    ForceEngine euler_engine = ForceEngine(euler_system);
    ParticleSet euler_over_time = euler_engine.evolve(dt, Method::direct, IntegrationMethod::euler, N_iter, 2); // we always use direct method here
    std::cout << "Euler: " << euler_over_time.size() << std::endl;
    
    // --- Leapfrog ---
    ParticleSet leapfrog_system = system; // deep copy
    ForceEngine leapfrog_engine = ForceEngine(leapfrog_system);
    ParticleSet leapfrog_over_time = leapfrog_engine.evolve(dt, Method::direct, IntegrationMethod::leapfrog, N_iter, 2); // we always use direct method here

    // --- RK2 ---
    ParticleSet rk2_system = system; // deep copy
    ForceEngine rk2_engine = ForceEngine(rk2_system);
    ParticleSet rk2_over_time = rk2_engine.evolve(dt, Method::direct, IntegrationMethod::rk2, N_iter, 2); // we always use direct method here

    // --- RK4 ---
    ParticleSet rk4_system = system; // deep copy
    ForceEngine rk4_engine = ForceEngine(rk4_system);
    ParticleSet rk4_over_time = rk4_engine.evolve(dt, Method::direct, IntegrationMethod::rk4, N_iter, 2); // we always use direct method here

    // --- Symplectic ---
    ParticleSet symplectic_system = system; // deep copy
    ForceEngine symplectic_engine = ForceEngine(symplectic_system);
    ParticleSet symplectic_over_time = symplectic_engine.evolve(dt, Method::direct, IntegrationMethod::symplectic, N_iter, 2); // we always use direct method here


    // Let's save everything
    ParticleSet::save("files/tests/integration_methods_notebook/euler.csv", euler_over_time);
    ParticleSet::save("files/tests/integration_methods_notebook/leapfrog.csv", leapfrog_over_time);
    ParticleSet::save("files/tests/integration_methods_notebook/rk2.csv", rk2_over_time);
    ParticleSet::save("files/tests/integration_methods_notebook/rk4.csv", rk4_over_time);
    ParticleSet::save("files/tests/integration_methods_notebook/symplectic.csv", symplectic_over_time);

    std::cout << "System:" << std::endl;
    system.get(0).display();
}




int main() {
    ParticleSet particles = loadData();

    computeExactForcesOnData(particles);
    treeMonopoleOnData(particles);
    treeQuadOnData(particles);
    return 0;
}
