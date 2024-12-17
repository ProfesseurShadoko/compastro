

#include <iostream>
#include <Eigen/Dense>
#include <string>

#include "particle.hpp"
#include "force_engine.hpp"
#include "tree.hpp"
#include "message.hpp"
#include "graphics.hpp"



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

void simulateMilkyWay() {
    ParticleSet particles = ParticleSet::load("files/tests/galaxy/milkyway.txt");
    ForceEngine engine(particles);

    double softening = 0.01; // 1% of R
    double crossing_time = 1.0;
    double dt = crossing_time / 50;
    int N_iter = 3000;
    int N_save = 1000; // -1 means all of them
    int N_skip = 10;

    ForceEngine::softening = softening;
    ForceEngine::openingAngle = 0.5;

    ForceEngine engine_initial(particles);

    
    ParticleSet over_time = engine.evolve(dt, Method::tree_mono, IntegrationMethod::symplectic, N_iter, N_save,N_skip);

    ParticleSet::save("files/tests/galaxy/milkyway_over_time.csv", over_time);
}


void testGraphics() {
    ParticleSet particles = ParticleSet::load("files/tests/galaxy/milkyway.txt");
    ForceEngine engine(particles);

    ForceEngine::softening = 0.03;
    double dt = engine.crossingTime() / 50;
    int N_iter = 10000;
    int N_save = 3000;
    int N_skip = 10;
    double radius = particles.radius()/1.4; // need to fix it here, later particles get ejected

    ParticleSet particles_over_time = engine.evolve(dt, Method::tree_mono, IntegrationMethod::symplectic, N_iter, N_save, N_skip, true);
    
    Window window(particles_over_time, radius);
    window.animate();
}

void evolveData() {
    ParticleSet particles = ParticleSet::load("files/data.txt");
    ForceEngine engine(particles);

    ForceEngine::softening = 0.01;
    double dt = engine.crossingTime() / 50;
    int N_iter = 100;
    int N_save = 10000;
    int N_skip = 1;
    double radius = 1;

    ParticleSet particles_over_time = engine.evolve(dt, Method::tree_mono, IntegrationMethod::symplectic, N_iter, N_save, N_skip, true);
    Window window(particles_over_time, radius);
    window.animate();
}

void computeForcesVariousEpsilons() {
    ParticleSet particles = ParticleSet::load("files/data.txt"); // 50000 of them, a lot
    double R = 1;
    int N = particles.size();
    double m = particles.get(0).mass;
    double M_tot = m * N; // we assume all have same mass, which is the case
    //double G = 1;

    // let's create our array of epsilons
    double min_epsilon = R / N * 0.1;
    double max_epsilon = R / 100 * 10;
    double eps;
    int n_epsilon = 100;
    std::vector<double> epsilons;
 
    std::cout << "Min epsilon: " << min_epsilon << std::endl;
    std::cout << "Max epsilon: " << max_epsilon << std::endl;
    if (max_epsilon < min_epsilon) throw("Not enough epsilon range");
    

    for (int i=0; i<n_epsilon; i++) {
        //eps = min_epsilon + i*(max_epsilon - min_epsilon) / (n_epsilon-1);
        // lets do log norm rather
        eps = min_epsilon * pow((max_epsilon / min_epsilon), (double)i/(n_epsilon - 1));
        epsilons.push_back(eps);
        //std::cout << "Current Epsilon: " << eps << std::endl;
    }

    // let's compute the forces for each one of them
    Method method = Method::tree_quad;
    ProgressBar bar = ProgressBar(n_epsilon);
    ForceEngine engine = ForceEngine(particles);

    std::vector<double> errors;
    double error;
    double fr;
    double f_theoretical;
    double a=0.08;
    Eigen::Vector3d r;
    /*
    # theoretical curve
    def f(r, a=0.08):
        return -M*m / (r+a)**2
    */
    for (int i=0; i<n_epsilon; i++) {
        error = 0;
        ForceEngine::softening = epsilons[i];
        std::vector<Eigen::Vector3d> forces = engine.computeForce(method);

        // let's compute the sqaured error as 
        for (int j=0; j<particles.size(); j++) {
            r = particles.get(j).position;
            f_theoretical = -M_tot * m / pow(r.norm()+a, 2);
            fr = (forces[j](0) * r(0) + forces[j](1) * r(1) + forces[j](2) * r(2))/r.norm();
            error += pow(f_theoretical - fr, 2);
        }
        bar.update();
        errors.push_back(error);
    }

    // let's dispaly the results:
    std::cout << "[";
    for (int i=0; i<n_epsilon; i++) {
        eps = epsilons[i];
        error = errors[i];
        std::cout << "[" << eps << "," << error << "],";
    }
    std::cout << "]" << std::endl;

}


void testPotentialComputation() {
    Particle p1 = Particle(
        Eigen::Vector3d(0, 0, 0),
        Eigen::Vector3d(1, 1, 0),
        1.
    );

    Particle p2 = Particle(
        Eigen::Vector3d(1, 0, 0.5),
        Eigen::Vector3d(0, 0, 0),
        1.
    );

    Particle p3 = Particle(
        Eigen::Vector3d(-1, -0.5, -0.3),
        Eigen::Vector3d(0.5, 0.5, 0.2),
        2.
    );


    ParticleSet particles = {p1, p2, p3};
    particles.com();
    ForceEngine engine(particles);

    double dt = 0.001;
    ParticleSet particles_over_time = engine.evolve(dt, Method::direct, IntegrationMethod::rk4, 10000, -1, 10);
    ParticleSet::save("files/tests/3body_random.csv", particles_over_time);

    Window window(particles_over_time, 10);
    window.animate();
}






int main() {
    //ParticleSet particles = loadData();

    //computeExactForcesOnData(particles);
    //treeMonopoleOnData(particles);
    //treeQuadOnData(particles);
    //simulateMilkyWay();
    //testGraphics();
    //evolveData();
    //computeForcesVariousEpsilons();

    testPotentialComputation();

    return 0;
}
