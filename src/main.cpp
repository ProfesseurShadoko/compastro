

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

    Timer timer2("Direct Solver (Potential)");
    timer2.start();
    engine.computePotential(Method::direct);
    timer2.stop();
    timer2.display();

    ParticleSet::saveForces("files/output/DIRECT_FORCES_eps=" + std::to_string(engine.softening) + "_data.csv", forces);
     ParticleSet::save("files/output/DIRECT_POTENTIAL_data.csv", particles);
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

    Timer timer2("Direct Solver (Potential)");
    timer2.start();
    engine.computePotential(Method::tree_mono);
    timer2.stop();
    timer2.display();


    ParticleSet::saveForces("files/output/TREE-MONO_FORCES_data.csv", forces_tree_mono);
    ParticleSet::save("files/output/TREE-MONO_POTENTIAL_data.csv", particles);

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

    Timer timer2("Direct Solver (Potential)");
    timer2.start();
    engine.computePotential(Method::tree_quad);
    timer2.stop();
    timer2.display();

    ParticleSet::saveForces("files/output/TREE-QUAD_FORCES_data.csv", forces_tree_quad);
    std::cout << "Force call count (quadrupole):" << Particle::popForceCallCounter() << std::endl;

    ParticleSet::save("files/output/TREE-QUAD_POTENTIAL_data.csv", particles);

}


void testIntegrationMethods() {
    Particle sun = Particle(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    double dt = 1e-3;
    ForceEngine::softening = 0;

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

    std::cout << "System RK4:" << std::endl;
    rk4_over_time.get(0).display();
}

void simulateMilkyWay() {
    ParticleSet particles = ParticleSet::load("files/tests/galaxy/milkyway.txt");
    ForceEngine engine(particles);

    double softening = 0.1; // 1% of R
    double crossing_time = 0.1;
    double dt = crossing_time / 50;
    int N_iter = 3000;
    int N_save = -1; // -1 means all of them
    int N_skip = 5;
    double initial_radius = particles.radius();

    ForceEngine::softening = softening;
    ForceEngine::openingAngle = 0.5;

    ForceEngine engine_initial(particles);

    
    ParticleSet over_time = engine.evolve(dt, Method::tree_mono, IntegrationMethod::symplectic, N_iter, N_save,N_skip);

    ParticleSet::save("files/tests/galaxy/milkyway_over_time.csv", over_time);

    Window window(over_time, initial_radius);
    window.animate();
}


void testGraphics() {
    ParticleSet particles = ParticleSet::load("files/tests/galaxy/milkyway.txt");
    particles = particles.slice(1000);
    ForceEngine engine(particles);

    ForceEngine::softening = 0.03;
    double dt = engine.crossingTime() / 50;
    int N_iter = 8000;
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

    //ForceEngine::softening = 0.01;
    double dt = engine.crossingTime() / 50;
    int N_iter = 100;
    int N_save = 10000;
    int N_skip = 1;
    double radius = 1;

    ParticleSet particles_over_time = engine.evolve(dt, Method::tree_quad, IntegrationMethod::rk4, N_iter, N_save, N_skip, true);
    ParticleSet::save("files/output/data_over_time.csv", particles_over_time);
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

void computeTimeComplexity() {

    ParticleSet particles = ParticleSet::load("files/data.txt");
    int N_max = 50000;
    int N_min = 10;
    int N_steps = 50;
    int N_repeat = 10;
    Timer timer;

    std::vector<int> Ns;
    // log scale since a complexity plot is in log log
    for (int i=0; i<N_steps; i++) {
        for (int _=0; _<N_repeat; _++) {
            Ns.push_back(N_min * pow((N_max / N_min), (double)i/(N_steps - 1)));
        }
    }
    Message("Data sizes initialized");

    // --- Direct ---
    Message("Direct Force Calculation");
    ProgressBar bar_direct = ProgressBar(Ns.size());
    std::vector<long> times_direct;
    std::vector<long> iters_direct;

    for (size_t n=0; n<Ns.size(); n++) {
        ParticleSet particles_sub = particles.slice(Ns[n]);
        ForceEngine engine = ForceEngine(particles_sub);
        timer.start();
        engine.computeForce(Method::direct);
        times_direct.push_back(timer.stop());
        iters_direct.push_back(Particle::popForceCallCounter());
        bar_direct.update();
    }
    Message("Direct Force Calculation Complete", "#");

    // --- Tree Monopole ---
    Message("Tree Monopole Force Calculation");
    ProgressBar bar_tree_mono = ProgressBar(Ns.size());
    std::vector<long> times_tree_mono;
    std::vector<long> iters_tree_mono;

    for (size_t n=0; n<Ns.size(); n++) {
        ParticleSet particles_sub = particles.slice(Ns[n]);
        ForceEngine engine = ForceEngine(particles_sub);
        timer.start();
        engine.computeForce(Method::tree_mono);
        times_tree_mono.push_back(timer.stop());
        iters_tree_mono.push_back(Particle::popForceCallCounter());
        bar_tree_mono.update();
    }

    // --- Tree Quad ---
    Message("Tree Quad Force Calculation");
    ProgressBar bar_tree_quad = ProgressBar(Ns.size());
    std::vector<long> times_tree_quad;
    std::vector<long> iters_tree_quad;

    for (size_t n=0; n<Ns.size(); n++) {
        ParticleSet particles_sub = particles.slice(Ns[n]);
        ForceEngine engine = ForceEngine(particles_sub);
        timer.start();
        engine.computeForce(Method::tree_quad);
        times_tree_quad.push_back(timer.stop());
        iters_tree_quad.push_back(Particle::popForceCallCounter());
        bar_tree_quad.update();
    }

    // --- Just Tree Construction --- //
    Message("Tree Construction");
    ProgressBar bar_tree_construction = ProgressBar(Ns.size());
    std::vector<long> times_tree_construction;
    std::vector<long> times_second_pass;

    for (size_t n=0; n<Ns.size(); n++) {
        ParticleSet particles_sub = particles.slice(Ns[n]);
        Octree tree = Octree(particles_sub.radius());
        timer.start();
        tree.insert(particles_sub);
        times_tree_construction.push_back(timer.stop());

        timer.start();
        tree.computeQuadrupoles();
        times_second_pass.push_back(timer.stop());
        bar_tree_construction.update();
    }

    // --- Save the results ---
    // into csv with columns N, direct, tree_mono, tree_quad, tree_construction, second_pass
    std::ofstream file("files/output/time_complexity.csv");
    file << "N,direct,tree_mono,tree_quad,tree_construction,second_pass,iter_direct,iter_mono,iter_quad" << std::endl;
    for (size_t n=0; n<Ns.size(); n++) {
        file << Ns[n] << "," << times_direct[n] << "," << times_tree_mono[n] << "," << times_tree_quad[n] << "," << times_tree_construction[n] << "," << times_second_pass[n] << "," << iters_direct[n] << "," << iters_tree_mono[n] << "," << iters_tree_quad[n] << std::endl;
    }
    file.close();
    Message("Time Complexity Computation Complete", "#");
    
}


void testPotentials() {
    Particle sun = Particle(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    double dt = 1e-3;
    ForceEngine::softening = 0;

    double e = 0.5; // eccentricity // in plane x-y
    Particle earth = Particle(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, sqrt(1 + e), 0), 1e-30); // we want the sun to stay still
    // initial velocity is sqrt(1 + e)

    int N_iter = 53315;
    ParticleSet system = {earth, sun};

    ForceEngine engine(system);
    engine.computePotential(Method::direct);
    engine.particles.get(0).display();

    // --- RK4 ---
    ParticleSet rk4_system = system; // deep copy
    ForceEngine rk4_engine = ForceEngine(rk4_system);
    ParticleSet rk4_over_time = rk4_engine.evolve(dt, Method::direct, IntegrationMethod::rk4, N_iter, 2); // we always use direct method here
    rk4_over_time.get(0).display();
    
    
};

void testEnergyConservation() {
    //evolveData();
    Particle p1 = Particle(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), 1);
    Particle p2 = Particle(Eigen::Vector3d(2, 0, 0), Eigen::Vector3d(0, -1, 0), 1);
    Particle p3 = Particle(Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(1, 0, 0), 1);
    Particle p4 = Particle(Eigen::Vector3d(4, 3, 0), Eigen::Vector3d(1, -1, 0), 1);
    ParticleSet particles = {p1, p2, p3, p4};
    particles.com();

    ForceEngine engine(particles);
    engine.computePotential(Method::direct);
    double energy = particles.totalEnergy();
    std::cout << "Total Energy: " << energy << std::endl;

    ParticleSet particles_over_time = engine.evolve(0.005, Method::direct, IntegrationMethod::rk4, 10000, 4, 20, true);
    Window window(particles_over_time, 10);
    window.animate();

    std::cout << "Total Energy: " << particles.totalEnergy() << std::endl;
    ParticleSet::save("files/output/direct_energy_conservation.csv", particles_over_time);
}


void computeForcesVariousOpeningAngles() {
    ParticleSet particles = ParticleSet::load("files/data.txt"); // 50000 of them, a lot
    //double R = 1;
    int N = particles.size();
    //double m = particles.get(0).mass;
    //double M_tot = m * N; // we assume all have same mass, which is the case
    //double G = 1;

    // let's create our array of opening angles
    double min_theta = 0.01;
    double max_theta = 1;
    int n_theta = 30;
    double theta;
    std::vector<double> thetas;
 
    std::cout << "Min theta: " << min_theta << std::endl;
    std::cout << "Max theta: " << max_theta << std::endl;
    if (max_theta < min_theta) throw("Not enough epsilon range");
    
    for (int i=0; i<n_theta; i++) {
        theta = min_theta + i*(max_theta - min_theta) / (n_theta-1);
        thetas.push_back(theta);
        //std::cout << "Current Epsilon: " << eps << std::endl;
    }

    // let's create the vectors to store the results
    std::vector<double> relative_errors_mono;
    std::vector<double> time_mono;
    std::vector<double> iters_mono;

    std::vector<double> relative_errors_quad;
    std::vector<double> time_quad;
    std::vector<double> iters_quad;

    

    // let's compute the forces for each one of them
    ForceEngine engine = ForceEngine(particles);
    Timer timer;
    int N_rep = 1;

    // let's compute reference with direct first
    Message("Computing the reference forces with direct summation");
    std::vector<Eigen::Vector3d> forces_direct;
    timer.start();
    for (int j=0; j<N_rep; j++) {
        forces_direct = engine.computeForce(Method::direct);
    }

    double time_direct = timer.stop() / N_rep;
    std::cout << "t_direct = " << time_direct << std::endl;
    std::cout << "f_count_direct_per_part = " << Particle::popForceCallCounter() / N << std::endl;

    ProgressBar bar = ProgressBar(n_theta * 2);
    for (size_t i=0; i<thetas.size(); i++) {
        bar.print("Computing forces for theta = " + std::to_string(thetas[i]));
        ForceEngine::openingAngle = thetas[i];

        // --- Monopole ---
        std::vector<Eigen::Vector3d> forces_mono;
        timer.start();
        for (int j=0; j<N_rep; j++) forces_mono = engine.computeForce(Method::tree_mono);
        time_mono.push_back(timer.stop() / N_rep);
        iters_mono.push_back((double)Particle::popForceCallCounter() / (N_rep * N));
        bar.update();

        // --- Quadrupole ---
        std::vector<Eigen::Vector3d> forces_quad;
        timer.start();
        for (int j=0; j<N_rep; j++) forces_quad = engine.computeForce(Method::tree_quad);
        time_quad.push_back(timer.stop() / N_rep);
        iters_quad.push_back((double)Particle::popForceCallCounter() / (N_rep * N));
        bar.update();

        // --- Compute the relative errors ---
        double error_mono = 0;
        double error_quad = 0;

        for (int j=0; j<particles.size(); j++) {
            error_mono += (forces_mono[j] - forces_direct[j]).norm() / forces_direct[j].norm();
            error_quad += (forces_quad[j] - forces_direct[j]).norm() / forces_direct[j].norm();
        }

        relative_errors_mono.push_back(error_mono / N);
        relative_errors_quad.push_back(error_quad / N);
    }


    // let's dispaly the results:
    std::cout << "theta_errm_errz_timem_timeq_iterm_iterq = [";
    for (int i=0; i<n_theta; i++) {
        theta = thetas[i];
        double error_mono = relative_errors_mono[i];
        double error_quad = relative_errors_quad[i];
        double time_m = time_mono[i];
        double time_q = time_quad[i];
        double iter_m = iters_mono[i];
        double iter_q = iters_quad[i];
        std::cout << "[" << theta << "," << error_mono << "," << error_quad << "," << time_m << "," << time_q << "," << iter_m << "," << iter_q << "],";
        
    }
    std::cout << "]" << std::endl;
    std::cout <<"Comment: force count is per particle." << std::endl;

}




int main() {
    //simulateMilkyWay();
    //evolveData();
    //testEnergyConservation();
    computeForcesVariousOpeningAngles();

    return 0;
}
