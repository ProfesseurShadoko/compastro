
#include "force_engine.hpp"
#include <iostream>
#include "tree.hpp"
#include "message.hpp"
#include <sstream> // scientific notation
#include <iomanip>

double ForceEngine::openingAngle = 0.5;
double ForceEngine::softening = 0.00113221;
bool ForceEngine::compute_potential = true;
std::vector<double> ForceEngine::totalEnergies;


/**
 * --------------------------
 * ! --- Public Methods --- !
 * --------------------------
 */

std::vector<Eigen::Vector3d> ForceEngine::computeForce(Method method) {

    switch (method) {
        case Method::direct:
            return directForce();
        
        case Method::tree_mono:
            return treeForce(false);

        case Method::tree_quad:
            return treeForce(true);
        
        case Method::direct_opt:
            return directForceOpt();

        default:
            throw std::runtime_error("Method not implemented.");
    }
}

std::vector<double> ForceEngine::computePotential(Method method) {
    if (!compute_potential) {
        Message("Potential computation is disabled, but called anyway!", "!");
    }
    //std::cout << "Computing potential" << std::endl;
    std::vector<double> potentials;
    
    switch (method) {
        case Method::direct:
            potentials = directPotential();
            break;
        
        case Method::direct_opt:
            potentials = directPotentialOpt();
            break;
        
        case Method::tree_mono:
            potentials =  treePotential(false);
            break;

        case Method::tree_quad:
            potentials = treePotential(true);
            break;

        default:
            throw std::runtime_error("Method not implemented.");
    }

    // let's assign the potentials to the particles
    double radius = (fixed_radius == -1) ? particles.radius() : fixed_radius;
    for (int i=0; i<particles.size(); i++) {

        if (method==Method::tree_mono || method==Method::tree_quad) {
            // if particle outside of  the radius of the tree, we don't want to update its potential since we didn't apply force on it => must keep constant velocity & constant energy
            if (particles.get(i).position.norm() > radius) continue;
        }
        
        
        particles.get(i).potentialEnergy = potentials[i];
    }
    //std::cout << "Potentials[0]: " << potentials[0] << std::endl;
    //std::cout << particles.get(0).potentialEnergy << std::endl;
    return potentials;
}


void ForceEngine::evolve(double dt, Method method, IntegrationMethod i_method) {
    particles.updateCurrentTime(dt);

    switch (i_method) {
        case IntegrationMethod::euler:
            euler(dt, method);
            break;
        
        case IntegrationMethod::leapfrog:
            leapfrog(dt,  method);
            break;
        
        case IntegrationMethod::rk2:
            rk2(dt, method);
            break;

        case IntegrationMethod::rk4:
            rk4(dt, method);
            break;

        case IntegrationMethod::symplectic:
            symplectic(dt, method);
            break;
        
        default:
            throw std::runtime_error("Integration Method not implemented.");
    }

    if (compute_potential) {
        computePotential(method); // no that we moved our particles, we need to update their potentials
        particles.totalEnergy(); // this assigns to each particle the total energy of the system in the envEnergy field
    }; 
}


ParticleSet ForceEngine::evolve(double dt, Method method, IntegrationMethod i_method, int N_iter, int N_save, int N_skip, bool cap_radius) {

    /**
     * -----------------------------
     * ! --- Print Information --- !
     * -----------------------------
     */
    Message("Evolving particles:");
    Message::print(" > Method: " + methodToString(method));
    Message::print(" > Integration Method: " + methodToString(i_method));
    Message::print(" > Particle Set Size: " + std::to_string(particles.size()));
    Message::print(" > Number of iterations: " + std::to_string(N_iter));
    Message::print(" > Number of particles to save: " + std::to_string(N_save));
    Message::print(" > Number of iterations to skip: " + std::to_string(N_skip));

    int  real_n_save = (N_save==-1) ? particles.size() : N_save;
    real_n_save = std::min(real_n_save, particles.size());
    std::stringstream ss;
    ss << std::scientific << std::setprecision(1) << ((double)N_iter * real_n_save / N_skip);

    Message::print(" > Total number of particles saved: " + cstr(ss.str()).yellow());

    std::string compute_str = (compute_potential ? cstr("True").green() : cstr("False").red());
    Message::print(" > Compute potentials: " + compute_str);
    Message::print(" > Opening Angle: " + std::to_string(openingAngle));
    Message::print(" > Softening: " + std::to_string(softening));
    Message::print(" > Crossing Time: " + std::to_string(crossingTime()));

    if (cap_radius) {
        fixed_radius = particles.radius() * 10;
        Message::print(" > Fixed Radius: " + std::to_string(fixed_radius));
    }

    if (compute_potential) {
        computePotential(method); // also assigns the potentials

        std::stringstream ss2;
        ss2 << std::scientific << std::setprecision(4) << particles.totalEnergy(); // also saves total energy <3
        Message::print(" > Initial Total Energy: " + cstr(ss2.str()).green());
    }
    
    
    //std::cout << "Initial particle" << std::endl;
    //particles.get(0).display();

    /**
     * ---------------------
     * ! --- Evolution --- !
     * ---------------------
     */

    ProgressBar bar(N_iter);
    ParticleSet particles_over_time = ParticleSet(particles).slice(0, real_n_save); // snapshot of the particles
    long long forceCallCounter = Particle::getForceCallCounter();
    Timer timer("Evolve function: (" + std::to_string((int) method) + ", " + std::to_string((int) i_method) + ")");

    timer.start();
    for (int i=0; i<N_iter; i++) {
        bar.update();
        evolve(dt, method, i_method); // potentials get assgned again here
        if (i % N_skip != 0) continue; // only one out of N_skip steps get saved
        
        particles_over_time.add(particles.slice(0, real_n_save));
    }
    timer.stop();

    forceCallCounter = Particle::getForceCallCounter() - forceCallCounter;
    Message::print(" > Force Call Counter: " + std::to_string(forceCallCounter));
    Message::print(" > Force Call Counter per Particle per Iteration: " + std::to_string((double) forceCallCounter / particles.size() / N_iter));
    Message::print(" > Final radius: " + std::to_string(particles.radius()));

    if (compute_potential) {
        std::stringstream ss3;
        ss3 << std::scientific << std::setprecision(4) << particles.totalEnergy();
        Message::print(" > Final Total Energy: " + cstr(ss3.str()).green());
    }
   

    timer.display();
    Message("Evolution complete!", "#");

    
    return particles_over_time;
}

double ForceEngine::crossingTime() const {
    double a = 0.0804;
    double R_half = (1 + sqrt(2)) * a;
    double M_tot = particles.get(0).mass * particles.size();
    return sqrt(
        pow(R_half, 3) / (M_tot)
    );
}


/**
 * -------------------------------------
 * ! --- Force Computation Methods --- !
 * -------------------------------------
 */

std::vector<Eigen::Vector3d> ForceEngine::directForce() {
    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));

    for (int i = 0; i < particles.size(); i++) {
        for (int j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue;
            }
            Eigen::Vector3d f = Particle::computeForce(particles.get(i), particles.get(j), softening);            
            forces[i] += f;
        }
    }
    return forces;
}

std::vector<double> ForceEngine::directPotential() {
    std::vector<double> potentials(particles.size(), 0);

    for (int i = 0; i < particles.size(); i++) {
        for (int j = 0; j < particles.size(); j++) {
            if (i == j) {
                continue;
            }
            double p = Particle::computePotential(particles.get(i), particles.get(j), softening);
            potentials[i] += p;
        }
    }
    return potentials;
}

std::vector<Eigen::Vector3d> ForceEngine::directForceOpt() {
    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));

    for (int i = 0; i < particles.size(); i++) {
        for (int j = i+1; j < particles.size(); j++) {
            Eigen::Vector3d f = Particle::computeForce(particles.get(i), particles.get(j), softening);
            forces[i] += f;
            forces[j] -= f;
        }
    }
    return forces;
}

std::vector<double> ForceEngine::directPotentialOpt() {
    std::vector<double> potentials(particles.size(), 0);

    for (int i = 0; i < particles.size(); i++) {
        for (int j = i+1; j < particles.size(); j++) {
            double p = Particle::computePotential(particles.get(i), particles.get(j), softening);
            potentials[i] += p;
            potentials[j] += p / particles.get(j).mass * particles.get(i).mass; // because computePotential returns the potential energy normalized by mass
        }
    }
    return potentials;
}


std::vector<Eigen::Vector3d> ForceEngine::treeForce(bool use_quad) {
    double radius = (fixed_radius == -1) ? particles.radius() : fixed_radius; // if particles get ejected it gets difficut for the tree // might want to fix a max radius

    std::vector<Eigen::Vector3d> forces(particles.size(), Eigen::Vector3d(0, 0, 0));
    Octree tree(radius); // doesn't care about outside of bunds particles
    tree.insert(particles);
    if (use_quad) {
        tree.computeQuadrupoles();
    }

    for (int i = 0; i < particles.size(); i++) {
        forces[i] = tree.getForce(particles.get(i), openingAngle);
    }
    
    return forces;
}

std::vector<double> ForceEngine::treePotential(bool use_quad) {
    double radius = (fixed_radius == -1) ? particles.radius() : fixed_radius; // if particles get ejected it gets difficut for the tree // might want to fix a max radius

    std::vector<double> potentials(particles.size(), 0);
    Octree tree(radius); // doesn't care about outside of bunds particles
    tree.insert(particles);
    if (use_quad) {
        tree.computeQuadrupoles();
    }

    for (int i = 0; i < particles.size(); i++) {
        potentials[i] = tree.getPotential(particles.get(i), openingAngle);
    }
    
    return potentials;
}


/**
 * -------------------------------
 * ! --- Integration Methods --- !
 * -------------------------------
 */

// forces get copied each time, this might not be optimal...

void ForceEngine::euler(double dt, Method method) {
    std::vector<Eigen::Vector3d> forces = computeForce(method);

    for (int i = 0; i < particles.size(); i++) {
        particles.get(i).position += dt * particles.get(i).velocity;
        particles.get(i).velocity += dt * forces[i] / particles.get(i).mass;
    }
}

void ForceEngine::symplectic(double dt,  Method method) {
    std::vector<Eigen::Vector3d> forces = computeForce(method);

    for (int i = 0; i < particles.size(); i++) {
        particles.get(i).velocity += dt * forces[i] / particles.get(i).mass;
        particles.get(i).position += dt * particles.get(i).velocity;
    }
}

void ForceEngine::leapfrog(double dt,  Method method) {
    std::vector<Eigen::Vector3d> forces = computeForce(method);

    for (int i = 0; i < particles.size(); i++) {
        particles.get(i).velocity += 0.5 * dt * forces[i] / particles.get(i).mass;
        particles.get(i).position += dt * particles.get(i).velocity;
        particles.get(i).velocity += 0.5 * dt * forces[i] / particles.get(i).mass;
    }
}



void ForceEngine::rk2(double dt, Method method) {
    std::vector<Eigen::Vector3d> forces = computeForce(method);
    std::vector<Eigen::Vector3d> k1_v, k1_r, k2_v, k2_r;

    // --- k1 ---
    for (int i=0; i<particles.size(); i++) {
        k1_v.push_back(forces[i] / particles.get(i).mass); // a
        k1_r.push_back(particles.get(i).velocity);         // v
    }

    // --- k2 ---
    // here we need to copy our particle set, since we want to move them in order to compute k2
    ParticleSet particles_copy = particles; // deep copy here

    // let's move our particles
    for (int i=0; i<particles.size(); i++) {
        particles_copy.get(i).velocity += dt * k1_v[i];
        particles_copy.get(i).position += dt * k1_r[i];
    }

    std::vector<Eigen::Vector3d> forces_k2 = ForceEngine(particles_copy).computeForce(method);

    for (int i=0; i<particles.size(); i++) {
        k2_v.push_back(forces_k2[i] / particles_copy.get(i).mass);
        k2_r.push_back(particles_copy.get(i).velocity);
    }

    // UPDATE
    for (int i=0; i<particles.size(); i++) {
        particles.get(i).position += dt * (k1_r[i] + k2_r[i]) / 2;
        particles.get(i).velocity += dt * (k1_v[i] + k2_v[i]) / 2;
    }
}

void ForceEngine::rk4(double dt, Method method) {
    std::vector<Eigen::Vector3d> forces = computeForce(method);
    std::vector<Eigen::Vector3d> k1_r, k1_v, k2_r, k2_v, k3_r, k3_v, k4_r, k4_v;

    // --- k1 ---
    for (int i = 0; i < particles.size(); i++) {
        k1_r.push_back(particles.get(i).velocity);                 
        k1_v.push_back(forces[i] / particles.get(i).mass);        
    }

    // --- k2 ---
    ParticleSet particles_copy = particles; // Deep copy for k2
    for (int i = 0; i < particles.size(); i++) {
        particles_copy.get(i).position += (dt / 2) * k1_r[i];
        particles_copy.get(i).velocity += (dt / 2) * k1_v[i];
    }

    std::vector<Eigen::Vector3d> forces_k2 = ForceEngine(particles_copy).computeForce(method);
    for (int i = 0; i < particles.size(); i++) {
        k2_r.push_back(particles_copy.get(i).velocity);           
        k2_v.push_back(forces_k2[i] / particles_copy.get(i).mass);
    }

    // --- k3 ---
    particles_copy = particles; // Reset particles_copy for k3
    for (int i = 0; i < particles.size(); i++) {
        particles_copy.get(i).position += (dt / 2) * k2_r[i];
        particles_copy.get(i).velocity += (dt / 2) * k2_v[i];
    }

    std::vector<Eigen::Vector3d> forces_k3 = ForceEngine(particles_copy).computeForce(method);
    for (int i = 0; i < particles.size(); i++) {
        k3_r.push_back(particles_copy.get(i).velocity);           
        k3_v.push_back(forces_k3[i] / particles_copy.get(i).mass); 
    }

    // --- k4 ---
    particles_copy = particles; // Reset particles_copy for k4
    for (int i = 0; i < particles.size(); i++) {
        particles_copy.get(i).position += dt * k3_r[i];
        particles_copy.get(i).velocity += dt * k3_v[i];
    }

    std::vector<Eigen::Vector3d> forces_k4 = ForceEngine(particles_copy).computeForce(method);
    for (int i = 0; i < particles.size(); i++) {
        k4_r.push_back(particles_copy.get(i).velocity);           // r' = v_end
        k4_v.push_back(forces_k4[i] / particles_copy.get(i).mass); // v' = a_end
    }

    // --- Update ---
    for (int i = 0; i < particles.size(); i++) {
        particles.get(i).position += (dt / 6.0) * (k1_r[i] + 2 * k2_r[i] + 2 * k3_r[i] + k4_r[i]);
        particles.get(i).velocity += (dt / 6.0) * (k1_v[i] + 2 * k2_v[i] + 2 * k3_v[i] + k4_v[i]);
    }
}



/**
 * -----------------
 * ! --- Utils --- !
 * -----------------
 */

std::string methodToString(Method method) {
    switch (method) {
        case Method::direct:
            return "direct";
        
        case Method::tree_mono:
            return "tree-monopole";

        case Method::tree_quad:
            return "tree-quadrupole";
        
        case Method::direct_opt:
            return "direct-optimized";

        default:
            throw std::runtime_error("Method not implemented.");
    }
}

std::string methodToString(IntegrationMethod method) {
    switch (method) {
        case IntegrationMethod::euler:
            return "euler";
        
        case IntegrationMethod::leapfrog:
            return "leapfrog";
        
        case IntegrationMethod::rk2:
            return "rk2";

        case IntegrationMethod::rk4:
            return "rk4";

        case IntegrationMethod::symplectic:
            return "symplectic";

        default:
            throw std::runtime_error("Method not implemented.");
    }
}