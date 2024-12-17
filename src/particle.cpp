
#include "particle.hpp"
#include <iostream>
#include <fstream>
#include <iomanip> // precision of csv files
#include "force_engine.hpp" // for softening


/**
 * ######################
 * ### PARTICLE CLASS ###
 * ######################
 */

long long Particle::forceCallCounter = 0; 
long long Particle::id_counter = 0;

Particle::Particle() {
    this->position = Eigen::Vector3d(0, 0, 0);
    this->velocity = Eigen::Vector3d(0, 0, 0);
    this->mass = 1;
    this->id = id_counter++;
}

Particle::Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass) {
    this->position = position;
    this->velocity = velocity;
    this->mass = mass;
    this->id = id_counter++;
}

Particle::Particle(const Particle& p) {
    this->position = p.position; // this gets effectively copied
    this->velocity = p.velocity;
    this->mass = p.mass;
    this->current_time = p.current_time;
    this->id = p.id;
}

Particle::Particle(std::initializer_list<double> init) {
    if (init.size() != 3) {
        throw std::runtime_error("Initializer list must have 3 elements.");
    }
    this->position = Eigen::Vector3d(*(init.begin()), *(init.begin() + 1), *(init.begin() + 2));
    this->velocity = Eigen::Vector3d(0, 0, 0);
    this->mass = 1;
    this->id = id_counter++;
}

void Particle::updateCurrentTime(double dt) {
    current_time += dt;
}

void Particle::display() {
    std::cout << "Particle <" + std::to_string(id) + ">:" << std::endl;
    std::cout << "\t- Position: " << position.transpose() << std::endl;
    std::cout << "\t- Velocity: " << velocity.transpose() << std::endl;
    std::cout << "\t- Mass: " << mass << std::endl;
    std::cout << std::endl;
}

int Particle::getId() const {
    return id;
}



Eigen::Vector3d Particle::computeForce(Particle& particle, Particle& p_attractor, double eps) {
    /*
    
    Eigen::Vector3d r = particle.position - p_attractor.position; // u_r goes from p_attractor to our particle
    return - (particle.mass * p_attractor.mass) / (pow(r.norm(), 2) + eps * eps) * r.normalized(); // F = -GM/(r+eps)^2
    // -u_r goes from particle to p_attractor => particle is moved towards p_attractor!*/ // actually this is incorrect
    forceCallCounter++;
    Eigen::Vector3d r = particle.position - p_attractor.position; // u_r goes from p_attractor to our particle
    double r_squared = r.squaredNorm() + eps * eps;              // r^2 + eps^2
    double r_cubed = sqrt(r_squared) * r_squared;                // (r^2 + eps^2)^(3/2)
    return - (particle.mass * p_attractor.mass) / r_cubed * r;   // F = -GM/(r^2 + eps^2)^(3/2) * r
}

double Particle::computePotential(Particle& particle, Particle& p_attractor, double eps) {
    double r = (particle.position - p_attractor.position).norm();
    return - p_attractor.mass / pow(pow(r, 2) + pow(eps, 2), 0.5);
}


long long Particle::getForceCallCounter() {
    return forceCallCounter;
}

long long Particle::popForceCallCounter() {
    long long out = getForceCallCounter();
    forceCallCounter = 0;
    return out;
}

/**
 * #########################
 * ### PARTICLESET CLASS ###
 * #########################
 */


// set functions
ParticleSet::ParticleSet() {
    this->particles = std::vector<Particle>();
}

ParticleSet::ParticleSet(std::vector<Particle> particles) {
    this->particles = particles; // here a copy is made
}

ParticleSet::ParticleSet(std::initializer_list<Particle> init) {
    this->particles = std::vector<Particle>(init);
}

ParticleSet::ParticleSet(const ParticleSet& ps) {
    this->particles = ps.particles; // this is a deep copy, because std::vector does deep copy on its own
}

Particle& ParticleSet::get(int i) {
    return particles[i];
}

ParticleSet ParticleSet::slice(int start, int end) {
    ParticleSet ps;

    if (end < start) {
        // this means add all of them
        end = particles.size();
    }

    for (int i = start; i < end; i++) {
        ps.add(get(i));
    }
    return ps;
}



void ParticleSet::add(Particle p) {
    particles.push_back(p);
}

void ParticleSet::add(ParticleSet ps) {
    for (int i = 0; i < ps.size(); i++) {
        add(ps.get(i)); // this makes a copy because ParticleSet::add makes a copy!
    }
}

void ParticleSet::add(std::vector<Particle> ps) {
    for (size_t i = 0; i < ps.size(); i++) {
        add(ps[i]);
    }
}

void ParticleSet::add(std::initializer_list<Particle> init) {
    for (auto p : init) {
        add(p);
    }
}


int ParticleSet::size() {
    return particles.size();
}

double ParticleSet::radius() {
    double radius = 0;
    for (size_t i = 0; i < particles.size(); i++) {
        double r = particles[i].position.cwiseAbs().maxCoeff();
        if (r > radius) {
            radius = r;
        }
    }
    return radius;
}

void ParticleSet::com() {
    double m_tot = 0;
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    Eigen::Vector3d system_vel = Eigen::Vector3d::Zero();

    for (int i=0; i<size(); i++) {
        m_tot += get(i).mass;
        com += get(i).mass * get(i).position;
        system_vel += get(i).mass * get(i).velocity;
    }

    com /= m_tot;
    system_vel /= m_tot;

    // let's move everyone
    for (int i=0; i<size(); i++) {
        get(i).position -= com;
        get(i).velocity -= system_vel;
    }

}


void ParticleSet::updateCurrentTime(double dt) {
    for (size_t i = 0; i < particles.size(); i++) {
        get(i).updateCurrentTime(dt);
    }
}


ParticleSet ParticleSet::load(std::string path) {
    ParticleSet particles = ParticleSet();
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }

    // each line contains index, mass, x, y, z, vx, vy, vz, eps, phi // but we do not care about eps and phi
    int index;
    double mass, x, y, z, vx, vy, vz, eps, phi;
    while (file >> index >> mass >> x >> y >> z >> vx >> vy >> vz >> eps >> phi) { // we actually don't care about espilon and phi
        Particle p(Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz), mass);
        particles.add(p);
    }

    return particles;
}

void ParticleSet::save(std::string path, ParticleSet particles) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }
    file << "index,mass,x,y,z,vx,vy,vz,eps,potential,t" << std::endl; // phi becomes now just a stupid name for time

    file << std::fixed << std::setprecision(30); // need high precision to evaluate rk2, who has precision up to 1e-15
    for (int i = 0; i < particles.size(); i++) {
        Particle p = particles.get(i);
        file << p.getId() << "," << p.mass << "," << p.position(0) << "," << p.position(1) << "," << p.position(2) << "," << p.velocity(0) << "," << p.velocity(1) << "," << p.velocity(2) << "," << ForceEngine::softening << "," << p.potentialEnergy << "," << p.current_time << std::endl;
    }
}


void ParticleSet::saveForces(std::string path, std::vector<Eigen::Vector3d> forces) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }
    file << "index,fx,fy,fz" << std::endl;

    file << std::fixed << std::setprecision(30);
    for (size_t i = 0; i < forces.size(); i++) {
        file << i << "," << forces[i](0) << "," << forces[i](1) << "," << forces[i](2) << std::endl;
    }
}


std::vector<ParticleSet> ParticleSet::split() {
    double last_time_added = -1;
    std::vector<ParticleSet> sets;

    for (int i = 0; i < size(); i++) {
        if (particles[i].current_time == last_time_added) {
            sets.back().add(get(i));
            continue;
        }
        if (particles[i].current_time < last_time_added) {
            throw std::runtime_error("Particles are not sorted in time.");
        }
        last_time_added = get(i).current_time;
        sets.push_back(ParticleSet());
        sets.back().add(get(i));
    }
    
    return sets;
}
