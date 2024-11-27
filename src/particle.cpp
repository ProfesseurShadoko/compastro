
#include "particle.hpp"
#include <iostream>
#include <fstream>


/**
 * ######################
 * ### PARTICLE CLASS ###
 * ######################
 */

Particle::Particle() {
    this->position = Eigen::Vector3d(0, 0, 0);
    this->velocity = Eigen::Vector3d(0, 0, 0);
    this->acceleration = Eigen::Vector3d(0, 0, 0);
    this->mass = 1;
}

Particle::Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass) {
    this->position = position;
    this->velocity = velocity;
    this->acceleration = Eigen::Vector3d(0, 0, 0);
    this->mass = mass;
}

Particle::Particle(const Particle& p) {
    this->position = p.position; // this gets effectively copied
    this->velocity = p.velocity;
    this->acceleration = p.acceleration;
    this->mass = p.mass;
}

void Particle::applyForce(Eigen::Vector3d force) {
    Eigen::Vector3d f = force / mass;
    acceleration += f;
}

void Particle::resetForce() {
    acceleration = Eigen::Vector3d(0, 0, 0);
}

void Particle::update(double dt, IntegrationMethod method) {
    switch (method)
    {
    case IntegrationMethod::euler:
        velocity += acceleration * dt;
        position += velocity * dt;/* code */
        break;
    case IntegrationMethod::leapfrog:
        velocity += acceleration * dt / 2;
        position += velocity * dt;
        velocity += acceleration * dt / 2;
        break;
    
    default:
        throw std::runtime_error("Integration method not implemented.");
        break;
    }
    
    resetForce();
}

void Particle::display() {
    std::cout << "Particle:" << std::endl;
    std::cout << "\t- Position: " << position.transpose() << std::endl;
    std::cout << "\t- Velocity: " << velocity.transpose() << std::endl;
    std::cout << "\t- Acceleration: " << acceleration.transpose() << std::endl;
    std::cout << "\t- Mass: " << mass << std::endl;
    std::cout << std::endl;
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

Particle& ParticleSet::get(int i) {
    return particles[i];
}

void ParticleSet::add(Particle p) {
    particles.push_back(p);
}

void ParticleSet::add(ParticleSet ps) {
    for (int i = 0; i < ps.size(); i++) {
        add(ps.get(i));
    }
}

void ParticleSet::add(std::vector<Particle> ps) {
    for (size_t i = 0; i < ps.size(); i++) {
        add(ps[i]);
    }
}


int ParticleSet::size() {
    return particles.size();
}

// particle functions
void ParticleSet::applyForces(std::vector<Eigen::Vector3d> forces) {
    // check the forces has the right size
    if (forces.size() != particles.size()) {
        throw std::runtime_error("Forces vector has wrong size.");
    }
    for (size_t i = 0; i < particles.size(); i++) {
        get(i).applyForce(forces[i]);
    }
}

void ParticleSet::resetForces() {
    for (size_t i = 0; i < particles.size(); i++) {
        get(i).resetForce();
    }
}

void ParticleSet::update(double dt, IntegrationMethod method) {
    for (size_t i = 0; i < particles.size(); i++) {
        get(i).update(dt, method);
    }
}


ParticleSet ParticleSet::load(std::string filename) {
    ParticleSet particles = ParticleSet();
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }

    // each line contains index, mass, x, y, z, vx, vy, vz, eps, phi // but we do not care about eps and phi
    int index;
    double mass, x, y, z, vx, vy, vz;
    while (file >> index >> mass >> x >> y >> z >> vx >> vy >> vz) {
        Particle p(Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz), mass);
        particles.add(p);
    }

    return particles;
}

void ParticleSet::save(std::string filename, ParticleSet particles) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }

    for (int i = 0; i < particles.size(); i++) {
        Particle p = particles.get(i);
        file << i << " " << p.mass << " " << p.position(0) << " " << p.position(1) << " " << p.position(2) << " " << p.velocity(0) << " " << p.velocity(1) << " " << p.velocity(2) << " 0 0" << std::endl;
    }
}
