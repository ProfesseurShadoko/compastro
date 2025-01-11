


#include "tree.hpp"
#include "force_engine.hpp" // for softening



/**
 * ##################
 * ### NODE CLASS ###
 * ##################
 */

Node::Node(Eigen::Vector3d position, double halfWidth) {
    this->position = position;
    this->halfWidth = halfWidth;
    this->centerOfMass = Eigen::Vector3d(0, 0, 0);
    this->totalMass = 0;
    this->particle = nullptr;
    this->quadrupole = Eigen::Matrix3d::Zero();

    for (int i = 0; i < 8; i++) {
        this->children[i] = nullptr;
    }
}


Node::~Node() { // the only important thing is to call delete on everyone => recursively destory children
    for (int i = 0; i < 8; i++) {
        if (children[i] != nullptr) {
            delete children[i];
            children[i] = nullptr;
        }
    }

    if (particle != nullptr) {
        delete particle;  // Assuming particle is allocated dynamically with `new`
        particle = nullptr;
    }
}

void Node::insert(Particle& particle) {
    // check whether the particle is in the node
    if (particle.position.x() <= position.x() - halfWidth || particle.position.x() > position.x() + halfWidth ||
        particle.position.y() <= position.y() - halfWidth || particle.position.y() > position.y() + halfWidth ||
        particle.position.z() <= position.z() - halfWidth || particle.position.z() > position.z() + halfWidth) {
        return;
    }

    // update the center of mass and total mass
    updateMass(particle);

    // if the node is a leaf
    if (isLeaf()) {
        // if the leaf is empty
        if (this->particle == nullptr) {
            this->particle = new Particle(particle); // make a copy an allocate memory, we don't want the particle to get erased
        } else {
            // if the leaf is occupied
            Particle* oldParticle = this->particle;
            this->particle = nullptr;
            createChildren();
            // insert the particle again
            children[childIndex(*oldParticle)]->insert(*oldParticle);
            children[childIndex(particle)]->insert(particle);
        }
    } else {
        // if the node is not a leaf
        children[childIndex(particle)]->insert(particle);
    }
}


void Node::updateMass(Particle& particle) {
    centerOfMass = (centerOfMass * totalMass + particle.position * particle.mass) / (totalMass + particle.mass);
    totalMass += particle.mass;
}

bool Node::isLeaf() {
    return children[0] == nullptr;
}

void Node::createChildren() {
    for (int i = 0; i < 8; i++) {
        Eigen::Vector3d childPosition = position;
        childPosition.x() += halfWidth * (i & 1 ? 0.5 : -0.5);
        childPosition.y() += halfWidth * (i & 2 ? 0.5 : -0.5);
        childPosition.z() += halfWidth * (i & 4 ? 0.5 : -0.5);
        children[i] = new Node(childPosition, halfWidth / 2);
    }
}

int Node::childIndex(Particle& particle) {
    return (particle.position.x() > position.x()) |
           ((particle.position.y() > position.y()) << 1) |
           ((particle.position.z() > position.z()) << 2);
}

//                    __
// force computation (:D)
//                    \/
#include <iostream>

Eigen::Vector3d Node::getForce(Particle& particle, double theta, bool useQuadrupoles) {
    
    /**
     * ############
     * ### LEAF ###
     * ############ => simply compute the usual force
     */
    if (isLeaf()) {
        
        if (this->particle == nullptr || *(this->particle) == particle) { // id comparison
            return Eigen::Vector3d(0, 0, 0); // no force if no particle / a particle doesn't apply a force on itself
        }
        return Particle::computeForce(particle, *this->particle, ForceEngine::softening); // this is for now first order, we will add quadrupole later
    }

    /**
     * ##################
     * ### NOT A LEAF ###
     * ################## => check opening angle condition
     */
    double distance = (centerOfMass - particle.position).norm();
    if (halfWidth / distance < theta) { 
        /**
         * small parenthesis on theta: what happens when we compare a particle to itself?
         * well centerOfMass is computed with the particle => particle is in the cube => distance to center should be less than halfwidth (except if center of mass is unluckily far away from cube center)
         * halfWidth / distance > 1 => setting theta to more than 1 will bring isues here! theta should be at most 1 for debug
         */

        // ### OPENING ANGLE CONDITION SATISFIED ### // => return approximationwith center of mass (and potentially quadrupole)
        Eigen::Vector3d r_tilde = particle.position - centerOfMass; // r_tilde = r-com
        Particle pseudoParticle(centerOfMass, Eigen::Vector3d::Zero(), totalMass); // pseudo particle at center of mass
        Eigen::Vector3d monopole_force = Particle::computeForce(particle, pseudoParticle, ForceEngine::softening); // F = -GM/r^2
        
        if (!useQuadrupoles) { // maybe we don't want to use quadrupoles yet
            return monopole_force;
        }   

        // compute the quadrupole force
        double r_tilde_norm = r_tilde.norm();
        double r_tilde5 = r_tilde_norm * r_tilde_norm * r_tilde_norm * r_tilde_norm * r_tilde_norm;
        double r_tilde7 = r_tilde5 * r_tilde_norm * r_tilde_norm;

        Eigen::Vector3d q_term = 2 * quadrupole * r_tilde / r_tilde5; // (2Qr / r^5)
        double rQr = r_tilde.transpose() * quadrupole * r_tilde;
        Eigen::Vector3d r_term = 5 * rQr * r_tilde / r_tilde7; // (5rQr * r / r^7)
        Eigen::Vector3d quadrupole_force = 0.5 * particle.mass * (q_term - r_term); // -GM(2Qr / r^5 - 5rQr * r / r^7) // - !!! <3 :( => (u/v)' = u'v - uv' / v^2 !!!!! et pas +
        return monopole_force + quadrupole_force;
    }

    // ### OPENING ANGLE CONDITION NOT SATISFIED ### // => return sum of forces of children
    Eigen::Vector3d force(0, 0, 0);
    for (int i = 0; i < 8; i++) {
        force += children[i]->getForce(particle, theta, useQuadrupoles);
    }
    return force;
}


double Node::getPotential(Particle& particle, double theta, bool useQuadrupoles) {
    
    /**
     * ############
     * ### LEAF ###
     * ############ => simply compute the usual potential
     */
    if (isLeaf()) {
        
        if (this->particle == nullptr || *(this->particle) == particle) { // id comparison
            return 0; // no potential if no particle / a particle doesn't apply a potential on itself
        }
        return Particle::computePotential(particle, *this->particle, ForceEngine::softening); // this is for now first order, we will add quadrupole later
    }

    /**
     * ##################
     * ### NOT A LEAF ###
     * ################## => check opening angle condition
     */
    double distance = (centerOfMass - particle.position).norm();
    if (halfWidth / distance < theta) { 
        
        /**
         * small parenthesis on theta: what happens when we compare a particle to itself?
         * well centerOfMass is computed with the particle => particle is in the cube => distance to center should be less than halfwidth (except if center of mass is unluckily far away from cube center)
         * halfWidth / distance > 1 => setting theta to more than 1 will bring isues here! theta should be at most 1 for debug
         */

        // ### OPENING ANGLE CONDITION SATISFIED ### // => return approximationwith center of mass (and potentially quadrupole)
        Eigen::Vector3d r_tilde = particle.position - centerOfMass; // r_tilde = r-com
        Particle pseudoParticle(centerOfMass, Eigen::Vector3d::Zero(), totalMass); // pseudo particle at center of mass
        double monopole_force = Particle::computePotential(particle, pseudoParticle, ForceEngine::softening); // F = -GM/r^2
        if (!useQuadrupoles) { // maybe we don't want to use quadrupoles yet
            return monopole_force;
        }   
        
        double rQr = r_tilde.transpose() * quadrupole * r_tilde;
        double r_tilde_norm = r_tilde.norm();
        double r_tilde5 = r_tilde_norm * r_tilde_norm * r_tilde_norm * r_tilde_norm * r_tilde_norm;
        double quadrupole_force = - 0.5 * rQr / r_tilde5; // - 1/2 (rQr / r^5)
        return monopole_force + quadrupole_force; // we choose to always return the normalized potential (like electrmagntism with charge). for actual potential, multiply by particle.mass
    }

    // ### OPENING ANGLE CONDITION NOT SATISFIED ### // => return sum of forces of children
    double potential = 0;
    for (int i = 0; i < 8; i++) {
        potential += children[i]->getPotential(particle, theta, useQuadrupoles);
    }
    return potential;
}



/**
 * ####################
 * ### OCTREE CLASS ###
 * ####################
 */

Octree::Octree(Eigen::Vector3d position, double halfWidth) {
    this->useQuadrupoles = false;
    this->root = new Node(position, halfWidth);
}

Octree::Octree(double halfWidth) {
    this->useQuadrupoles = false;
    this->root = new Node(Eigen::Vector3d(0, 0, 0), halfWidth);
}

Octree::~Octree() {
    delete root;
}   

void Octree::insert(ParticleSet particles) {
    this->clear(); // useless for the first iteration, necessary after though
    for (int i = 0; i < particles.size(); i++) {
        root->insert(particles.get(i));
    }
}

void Octree::clear() {
    Eigen::Vector3d position = root->position;
    double halfWidth = root->halfWidth;
    delete root;
    root = new Node(position, halfWidth);
}

Eigen::Vector3d Octree::getForce(Particle& particle, double openingAngle) {
    return root->getForce(particle, openingAngle, useQuadrupoles);
}

double Octree::getPotential(Particle& particle, double openingAngle) {
    return root->getPotential(particle, openingAngle, useQuadrupoles);
}






/**
 * ##############
 * ### TO CSV ###
 * ##############
 */

void Octree::save(std::string path, Octree& tree) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file.");
    }
    // set column names separated by commas (csv format)
    file << "x,y,z,halfWidth,totalMass,centerOfMass_x,centerOfMass_y,centerOfMass_z" << std::endl;
    Node::save(file, tree.root);
}

void Node::save(std::ofstream& file, Node* node) {
    file << node->position.x() << "," << node->position.y() << "," << node->position.z() << "," << node->halfWidth << "," << node->totalMass << "," << node->centerOfMass.x() << "," << node->centerOfMass.y() << "," << node->centerOfMass.z() << std::endl;
    if (!node->isLeaf()) {
        for (int i = 0; i < 8; i++) {
            Node::save(file, node->children[i]);
        }
    }
}


/**
 * ###################
 * ### QUADRUPOLES ###
 * ###################
 */

ParticleSet Octree::getParticles() {
    return root->getParticles();
}

ParticleSet Node::getParticles() {
    ParticleSet particles;
    if (isLeaf()) {
        if (particle != nullptr) {
            particles.add(*particle);
        }
    } else {
        for (int i = 0; i < 8; i++) {
            ParticleSet childParticles = children[i]->getParticles();
            particles.add(childParticles);
        }
    }
    return particles;
}



std::vector<Particle*> Node::secondPass() {
    std::vector<Particle*> particles;
    if (isLeaf()) { // if empty, quadrupole is zero. if not empty, center of mass is particle, and thus quadurpole is also zero!
        if (particle != nullptr) {
            particles.push_back(particle); // we need to return all sub particles to the parent node for recursion
        }
        // quadrupole = Eigen::Matrix3d::Zero(); // is already initialized to zero for everyone
        return particles;
    }

    for (int i = 0; i < 8; i++) {
        std::vector<Particle*> childParticles = children[i]->secondPass(); // compute children quadrupoles and collects all particles below the child.
        particles.insert(particles.end(), childParticles.begin(), childParticles.end()); // concateante!
    }
    // all quadrupoles are computed below / we have collected all particles below => we can compute the quadrupole
    // quadrupole = Eigen::Matrix3d::Zero();
    for (size_t k = 0; k < particles.size(); k++) {
        Eigen::Vector3d r_tilde = particles[k]->position - centerOfMass; // this is s-xk
        quadrupole += particles[k]->mass * (3 * r_tilde * r_tilde.transpose() - r_tilde.squaredNorm() * Eigen::Matrix3d::Identity());
    }
    
    return particles;
}

void Octree::computeQuadrupoles() {
    if (useQuadrupoles) {
        throw std::runtime_error("Quadrupoles have already been computed.");
    }
    root->secondPass();
    this->useQuadrupoles = true;
}

