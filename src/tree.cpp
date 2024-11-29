


#include "tree.hpp"



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
Eigen::Vector3d Node::getForce(Particle& particle, double theta) {
    
    // if the node is a leaf
    if (isLeaf()) {
        
        if (this->particle == nullptr || this->particle->position == particle.position) {
            return Eigen::Vector3d(0, 0, 0); // no force if no particle / a particle doesn't apply a force on itself
        }
        return Particle::computeForce(particle, *this->particle); // this is for now first order, we will add quadrupole later
    }

    // if the node is not a leaf
    double distance = (centerOfMass - particle.position).norm();
    if (halfWidth / distance < theta) {
        Particle pseudoParticle = Particle(centerOfMass, Eigen::Vector3d(0, 0, 0), totalMass);
        return Particle::computeForce(particle, pseudoParticle);
    }

    Eigen::Vector3d force(0, 0, 0);
    for (int i = 0; i < 8; i++) {
        force += children[i]->getForce(particle, theta);
    }
    return force;
}



/**
 * ####################
 * ### OCTREE CLASS ###
 * ####################
 */

double Octree::openingAngle = 0.5; // default value for the opening angle ~ good accuracy here, for faster go to 1.0


Octree::Octree(Eigen::Vector3d position, double halfWidth) {
    this->root = new Node(position, halfWidth);
}

Octree::Octree(double halfWidth) {
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

Eigen::Vector3d Octree::getForce(Particle& particle) {
    return root->getForce(particle, Octree::openingAngle);
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
            particles.push_back(particle);
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
        Eigen::Vector3d r = particles[k]->position - centerOfMass; // this is s-xk
        quadrupole += particles[k]->mass * (3 * r * r.transpose() - r.squaredNorm() * Eigen::Matrix3d::Identity());
    }
    
    return particles;
}

void Octree::computeQuadrupoles() {
    root->secondPass();
}

