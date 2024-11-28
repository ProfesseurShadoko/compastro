


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
    // update the center of mass and total mass
    updateMass(particle);

    // if the node is a leaf
    if (isLeaf()) {
        // if the leaf is empty
        if (this->particle == nullptr) {
            this->particle = &particle;
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



/**
 * ####################
 * ### OCTREE CLASS ###
 * ####################
 */

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