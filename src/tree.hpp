

#pragma once 

#include <Eigen/Dense>
#include <vector>
#include "particle.hpp"
#include <string>
#include <fstream>


/**
 * ############
 * ### NODE ###
 * ############ 
 */

class Node {
    public:
    Eigen::Vector3d position; // position of the center of the cube
    Eigen::Vector3d centerOfMass; // center of mass of the particles in the node
    double totalMass;
    double halfWidth;
    Node* children[8]; // 8 children nodes
    Particle* particle; // pointer to the particle in the node (must be allowed to be null, thus we need a pointer here!)


    Node(Eigen::Vector3d position, double halfWidth);

    /**
     * Destructor: Recursively delete all children and set all pointers to null.
     */
    ~Node();


    /**
     * Here is the procedure:
     * Update total mass and center of mass
     * 
     * 1. If the node is a leaf (no children):
     *    a. Empty Leaf: Put particle there
     *    b. Occupied Leaf: Create 8 children / do 2. (call insert again)
     * 2. If the node is not a leaf (or not anymore):
     *    a. Find the correct child / take this child / insert on it
     */
    void insert(Particle& particle);

    /**
     * Here is the procedure:
     * 1. If the node is a leaf:
     *   a. Return the force of the particle in the node
     * 2. If the node is not a leaf:
     *   a. If the opening angle condition is satisfied, return the force of the center of mass
     *   b. If the opening angle condition is not satisfied, return the sum of the forces of the children
     */
    Eigen::Vector3d getForce(Particle& particle, double theta);

    static void save(std::ofstream& file, Node* node);

    /**
     * get recusrsively the particles in the tree
     */
    ParticleSet getParticles();


private:
    /**
     * int childIndex = (particle.x > node.center.x) |
                 ((particle.y > node.center.y) << 1) |
                 ((particle.z > node.center.z) << 2); // magic <3
     */
    int childIndex(Particle& particle);

    /**
     * Checks whether the first child is a null pointer
     */
    bool isLeaf();

    /**
     * Creates the 8 children of the node
     */
    void createChildren();

    /**
     * Updates the center of mass and total mass of the node
     */
    void updateMass(Particle& particle);
};


// for the data, the depth of the tree is 14 <3
class Octree {
    public:
    Node* root;
    static double openingAngle;

    Octree(Eigen::Vector3d position, double halfWidth);
    Octree(double halfWidth);
    ~Octree();

    /**
     * Calls *clear()*
     * Inserts a set of particles into the tree.
     * All particles we get copied, I don't understand why but if i don't do this my pointers get deleted.
     */
    void insert(ParticleSet particles);

    /**
     * get Force on a particle by descending the tree until opening angle condition is satisfied.
     */
    Eigen::Vector3d getForce(Particle& particle);

    /**
     * get all particles in the tree (for debugging purposes)
     */
    ParticleSet getParticles();

public:
    /**
     * Saves recursively the nodes into a csv file
     */
    static void save(std::string path, Octree& tree);

private:
    /**
     * Resets the Tree to its initial state
     */
    void clear();

};

