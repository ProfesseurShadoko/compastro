

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

    Eigen::Matrix3d quadrupole; // quadrupole moment of the particles in the node
    /**
     * reminder on the quadrupol expression:
     * Q_ij = sum(m_k * (3 * r_k,i * r_k,j - |r_k|^2 * delta_ij)) # r = distance of particle to center of mass
     */


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
     * 
     * Returns the force acting on the input particle.
     */
    Eigen::Vector3d getForce(Particle& particle, double theta, bool useQuadrupoles);


    double getPotential(Particle& particle, double theta, bool useQuadrupoles);

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
     * # first pass!
     * Updates the center of mass and total mass of the node.
     */
    void updateMass(Particle& particle);
    
public:
     /**
      * # second pass!
     * Qij = sum(mk[3(s-xk)_i * (s-xk)_j - |s-xk|^2 * delta_ij]) with s the position of the center of mass => you need to have it already computed, thus the necessity of a second pass
     * 
     * The second pass will compute all quadrupoles of the nodes in the tree. Let's go step by step in depth:
     * 1. First quadrupole is computed with all N particles => O(N)
     * 2. Next 8 quadrupoles are computed with N/8 particles each approx => 8*N/8 = N => O(N)
     * 3. Next 64 quadrupoles are computed with N/64 particles each approx => 64*N/64 = N => O(N)
     * 4. ...
     * log_8(N). Here we stop => in total O(N log N) <3
     * 
     * Note: for each specific particle, no issue that the resulting force gets computed with a quadrupole involving that exact particle. Indeed, if we go thrgouh the nodes that contain the specific particle, we will go all the way done to the leaf. Thus, no issue here.
     * 
     * The second pass is called recursively by the previous node. However, the previous node needs all particles below him. Thus, we need to return a vector of particles.
     */
    std::vector<Particle*> secondPass();
};


// for the data, the depth of the tree is 14 <3
class Octree {
    public:
    Node* root;

private:
    bool useQuadrupoles;
public:


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
    Eigen::Vector3d getForce(Particle& particle, double openingAngle);

    /**
     * get Potential on a particle by descending the tree until opening angle condition is satisfied.
     */
    double getPotential(Particle& particle, double openingAngle);

    /**
     * get all particles in the tree (for debugging purposes)
     */
    ParticleSet getParticles();

    /**
     * Precompute all quadrupoles in the tree # second pass. If this function is called, the quadrupoles will be used in the force computation.
     */
    void computeQuadrupoles();

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

