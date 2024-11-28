

#include <Eigen/Dense>
#include <vector>
#include "particle.hpp"


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



class Octree {
    public:
    Node* root;

    Octree(Eigen::Vector3d position, double halfWidth);

    /**
     * Calls *clear()*
     * Inserts a set of particles into the tree.
     */
    void insert(ParticleSet particles);

private:
    /**
     * Resets the Tree to its initial state
     */
    void clear();

};