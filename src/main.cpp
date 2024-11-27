
#include "message.hpp"
#include <Eigen/Dense>

#include <iostream>
#include "particle.hpp"
#include "force_engine.hpp"

int tutoTimer() {
    Timer timer("Summing integers up to 1e6");
    timer.start();
    long long s=0;
    for (long long i = 0; i < 1e5; i++) { // 1e10 is too large for long long
        s+=i; 
    }
    timer.stop();
    timer.display();
    Message::print("Sum is: " + std::to_string(s));
    Message("Hello, World!", "#");
    return 0;
}


int tutoEigen() {
    /**
     * ################
     * ###  VECTORS ###
     * ################
     */

    Eigen::Vector3d v1(1, 2, 3); // equivalent to Eigen::VectorXd v1(3); v1 << 1, 2, 3; // allows you to have a dinamic dimension
    std::cout << "3D Vector:\n" << v1 << std::endl;

    /**
     * ################
     * ### MATRICES ###
     * ################
     */

    Eigen::Matrix<double, 2, 3> mat;
    mat << 1, 2, 3,
           4, 5, 6;
    std::cout << "2x3 Matrix:\n" << mat << std::endl;

    /**
     * ################
     * ### PRODUCTS ###
     * ################
     */
    Eigen::Vector2d result = mat * v1; // this works with vectors, matrices, reals
    std::cout << "Matrix-Vector Product:\n" << result << std::endl;

    /**
     * ###################
     * ### OTHER STUFF ###
     * ###################
     */
    Eigen::Matrix<double, 3, 3> mat2 = mat.transpose() * mat;
    Eigen::Matrix<double, 3, 3> mat3 = mat2.inverse();

    /**
     * #####################
     * ### INSTANTIATION ###
     * #####################
     */
    Eigen::Matrix3d zero_matrix = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d random_matrix = Eigen::Matrix3d::Random();
    Eigen::Matrix3d identity_matrix = Eigen::Matrix3d::Identity();
    Eigen::Vector3d ones_vector = Eigen::Vector3d::Ones();
    Eigen::Vector3d zero_vector = Eigen::Vector3d::Zero();
    zero_vector(0) = 1;
    Message("Done!", "#");

    //avoid some warnings:
    std::cout << mat3(0,0) << std::endl;
    std::cout << zero_matrix(0,0) << std::endl;
    std::cout << random_matrix(0,0) << std::endl;
    std::cout << identity_matrix(0,0) << std::endl;
    std::cout << ones_vector(0) << std::endl;


    return 0;

}





int testDirectSolver() {
    double m = 1;
    Particle p1(Eigen::Vector3d(0.97, -0.243, 0), Eigen::Vector3d(0.4662, 0.4324, 0), m); // these are values for stable 3 body problem
    Particle p2(Eigen::Vector3d(-0.97, 0.243, 0), Eigen::Vector3d(0.4662, 0.4324, 0), m);
    Particle p3(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(-0.9324, -0.8647, 0), m);

    ParticleSet particles;
    ParticleSet particles_over_time;

    particles.add({p1, p2, p3});
    particles_over_time.add(particles); // a copy is made when particle is added! this is good since we want to store them at each time step!
   

    ForceEngine engine(particles);
    double dt = 1e-2;

    for (int i = 0; i < 5000; i++) {
        std::vector<Eigen::Vector3d> forces = engine.computeForce(Method::direct);
        particles.applyForces(forces);
        particles.update(dt, IntegrationMethod::euler);
        particles_over_time.add(particles);
    }
    ParticleSet::save("files/direct_N_part_trajectory.txt", particles_over_time);
    
    return 0;
}

int main() {
    //tutoParticles();
    testDirectSolver();
    return 0;
}
