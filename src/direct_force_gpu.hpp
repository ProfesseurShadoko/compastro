#include "message.hpp"
#include "particle.hpp"
#include <CL/cl.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>


class DirectForceGPU {
    private:
        size_t dataSize;
        cl_uint numPlatforms;
        cl_platform_id platform;
        cl_device_id device;
        cl_int err;
        cl_context context;
        cl_command_queue queue;
        cl_mem bufferX;
        cl_mem bufferY;
        cl_mem bufferZ;
        cl_mem bufferMass;
        cl_mem bufferForceX;
        cl_mem bufferForceY;
        cl_mem bufferForceZ;
        std::string kernelCode;
        cl_program program;
        cl_kernel kernel;
        size_t globalSize;
        
        ParticleSet& particles;
    
    public:
        DirectForceGPU(ParticleSet& particles);
        ~DirectForceGPU();	

        void setBuffers();
        void release();
        std::vector<Eigen::Vector3d> computeForces();

        static void checkDevices();
        static void checkError(cl_int err, std::string operation);
        static std::string loadKernel(std::string filename);
};