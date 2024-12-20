#include "direct_force_gpu.hpp"
#include "force_engine.hpp"



DirectForceGPU::DirectForceGPU(ParticleSet& particles) : particles(particles) {
    // prepare gpu and buffers and all
    checkDevices();
    Message("Preparing GPU...");
   
    
    // 1. Get platform and device
    err = clGetPlatformIDs(1, &platform, &numPlatforms);
    checkError(err, "clGetPlatformIDs");

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, nullptr);
    checkError(err, "clGetDeviceIDs");

    // 2. Create Context and queue
    context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &err);
    checkError(err, "clCreateContext");

    
    cl_queue_properties properties[] = {0}; // Default properties
    queue = clCreateCommandQueueWithProperties(context, device, properties, &err);
    checkError(err, "clCreateCommandQueueWithProperties");

    // 3. Create buffers
    dataSize = particles.size() * sizeof(float);

    bufferX = clCreateBuffer(context, CL_MEM_READ_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer X");

    bufferY = clCreateBuffer(context, CL_MEM_READ_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer Y");

    bufferZ = clCreateBuffer(context, CL_MEM_READ_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer Z");

    bufferMass = clCreateBuffer(context, CL_MEM_READ_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer Mass");

    bufferForceX = clCreateBuffer(context, CL_MEM_WRITE_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer ForceX");

    bufferForceY = clCreateBuffer(context, CL_MEM_WRITE_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer ForceY");

    bufferForceZ = clCreateBuffer(context, CL_MEM_WRITE_ONLY, dataSize, nullptr, &err);
    checkError(err, "clCreateBuffer ForceZ");

    // 4. Load kernel
    kernelCode = loadKernel("src/kernels/direct_force_gpu.cl");
    const char *kernelSource = kernelCode.c_str();
    program = clCreateProgramWithSource(context, 1, &kernelSource, nullptr, &err);
    checkError(err, "clCreateProgramWithSource");

    size_t log_size;
    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

    err = clBuildProgram(program, 1, &device, nullptr, nullptr, nullptr); // compile kernel code
    if (err != CL_SUCCESS) {
        size_t logSize;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &logSize);
        std::vector<char> buildLog(logSize);
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logSize, buildLog.data(), nullptr);
        std::cerr << "Build log:\n" << buildLog.data() << std::endl;
        checkError(err, "clBuildProgram");
    }

    // 5. Create kernel
    kernel = clCreateKernel(program, "computeForces", &err);
    checkError(err, "clCreateKernel");

    // 6. Set kernel arguments
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &bufferX);
    checkError(err, "clSetKernelArg X");

    err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &bufferY);
    checkError(err, "clSetKernelArg Y");

    err = clSetKernelArg(kernel, 2, sizeof(cl_mem), &bufferZ);
    checkError(err, "clSetKernelArg Z");

    err = clSetKernelArg(kernel, 3, sizeof(cl_mem), &bufferMass);
    checkError(err, "clSetKernelArg Mass");

    err = clSetKernelArg(kernel, 4, sizeof(cl_mem), &bufferForceX);
    checkError(err, "clSetKernelArg ForceX");

    err = clSetKernelArg(kernel, 5, sizeof(cl_mem), &bufferForceY);
    checkError(err, "clSetKernelArg ForceY");

    err = clSetKernelArg(kernel, 6, sizeof(cl_mem), &bufferForceZ);
    checkError(err, "clSetKernelArg ForceZ");

    // additional params
    int numParticles = particles.size();
    err = clSetKernelArg(kernel, 7, sizeof(int), &numParticles);
    checkError(err, "clSetKernelArg numParticles");

    float softening = ForceEngine::softening;
    err = clSetKernelArg(kernel, 8, sizeof(float), &softening);  // 8 is the next argument index
    checkError(err, "clSetKernelArg softening");

    // 7. Set global size
    globalSize = particles.size();
    Message("GPU ready","#");    
}

DirectForceGPU::~DirectForceGPU() {
    clReleaseMemObject(bufferX);
    clReleaseMemObject(bufferY);
    clReleaseMemObject(bufferZ);
    clReleaseMemObject(bufferMass);
    clReleaseMemObject(bufferForceX);
    clReleaseMemObject(bufferForceY);
    clReleaseMemObject(bufferForceZ);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);
    Message("GPU released","#");
}


/**
 * ----------------------
 * !-- COMPUTE FORCES --!
 * ----------------------
 */

std::vector<Eigen::Vector3d> DirectForceGPU::computeForces() {
    setBuffers();

    err = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &globalSize, nullptr, 0, nullptr, nullptr);
    checkError(err, "clEnqueueNDRangeKernel");

    // Read results
    std::vector<float> forceX(particles.size()), forceY(particles.size()), forceZ(particles.size());
    err = clEnqueueReadBuffer(queue, bufferForceX, CL_TRUE, 0, dataSize, forceX.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueReadBuffer ForceX");
    err = clEnqueueReadBuffer(queue, bufferForceY, CL_TRUE, 0, dataSize, forceY.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueReadBuffer ForceY");
    err = clEnqueueReadBuffer(queue, bufferForceZ, CL_TRUE, 0, dataSize, forceZ.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueReadBuffer ForceZ");

    // translate results
    std::vector<Eigen::Vector3d> forces(particles.size());
    for (int i = 0; i < particles.size(); i++) {
        forces[i] = Eigen::Vector3d(forceX[i], forceY[i], forceZ[i]);
    }

    return forces;
}



void DirectForceGPU::setBuffers() {
    // prepare data as flat arrays
    std::vector<float> posX(particles.size()), posY(particles.size()), posZ(particles.size()), mass(particles.size());
    for (int i = 0; i < particles.size(); i++) {
        posX[i] = particles.get(i).position(0);
        posY[i] = particles.get(i).position(1);
        posZ[i] = particles.get(i).position(2);
        mass[i] = particles.get(i).mass;
    }

    std::cout << "posX: " << posX[0] << "," << posX[1] << std::endl;

    // write data to buffers
    err = clEnqueueWriteBuffer(queue, bufferX, CL_TRUE, 0, dataSize, posX.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer X");
    
    err = clEnqueueWriteBuffer(queue, bufferY, CL_TRUE, 0, dataSize, posY.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer Y");

    err = clEnqueueWriteBuffer(queue, bufferZ, CL_TRUE, 0, dataSize, posZ.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer Z");

    err = clEnqueueWriteBuffer(queue, bufferMass, CL_TRUE, 0, dataSize, mass.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer Mass");

    // Initialize force buffers to zero // might be useless
    std::vector<float> zeroForceX(particles.size(), 0.0f);
    std::vector<float> zeroForceY(particles.size(), 0.0f);
    std::vector<float> zeroForceZ(particles.size(), 0.0f);

    err = clEnqueueWriteBuffer(queue, bufferForceX, CL_TRUE, 0, dataSize, zeroForceX.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer ForceX");

    err = clEnqueueWriteBuffer(queue, bufferForceY, CL_TRUE, 0, dataSize, zeroForceY.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer ForceY");

    err = clEnqueueWriteBuffer(queue, bufferForceZ, CL_TRUE, 0, dataSize, zeroForceZ.data(), 0, nullptr, nullptr);
    checkError(err, "clEnqueueWriteBuffer ForceZ");


    //Message("Buffers set","#");
}





/**
 * ------------------------
 * !-- STATIC FUNCTIONS --!
 * ------------------------ 
 */

void DirectForceGPU::checkError(cl_int err, std::string operation) {
    if (err != CL_SUCCESS) {
        Message("Error during [" + operation + "]: " + std::to_string(err), "!");
        exit(EXIT_FAILURE);
    }
}

void DirectForceGPU::checkDevices() {
    Message("Checking GPU...");
    cl_uint numPlatforms;
    clGetPlatformIDs(1, nullptr, &numPlatforms);
    Message::print(" > Number of platforms: " + std::to_string(numPlatforms));
    //std::cout << "Number of platforms: " << numPlatforms << std::endl; // 1

    std::vector<cl_platform_id> platforms(numPlatforms);
    clGetPlatformIDs(numPlatforms, platforms.data(), nullptr);

    for (cl_uint i = 0; i < numPlatforms; i++) {
        cl_platform_id platform = platforms[i];
        cl_uint numDevices;
        clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, 0, nullptr, &numDevices);
        Message::print(" > Number of devices: " + std::to_string(numDevices));

        std::vector<cl_device_id> devices(numDevices);
        clGetDeviceIDs(platform, CL_DEVICE_TYPE_ALL, numDevices, devices.data(), nullptr);

        for (cl_uint j = 0; j < numDevices; j++) {
            cl_device_id device = devices[j];
            size_t nameSize;
            clGetDeviceInfo(device, CL_DEVICE_NAME, 0, nullptr, &nameSize);
            std::vector<char> name(nameSize);
            clGetDeviceInfo(device, CL_DEVICE_NAME, nameSize, name.data(), nullptr);
            Message::print(" > Device name: " + std::string(name.data()));
        }
    }
}


std::string DirectForceGPU::loadKernel(std::string filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file '" << filename << "'" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string kernelCode((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    return kernelCode;
}





