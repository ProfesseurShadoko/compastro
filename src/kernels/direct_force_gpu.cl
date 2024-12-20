
// Simulated atomic add for floats in OpenCL 1.2
float atomic_add_float(volatile __global float *ptr, float value) {
    union {
        unsigned int i;
        float f;
    } old, newVal;

    do {
        old.f = *ptr; // Read the current value as a float
        newVal.f = old.f + value; // Add the new value
    } while (atomic_cmpxchg((volatile __global unsigned int *)ptr, old.i, newVal.i) != old.i);

    return old.f; // Return the old value (optional)
}


__kernel void computeForces(
    __global const float *posX,  // X positions of particles
    __global const float *posY,  // Y positions of particles
    __global const float *posZ,  // Z positions of particles
    __global const float *mass,  // Mass of particles
    __global float *forceX,      // Output X forces
    __global float *forceY,      // Output Y forces
    __global float *forceZ,      // Output Z forces
    const int numParticles,      // Total number of particles
    const float softening        // Softening parameter
) {
    // Get global IDs
    int i = get_global_id(0);  // Particle i
    int j = get_global_id(1);  // Particle j

    // Check bounds
    if (i >= numParticles || j >= numParticles || i == j) return;

    // Load particle data
    float xi = posX[i], yi = posY[i], zi = posZ[i];
    float xj = posX[j], yj = posY[j], zj = posZ[j];
    float mj = mass[j];
    float mi = mass[i];

    // Compute softened distance
    float dx = xj - xi;
    float dy = yj - yi;
    float dz = zj - zi;
    float r2 = dx * dx + dy * dy + dz * dz + softening * softening; // Add softening squared
    float invR = rsqrt(r2);         // 1 / sqrt(r2)
    float invR3 = invR * invR * invR;

    // Compute force magnitude
    float f = mi * mj * invR3;

    // Compute force components
    float fx = f * dx;
    float fy = f * dy;
    float fz = f * dz;

    // Atomic addition to avoid race conditions
    atomic_add_float(&forceX[i], fx);
    atomic_add_float(&forceY[i], fy);
    atomic_add_float(&forceZ[i], fz);
}
