#pragma once
#include <vector>
#include "Particle.h"

class SPHsolver
{
public:
    SPHsolver();
    void InitSPH();
    void step();

    std::vector<Particle>& GetParticles()  { return particles; }
    std::vector<Particle> particles;

    float h; // kernel radius
    float rest_dens; // rest density
    float gas_const; // gas constant
    float visc; // viscosity constant
    float mass; // particle mass
    float dt; // time step
    float hsq; // radius squared for optimization
    float sigma;//surface tension parameter
    float delta;
    int tag;

    // smoothing kernels  
    float poly6;
    float spiky_grad;
    float visc_lap;

    // simulation parameters 
    float eps;
    float bound_damp;
    static int window_width;
    static int window_height;
    static double view_width;
    static double view_height;

    // interaction parameters 
    int max_particles;
    int dam_particles;
    int block_particles;

    void ComputeNonPressureForces();
    void InitPressureForces();
    void InitPrediction();
    void PredictVelocityPosition();
    void PredictDensity();
    void UpdatePressure();
    void ComputePressureForces();
    int Check(int tag);
    void Integrate();
};
