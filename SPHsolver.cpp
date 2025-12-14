#include "SPHsolver.h"
#include <cmath>
#include <iostream>


int SPHsolver::window_width = 800;
int SPHsolver::window_height = 600;
double SPHsolver::view_width = 1.5 * SPHsolver::window_width;
double SPHsolver::view_height = 1.5 * SPHsolver::window_height;

SPHsolver::SPHsolver()
{
    h = 16.f;
    rest_dens = 1000.f;
    visc = 200.f;
    mass = 2.5f;
    dt = 0.005f;
    hsq = h * h;
    sigma = 00000.f;
    bound_damp = -0.5f;
    max_particles = 2000;
    dam_particles = 1000;
    block_particles = 250;
    poly6 = 4.f / (M_PI * pow(h,8.f));
    spiky_grad = -10.f / (M_PI * pow(h,5.f));
    visc_lap = 40.f / (M_PI * pow(h,5.f));
    eps = h;
    delta = 1.f / (rest_dens * rest_dens * dt * dt);
}

void SPHsolver::InitSPH()
{
    std::cout << "Initializing Dam Break with " << dam_particles << " particles." << std::endl;
    for (float y = eps; y < view_height - eps * 2.f; y += h)
    {
        for (float x = view_width / 4; x <= view_width / 2; x += h)
        {
            if (particles.size() < dam_particles)
            {
                float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                particles.push_back(Particle(x + jitter, y));
            }
            else
            {
                return;
            }
        }
    }
}

void SPHsolver::step()
{
    ComputeNonPressureForces();
    InitPressureForces();
    InitPrediction();

    for(int iter =0,tag=1; iter <10 && tag!= 0  ; iter++)
    {
        PredictVelocityPosition();
        PredictDensity();
        UpdatePressure();
        ComputePressureForces();
        Check(tag);
    }
    Integrate();
}

void SPHsolver::InitPressureForces()
{
    for (auto& p: particles)
    {
        p.p = 0.f;
        p.f = Eigen::Vector2d(0.f, 0.f);
    }
}

void SPHsolver::InitPrediction()
{
    for (auto& p : particles)
    {
        p.v_pred = p.v + dt * p.f_nonpressure / mass;
        p.x_pred = p.x + dt * p.v_pred;
    }
}


void SPHsolver::ComputeNonPressureForces()
{   
    for (auto& pi : particles)// compute density first
    {
        pi.rho = 0.f;
        for (auto& pj : particles)
        {
            Eigen::Vector2d rij = pj.x - pi.x;
            float r2 = rij.squaredNorm();
            if (r2 < hsq)
            {
                pi.rho+= mass * poly6 * pow(hsq - r2, 3.f);
            }
        }
    }

    for (auto& pi : particles)
    {
        Eigen::Vector2d fvisc(0.f, 0.f);

        for (auto& pj : particles)
        {
            if (&pi == &pj) continue;
            Eigen::Vector2d rij = pj.x - pi.x;
            float r = rij.norm();
            if (r < h)
            {
                fvisc += visc * mass * (pj.v - pi.v) / pj.rho * visc_lap * (h - r);
            }
        }

        Eigen::Vector2d G(0.f, -10.f);
        //Eigen::Vector2d fgrav =  mass * G ;
        Eigen::Vector2d fgrav = mass * G /pi.rho;
        pi.f_nonpressure =  fvisc + fgrav;
    }
}

void SPHsolver::PredictVelocityPosition()
{
    for(auto& p :particles)
    {
        Eigen::Vector2d a =  p.f / mass;
        p.v_pred += dt * a;
        p.x_pred += dt * p.v_pred;

        // boundary conditions
        if (p.x_pred(0) - eps < 0.f) // left wall
        {
            p.v_pred(0) *= bound_damp;
            p.x_pred(0) = eps;
        }
        if (p.x_pred(0) + eps > view_width) // right wall
        {
            p.v_pred(0) *= bound_damp;
            p.x_pred(0) = view_width - eps;
        }
        if (p.x_pred(1) - eps < 0.f) // bottom wall
        {
            p.v_pred(1) *= bound_damp;
            p.x_pred(1) = eps;
        }
        if (p.x_pred(1) + eps > view_height) // top wall
        {
            p.v_pred(1) *= bound_damp;
            p.x_pred(1) = view_height - eps;
        }
    }
}

void SPHsolver::PredictDensity()
{
    for(auto& pi :particles)
    {
        pi.rho_pred = 0.f;
        for(auto& pj :particles)
        {
            Eigen::Vector2d rij = pj.x_pred - pi.x_pred;
            float r2 = rij.squaredNorm();

            if(r2 < hsq)
            {
                pi.rho_pred += mass * poly6 * pow(hsq - r2,3.f);
            }
        }
        pi.rho_error = pi.rho_pred - rest_dens;
    }
}

void SPHsolver::UpdatePressure()
{
    for(auto& p : particles)
    {      
        p.p += delta * p.rho_error;
    }
}

void SPHsolver::ComputePressureForces()
{
    for(auto& pi : particles)
    {
        Eigen::Vector2d fpress(0.f,0.f);

        for(auto& pj : particles)
        {
            if(&pi == &pj) continue;
            Eigen::Vector2d rij = pj.x_pred - pi.x_pred;
            float r = rij.norm();
            if(r < h && r > 1e-6f)
            {
              fpress += -rij.normalized() * /*mass **/ mass * (pi.p/(pi.rho_pred * pi.rho_pred) + pj.p/(pj.rho_pred * pj.rho_pred))  * spiky_grad * pow(h - r, 3.f);
            }
        }
        pi.f = fpress;
    }
}

int SPHsolver::Check(int tag)
{
    float max_error = 0.f;
    for(auto& p : particles)
    {
        max_error = std::max(max_error, std::abs(p.rho_error));
    }
    return (max_error < 0.01) ? 0 : 1;
}


void SPHsolver::Integrate()
{
    for (auto& p : particles)
    {
        p.v = p.v_pred;
        p.x = p.x_pred;

        // boundary conditions
        if (p.x(0) - eps < 0.f) // left wall
        {
            p.v(0) *= bound_damp;
            p.x(0) = eps;
        }
        if (p.x(0) + eps > view_width) // right wall
        {
            p.v(0) *= bound_damp;
            p.x(0) = view_width - eps;
        }
        if (p.x(1) - eps < 0.f) // bottom wall
        {
            p.v(1) *= bound_damp;
            p.x(1) = eps;
        }
        if (p.x(1) + eps > view_height) // top wall
        {
            p.v(1) *= bound_damp;
            p.x(1) = view_height - eps;
        }
    }
}
