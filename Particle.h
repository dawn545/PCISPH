#pragma once
#include<eigen3/Eigen/Dense>
using namespace Eigen;

class Particle 
{
    public:
        Vector2d x;//position
        Vector2d v;//velocity
        Vector2d f;//force
        Vector2d x_pred;
        Vector2d v_pred;
        Vector2d f_nonpressure;

        float rho_pred;
        float rho_error;
        float rho;//density
        float p;//pressure
        Particle(float x_,float y_):x(x_,y_),v(0.f,0.f),f(0.f,0.f),rho(0.f),p(0.f){}
};