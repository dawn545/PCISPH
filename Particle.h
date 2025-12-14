#pragma once
#include<eigen3/Eigen/Dense>
using namespace Eigen;

class Particle 
{
    public:
        Vector2d x;           // 当前位置
        Vector2d v;           // 当前速度
        Vector2d f;           // 压力力
        Vector2d x_pred;      // 预测位置
        Vector2d v_pred;      // 预测速度
        Vector2d f_nonpressure; // 非压力力（重力+粘性）

        float rho;            // 当前密度
        float rho_pred;       // 预测密度
        float rho_error;      // 密度误差
        float p;              // 压力
        
        Particle(float x_, float y_)
            : x(x_, y_), v(0.f, 0.f), f(0.f, 0.f), 
              x_pred(x_, y_), v_pred(0.f, 0.f),
              rho(0.f), rho_pred(0.f), rho_error(0.f), p(0.f) {}
};
