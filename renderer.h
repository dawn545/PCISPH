#pragma once
#include "SPHsolver.h"

class Renderer
{
    public:
        Renderer(SPHsolver& solver);
        void render();
    private:
        SPHsolver& solver;
};