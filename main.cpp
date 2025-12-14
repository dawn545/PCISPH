#include"SPHsolver.h"
#include"renderer.h"
#include<GL/glut.h>
#include<ctime>
#include<iostream>

SPHsolver* sphSolver;
Renderer* renderer;

void display()
{
    renderer->render();
}

void Update()
{
    sphSolver->step();
    glutPostRedisplay();
}

void keyboard(unsigned char key,  __attribute__((unused)) int x, __attribute__((unused)) int y)
{
    switch(key)
    {
        case ' ':
            if(sphSolver->GetParticles().size() >= sphSolver->max_particles)
            {
                std::cout << "maximum particles reached!" << std::endl;
            }
            else
            {
                unsigned int placed =0;
                // Create particles in grid pattern
                for(float y=SPHsolver::view_height/1.5f-SPHsolver::view_height/5.f;y<SPHsolver::view_height/1.5f +SPHsolver::view_height/5.f;y +=sphSolver->h*0.95f)
                {
                    for(float x=SPHsolver::view_width/2.f -SPHsolver::view_width/5.f;x<SPHsolver::view_width/2.f +SPHsolver::view_width/5.f;x +=sphSolver->h*0.95f)
                    {
                        if(placed++ < sphSolver->block_particles && sphSolver->GetParticles().size()<sphSolver->max_particles)
                        {
                            sphSolver->GetParticles().push_back(Particle(x,y));
                        }
                    }
                }
            }
            break;
        case 'r':
        case 'R':
            sphSolver->GetParticles().clear();
            sphSolver->InitSPH();
            break;
    }
}

void InitGL()
{
    glClearColor(0.9f,0.9f,0.9f,1);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(sphSolver->h/2.f);
    glMatrixMode(GL_PROJECTION);
}

int main(int argc,char** argv)
{   
    sphSolver = new SPHsolver();
    renderer = new Renderer(*sphSolver);
    srand(time(NULL));
    glutInitWindowSize(sphSolver->window_width,sphSolver->window_height);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInit(&argc,argv);
    glutCreateWindow("SPH Fluid Simulation");
    glutDisplayFunc(display);
    glutIdleFunc(Update);
    glutKeyboardFunc(keyboard);

    InitGL();
    sphSolver->InitSPH();

    glutMainLoop();
    return 0;

}