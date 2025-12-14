#include "renderer.h"
#include<GL/glut.h>

Renderer::Renderer(SPHsolver& s):solver(s){}

void Renderer::render()
{
    
	glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glOrtho(0, SPHsolver::view_width, 0, SPHsolver::view_height, 0, 1);

	glColor4f(0.2f, 0.6f, 1.f, 1);
	glBegin(GL_POINTS);
	for (auto &p : solver.GetParticles())
	{
		glVertex2f(p.x(0), p.x(1));
	}
	glEnd();

	glutSwapBuffers();
 
}