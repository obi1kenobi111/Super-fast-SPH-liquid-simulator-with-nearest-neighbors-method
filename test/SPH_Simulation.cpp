#define _CRT_SECURE_NO_WARNINGS

// Every _1_ second(s) we calculate fps
#define TIME_CHUNK 1

#include <ctime>

#include "stdlib.h"
#include "math.h"
#include <iostream>
#include <stdio.h> 
#include <time.h> 

//#include <GL\glut.h>
#include "glut.h"

#include "helper.h"
#include "Solver.h"

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT WINDOW_WIDTH

char* window_title;

float window_width = WINDOW_WIDTH;
float window_height = WINDOW_HEIGHT;


int simPause = 1;
int drawCells = 1;

vec2f mousePos;
std::shared_ptr<Solver> sphSolver;

/// FPS ///

std::clock_t last_time = clock();
unsigned int frames_rendered = 0;

void initialFormation(int whSize)
{
	int particleCtr = 0;

	/*
	for (int x = 0; x < whSize; x++)
		for (int y = 0; y < whSize; y++)
		{
			sphSolver->addParticle(vec2f(x*1.f * GRID_TRANSFORM, (y*1.f + GRID_WIDTH / 2.f) * GRID_TRANSFORM));
			particleCtr++;
		}

	for (int x = 0; x < whSize / 2; x++) //квадрат поменьше
		for (int y = 0; y < whSize / 2; y++)
		{
			sphSolver->addParticle(vec2f((x*1.f + 0.5f) * GRID_TRANSFORM, y*1.f * GRID_TRANSFORM));
			particleCtr++;
		}
	*/

	/*
	for (int x = 0; x < whSize; x++) { // whSize/2
		for (int y = 0; y < whSize/2; y++) // whSize
		{
			//sphSolver->addParticle(vec2f(x * 0.5f * GRID_TRANSFORM, y * 0.5f * GRID_TRANSFORM)); // h = 0.02 N =50 m = 0.05 вылетает с вязкостью // 30 частиц
			// 
			//sphSolver->addParticle(vec2f( x * 2.0f * GRID_TRANSFORM, y * 2.0f * GRID_TRANSFORM)); // N  = 20  h = 0.06 m = 0.6 // число частиц 40 whSize/2 whSize ( СУПЕР СЦЕНА)
			//sphSolver->addParticle(vec2f( x * 2.0f * GRID_TRANSFORM, y * 2.0f * GRID_TRANSFORM)); // N  = 20  h = 0.07 m = 1.3 // число частиц 45 whSize/2 whSize dt = 0.05 ( СУПЕР СЦЕНА 2)
			// sphSolver->addParticle(vec2f(x * 2.0f * GRID_TRANSFORM, y * 2.0f * GRID_TRANSFORM)); // N  = 25 h = 0.06 m = 1.3 // число частиц 45 whSize/2 whSize dt = 0.05 ( СУПЕР СЦЕНА 2) 40 FPS
			//sphSolver->addParticle(vec2f(x * 2.0f * radius, y * 2.0f * radius));// расположены на расстоянии двух радиусов
			sphSolver->addParticle(vec2f(x * 2.0f * radius, y * 2.0f * radius));
			//sphSolver->addParticle(vec2f( x * 2.0f * GRID_TRANSFORM, y * 2.0f * GRID_TRANSFORM)); // N  = 11 h = 0.06 m = 0.3
			//sphSolver->addParticle(vec2f(x * 1.0f * GRID_TRANSFORM,  (y * 1.0f) * GRID_TRANSFORM)); // 50 частиц  N  = 50 h = 0.06 m = 0.3 стоб (ОТЛАДКА ПОСЛЕДНЕГО АЛГОРИТМА)
			//sphSolver->addParticle(vec2f(x * 0.5f * GRID_TRANSFORM ,  (y * 0.5f) * GRID_TRANSFORM)); //N = 100 m = 0.06 h = 0.03
			//sphSolver->addParticle(vec2f(x * 0.7f * GRID_TRANSFORM, 0.05f + y * 0.7f * GRID_TRANSFORM)); //для 10 000  частиц
			particleCtr++;
		}
	}
	*/

	//Тут dt сделали меньше 0.001  +++ %5 пересчет соседей
	for (int x = 1; x < whSize / 2; x++) { // whSize/2
		for (int y = 1; y < whSize / 2; y++) // whSize
		{
			sphSolver->addParticle(vec2f(x * 2.0f * radius, y * 2.0f * radius)); //число частиц 50
			particleCtr++;
		}
	}

	//ГРАНИЧНЫЕ УСЛОВИЯ 50 частиц было
	for (int i = 0; i <= 100; i++) { // 50 = 1/2*radius
		for (int j = 0; j <= 100; j++) {

			if (j == 0 || i == 0 || i == 100 || j == 100) {
				sphSolver->addWall(vec2f(i * 2.0f * radius, j * 2.0f * radius)); //число частиц 50
				particleCtr++;
			}
		}
	}

	//h =0.05 m = 1.5 число частиц 50 N  = 22 dt = 0.03 radius = 1/100 Масса границ * 2.5
	/*
	for (int x = 1; x < whSize/2; x++) { // whSize/2
		for (int y = 1; y < whSize/2; y++) // whSize
		{

			sphSolver->addParticle(vec2f(0.25 + x * 2.0f * radius, y * 2.0f * radius)); //число частиц 50
			particleCtr++; 
		}
	}

	//ГРАНИЧНЫЕ УСЛОВИЯ 50 частиц было
	for (int i = 0; i <= 50; i++) { // 50 = 1/2*radius
		for (int j = 0; j <= 50; j++) {

			if (j == 0 || i == 0 || i == 50 || j == 50) {
				sphSolver->addWall(vec2f(i * 2.0f * radius, j * 2.0f * radius)); //число частиц 50
				particleCtr++;
			}
		}
	}
	*/
	std::cout << "Number of Particles: " << particleCtr << std::endl;
}

void init()
{
	sphSolver.reset(new Solver());
	sphSolver->initialize();


	initialFormation(60); //создаем частицы

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, sphSolver->worldSize_width, 0.0, sphSolver->worldSize_height);

	window_title = (char*)malloc(sizeof(char) * 50);

	glutIgnoreKeyRepeat(false);

	std::cout << std::endl << "\tKeybinds: " << std::endl;
	std::cout << "n	- Disable cell visualization" << std::endl;
	std::cout << "Space   - Resume/Pause" << std::endl;
	std::cout << "x	- Place wall particle" << std::endl;
	std::cout << "LMB	- Heavy particle" << std::endl;
	std::cout << "RMB	- Ring of particles" << std::endl;

	std::cout << "Simulation is paused - Press 'space' to unpause" << std::endl;
}

void render()
{
	
	/*
	if (drawCells)
		for (int i = 0; i <grid.size() - 1; i++) {
			for (int j = 0; j < grid.size() - 1; j++) {

				glColor3ub(grid[i][j].position.x* 255.f, grid[i][j].position.y * 255.f, 144.f);
				//glColor3ub(255.f/i*i, 255.f/j*j, 144.f);

				glBegin(GL_POLYGON);
				//glVertex2f(c->position.x * GRID_TRANSFORM, c->position.y * GRID_TRANSFORM);
				//glVertex2f((c->position.x + 1) * GRID_TRANSFORM, c->position.y * GRID_TRANSFORM);
				//glVertex2f((c->position.x + 1) * GRID_TRANSFORM, (c->position.y + 1) * GRID_TRANSFORM);
				//glVertex2f(c->position.x * GRID_TRANSFORM, (c->position.y + 1) * GRID_TRANSFORM);
				glVertex2f(grid[i][j].position.x, grid[i][j].position.y);
				glVertex2f(grid[i+ 1][j].position.x, grid[i+1][j].position.y);
				glVertex2f(grid[i + 1][j + 1].position.x, grid[i + 1][j + 1].position.y);
				glVertex2f(grid[i][j + 1].position.x, grid[i][j + 1].position.y);
				glEnd();

			}
		}
	*/

	//рисуем сетку
	if (drawCells){
		for (int i = 0; i < sphSolver->N + 1; i++) {
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < sphSolver->N + 1; j++) {

				glColor3f(0, 0.2, 0.8);
				glVertex2f(sphSolver->grid[i * (sphSolver->N + 1) + j].position.x, sphSolver->grid[i * (sphSolver->N + 1) + j].position.y);
			}
			glEnd();
		}

		for (int j = 0; j < sphSolver->N + 1; j++) {
			glBegin(GL_LINE_STRIP);
			for (int i = 0; i < sphSolver->N + 1; i++) {

				glColor3f(0, 0.2, 0.8);
				glVertex2f(sphSolver->grid[i * (sphSolver->N + 1) + j].position.x, sphSolver->grid[i * (sphSolver->N + 1) + j].position.y);
			}
			glEnd();
		}

		//glColor3f(0, 1, 1);
		//glPointSize(5);
		//glBegin(GL_POINTS);
		//glVertex2f(grid[11*(20+1) + 11].position.x, grid[11*(20+1) + 11].position.y);
		//glEnd();
	}
			
		
	

	for (int x = 0; x < sphSolver->currentParticles; x++)
	{
		// Draw
		glPointSize(8);
		glBegin(GL_POINTS);
		glColor3f(0.372549f, 0.623529f, 0.623529f);
		//glColor3f(0.92549f, 0.93529f, 0.623529f);

		if (sphSolver->particles[x]->pType == STATIC) {
			glColor3f(0.672549f, 0.623529f, 0.023529f);
		}
		glVertex2d(sphSolver->particles[x]->position.x, sphSolver->particles[x]->position.y);
		glEnd();
	}

	//отрисовка соседей  частицы
	
	
	for (int x = 0; x < sphSolver->id_neigbors1.size(); x++)
	{
		
		glPointSize(8);
		glBegin(GL_POINTS);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex2d(sphSolver->particles[sphSolver->id_neigbors1[x]]->position.x, sphSolver->particles[sphSolver->id_neigbors1[x]]->position.y);
		
		glEnd();
	}
	
}

void display()
{
	glViewport(0, 0, window_width, window_height);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glClearColor(1.0f, 1.0f, 1.0f, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.0f, 0.0f, 1.0f);

	if (simPause == 0) {
		
			sphSolver->update();
		
	}

	render();

	glutSwapBuffers();


	if ((clock() - last_time) / CLOCKS_PER_SEC >= TIME_CHUNK)
	{
		memset(window_title, 0, 50);
		sprintf(window_title, "SPH Solver - Hunter Werenskjold - FPS: %d", frames_rendered / TIME_CHUNK);
		glutSetWindowTitle(window_title);

		frames_rendered = 0;
		last_time = clock();
	}
	frames_rendered++;

}

void idle()
{
	glutPostRedisplay();
}


void reshape(int width, int height) {
	glutReshapeWindow(window_width, window_height);
}

vec2f lastWallParticle;

void process_keyboard(unsigned char key, int x, int y)
{
	if (key == ' ')
	{
		simPause = 1 - simPause;
	}

	if (key == 'x')
	{
		if (sqrtf(pow(mousePos.x - lastWallParticle.x, 2.f) + pow(mousePos.y - lastWallParticle.y, 2.f)) >= 30)
		{
			sphSolver->addWall(vec2f((mousePos.x / window_width) * sphSolver->worldSize_width, (((-1 * mousePos.y + window_height) / window_height)) * sphSolver->worldSize_height));
			lastWallParticle = mousePos;
		}
	}

	if (key == 'n')
	{
		drawCells = 1 - drawCells;
	}

}


void process_mouse_movement(int x, int y)
{
	mousePos.x = x;
	mousePos.y = y;
}


void process_mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		sphSolver->addParticle(DEFAULT_MASS * 2, vec2f((mousePos.x / window_width) * sphSolver->worldSize_width, (((-1 * mousePos.y + window_height) / window_height)) * sphSolver->worldSize_height));
	}

	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		for (float i = 0; i < 360; i = i + 40)
		{
			sphSolver->addParticle(vec2f((mousePos.x / window_width) * sphSolver->worldSize_width + 0.02 * unitVecFromDeg(i).x, (((-1 * mousePos.y + window_height) / window_height)) * sphSolver->worldSize_height + 0.02 * unitVecFromDeg(i).y));
		}
	}
}

int main(int argc, char** argv)
{
	
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("SPH Solver - H Werenskjold - FPS: ");

	init();
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(process_keyboard);
	glutPassiveMotionFunc(process_mouse_movement);
	glutMouseFunc(process_mouse);
	glutMainLoop();
	

	

	/*
	sphSolver = new Solver();
	sphSolver->initialize();
	initialFormation(50); //создаем частицы

	for (int i = 0; i < 100; i++) {

		
		sphSolver->FindParticles();
		

		
		sphSolver->calculateDensityPressure3();
		sphSolver->Solve3();
		//clock_t  start = clock();
		sphSolver->ClearCellsData();
		//clock_t  end = clock();

		
		//double seconds = (double)(end - start) / CLOCKS_PER_SEC;
		//std::cout << "seconds =  " << seconds << std::endl;
		std::cout << "seconds =  " << sphSolver->timefun<< std::endl;
		sphSolver->timefun = 0;
	}
	*/
	
	
	return 0;
}
