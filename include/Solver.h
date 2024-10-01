#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "helper.h"

#define MAX_PARTICLES 10000

// Constants
#define h 0.025 // 0.04 (default)	 0.06		
#define DEFAULT_MASS  0.3 // 0.05 (default) 0.01 (ћ≈ЌяЋќ—№)/ 0.3
#define REST_DENSITY 1000.0	// 1000.0 (default)
#define h2 h*h		
#define dt 0.001	// 0.003 (default)
#define PI 3.1415

// Simulation
#define BOUNDARY 0.005		// 0.005 (default) (ћен€лось) 0.0005		
#define BOUND_DAMPINING -0.8	//

// Grid 
#define GRID_XY(X,Y) (Y*GRID_WIDTH) + X
#define GRID_WIDTH 100				
#define GRID_SIZE GRID_WIDTH * GRID_WIDTH 	// Number of cells
#define GRID_TRANSFORM 1.0/GRID_WIDTH		// This is a basis, so multiply by cell number to get the location of the cell

//#define B  10.0; //газова€ посто€нна€

//added
#define radius 0.5/GRID_WIDTH //радиус частиц, которые расположены на рассто€нии пор€дка радиуса

#define mass_const REST_DENSITY*radius*radius*radius*8.0;

#define B 10.0 //газова€ посто€нна€

////////////
#define MonaganFunConst 5.0f/(PI * h2)
#define MonaganDiffConst 4.0f*3.0f*5.0f/(PI* h2*h)


enum particleType
{
	DYNAMIC,
	STATIC
};

struct Cell { //структура узла
	std::vector<int> particle_id;
	std::vector<int> particle_neigb_id;
	vec2f position;
	int NumParticles = 0;
};

struct CellWithParticles { //NEW
	int i;
	int j;
};

struct Particle
{
	// Use for dynamic particles
	Particle(vec2f position, vec2f velocity, int id)
	{
		identifier = id;

		pType = DYNAMIC;

		mass = DEFAULT_MASS;
		density = REST_DENSITY;
		pressure = 0.f;

		this->position = position;
		this->velocity = velocity;
	}

	// Use for static 
	Particle(vec2f position, int id)
	{
		identifier = id;

		pType = STATIC;

		mass = DEFAULT_MASS;
		density = REST_DENSITY;
		pressure = 0.f;

		this->position = position;
	}

	unsigned int identifier;

	// Dynamic or static
	particleType pType;

	float mass;
	float density;
	float pressure;

	vec2f position;
	vec2f velocity;

};


class Solver
{
public:
	Solver( ) {


	};

	~Solver() {};

	void initialize();
	void update();
	void addParticle(vec2f pos, vec2f vel = vec2f(0.0f, 0.0f));
	void addParticle(float mass, vec2f pos, vec2f vel = vec2f(0.0f, 0.0f));

	void addWall(vec2f pos);


	Particle* particles[MAX_PARTICLES]; 

	float worldSize_width;
	float worldSize_height;

	int currentParticles = 0;

private:

	void BoundaryConditions(const int &i);
	
	
	float MonaganFun(const float& q);
	float DiffMonaganFun(const float& q);
	float GaussFun(float q2);
	float DiffGaussFun(float q2);



	vec2f gravity;
	void CreateGrid();
	
	void FindNeighbors(const int &ii, const int &jj, std::vector<int> &id);//поиск ближайших соседей, получаем индексы €чейки частицы

public:
	void FindParticles();
	void FindNeigbor();
	void ClearCellsData();
	//бежим по частицам
	void calculateDensityPressure3();//бежим по €чейкам
	void Solve3();

	static const int N = 40; //число €чеек дл€ квадратной решетки // 25 // 40
	//N = 11 h  = 0.09 N = 1/h; m = 0.1;
	//N = 25 h  = 0.045 m  = 0.35
	
public:
	Cell grid[(N+1)*(N+1)]; // размерность -  число узлов сетки

	//ƒл€ отладки
	std::vector<int> id_neigbors1;

	//NEW
	std::vector<CellWithParticles> cellsWithParticles;

	//NEW
	float dr; //рассто€ние между частицами
	int SolverIterations = 0;

	
};
#pragma once
