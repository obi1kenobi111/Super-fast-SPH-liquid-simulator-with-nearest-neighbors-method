#include "Solver.h"
#include <iostream>
#include <stdio.h> 
#include <time.h> 

#define pow7(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x))

void Solver::initialize()
{
	worldSize_width = 1.0f;
	worldSize_height = 1.0f;

	gravity.x = 0.0;
	gravity.y = -4.5;

	CreateGrid();
	//масса частицы  = это площадь €чейки * плотность
}

void Solver::update()
{
	//в другом солвере 0.07 и 0.05 дл€ 990 частиц
	//

	time_t a = clock();

	if (SolverIterations % 5 == 0)
	FindParticles(); //0.001c

	time_t b = clock();

	if (SolverIterations % 5 == 0)
	FindNeigbor();

	time_t c = clock();
	calculateDensityPressure3(); // на три кадра быстрее считает
	time_t d = clock();
	Solve3();//
	time_t e = clock();

	if (SolverIterations%5 == 0) 
	ClearCellsData(); //0.001c

	time_t f = clock();

	/*
	std::cout <<" ClearCellsData =  "<<((double)((double)(f - e) / CLOCKS_PER_SEC)) <<
		"  Solve3 = " <<((double)((double)(e - d) / CLOCKS_PER_SEC)) << 
		"  calcDencity&Pressure = " <<((double)((double)(d - c) / CLOCKS_PER_SEC)) <<
		"S3&&CDP = "<< ((double)((double)(e - c) / CLOCKS_PER_SEC))<<
		"  FindNeigbor = " <<((double)((double)(c - b) / CLOCKS_PER_SEC)) <<
		"  FindParticles =  " <<((double)((double)(b - a) / CLOCKS_PER_SEC)) << std::endl;
    */

	SolverIterations++;
	//std::cout << "Solver Iterations  = " << SolverIterations << std::endl;
}



void Solver::calculateDensityPressure3()
{
	int Counter = 0;
	int index;
	//float r2;


		for (int k = 0; k < cellsWithParticles.size(); k++) { 

		     index = cellsWithParticles[k].i * (N + 1) + cellsWithParticles[k].j;
			
			for (int i = 0; i < grid[index].particle_id.size(); i++) { 

				int p_i = grid[index].particle_id[i];//ADDED

				if (p_i == 56) {
					id_neigbors1.clear();
				}
				
				particles[p_i]->density = 0.0f;

				for (int j = 0; j < grid[index].particle_neigb_id.size(); j++) {

					int p_j = grid[index].particle_neigb_id[j]; //ADDED
					vec2f r = particles[p_i]->position
						- particles[p_j]->position;
		
					float r2 = r * r;

					if (r2 >= h2 || r2 < 1e-12)
					{
						continue;
					}
					
					if (p_i == 56) {
						id_neigbors1.push_back(p_j);
						//Counter++;
					}
					

					float q = sqrt(r2) / h;
					
					particles[p_i]->density += particles[p_j]->mass * MonaganFun(q);
					//std::cout << MonaganFun(q) << std::endl;

				}
				
				 particles[p_i]->pressure = pow7(particles[p_i]->density / REST_DENSITY)*B -  B; 
	
			}

			
			
		}
}

	
   






void Solver::Solve3()
{
	//float r2;
	int index;
	//float viscosity;

	for (int k = 0; k < cellsWithParticles.size(); k++)
	{ 
		
		    index = cellsWithParticles[k].i * (N + 1) + cellsWithParticles[k].j;

			for (int i = 0; i < grid[index].particle_id.size(); i++) {

				if (particles[grid[index].particle_id[i]]->pType == STATIC)
					continue;

				int p_i = grid[index].particle_id[i]; //ADDED

				for (int j = 0; j < grid[index].particle_neigb_id.size(); j++) {

					int p_j = grid[index].particle_neigb_id[j]; //ADDED помогло

					if(p_i == p_j)
						continue;

					vec2f r0 = particles[p_i]->position - particles[p_j]->position;
					float r2 = r0 * r0;

					if (r2 < h2 && r2 > 1e-12)
					{

						float q = sqrt(r2) / h;
					

						vec2f v0 = particles[p_i]->velocity - particles[p_j]->velocity;
#define alpha -10.f
	
						float viscosity = r0 * v0 * alpha * h / ((particles[p_j]->density + particles[p_i]->density) * (r2 + 0.01f * h2));

						
						particles[p_i]->velocity -= r0 * dt * particles[p_j]->mass *
							(particles[p_i]->pressure / (particles[p_i]->density * particles[p_i]->density)
								+ particles[p_j]->pressure / (particles[p_j]->density * particles[p_j]->density)
								+ viscosity) * DiffMonaganFun(q);
											


					}

				}

				particles[p_i]->velocity += gravity * dt; //расчет скорости
				particles[p_i]->position += particles[p_i]->velocity * dt; //расчет координат

				
				if (particles[p_i]->pType == STATIC) {
					continue;
				}
				else {
					BoundaryConditions(p_i);
				}
				

			}
	
	}

}

float Solver::MonaganFun(const float& q)
{
	float W = 0;
	 W = MonaganFunConst * (-3.0f * q * q * q * q + 8.0f * q * q * q - 6.0f * q * q + 1.0f);
	return W;
}

float Solver::DiffMonaganFun(const float& q) //
{
	float W = 0;
	W = MonaganDiffConst* (-q*q + 2.0f*q - 1.0f); //Ќ”∆≈Ќ Ћ»  ќ–≈Ќ№ HSQ делить
	//W =-q * q + 2.0f * q - 1.0f;
	//W *= MonaganDiffConst;
	return W;
}

float Solver::GaussFun(float q2)
{
	return (2.0f/(h2*PI))*(1.5 - q2)*exp(-q2);
}

float Solver::DiffGaussFun(float q2)
{
	return (2.0f / (h2 * PI))*(-2.0*exp(-q2)*(2.5 - q2))/h; //нцмно умножить на вектор r
}


void Solver::BoundaryConditions(const int &i)
{

		if (particles[i]->position.x >= worldSize_width)
		{
			particles[i]->velocity.x = particles[i]->velocity.x * BOUND_DAMPINING;
			particles[i]->position.x = worldSize_width - BOUNDARY;
		}

		if (particles[i]->position.x <= 0.f)
		{
			particles[i]->velocity.x = particles[i]->velocity.x * BOUND_DAMPINING;
			particles[i]->position.x = BOUNDARY;
		}

		if (particles[i]->position.y >= worldSize_height)
		{
			particles[i]->velocity.y = particles[i]->velocity.y * BOUND_DAMPINING;
			particles[i]->position.y = worldSize_height - BOUNDARY;
		}

		if (particles[i]->position.y <= 0.f)
		{
			particles[i]->velocity.y = particles[i]->velocity.y * BOUND_DAMPINING;
			particles[i]->position.y = BOUNDARY;
		}
}


void Solver::addParticle(vec2f pos, vec2f vel)
{
	if (currentParticles <= MAX_PARTICLES)
	{
		particles[currentParticles] = new Particle(pos, vel, currentParticles);
		particles[currentParticles]->mass = DEFAULT_MASS;
		currentParticles++;
	}
}

void Solver::addParticle(float mass, vec2f pos, vec2f vel)
{
	if (currentParticles <= MAX_PARTICLES)
	{
		particles[currentParticles] = new Particle(pos, vel, currentParticles);
		particles[currentParticles]->mass = mass;
		currentParticles++;
	}
}

void Solver::addWall(vec2f pos)
{
	if (currentParticles <= MAX_PARTICLES)
	{
		
			particles[currentParticles] = new Particle(pos, currentParticles);
			particles[currentParticles]->mass = DEFAULT_MASS * 4.5f; // * 2 // 2.5
			currentParticles++;
		
		
	}
}



void Solver::CreateGrid()
{

	for (int i = 0; i < N + 1; i++) { //with

		for (int j = 0; j < N + 1; j++) { //height

			vec2f cell_pos;
			cell_pos.x = i * worldSize_width / N; 
			cell_pos.y = j * worldSize_width / N;

			grid[i*(N+1) + j].position = cell_pos;
			// ќќ–ƒ»Ќј“џ „ј—“»÷ Ќ”∆Ќќ ƒќћЌќ∆ј“№ Ќј N
		}
	}
	
}

void Solver::FindParticles() 
{

	CellWithParticles Cells_id;
	int index;

	for (int k = 0; k < currentParticles; k++) { 

		int i = (int)(std::round(particles[k]->position.x * N));
		int j = (int)(std::round(particles[k]->position.y * N));
		index = i*(N + 1) + j;

		grid[index].particle_id.push_back(k);
		
		grid[index].NumParticles++;
		if (grid[index].NumParticles == 1) { 
			Cells_id.i = i;
			Cells_id.j = j;
			cellsWithParticles.push_back(Cells_id); 
		}
		
	}
	
}

void Solver::FindNeigbor()
{
	std::vector<int> id;
	int ii, jj;
	for (int k = 0; k < cellsWithParticles.size(); k++) {


		ii = cellsWithParticles[k].i;
		jj = cellsWithParticles[k].j; //индексы узла

		FindNeighbors(ii, jj, id); //ищем соседей к данному узлу сетки
		grid[ii * (N + 1) + jj].particle_neigb_id = id;
		id.clear();
	}

}

void Solver::FindNeighbors( const int &ii, const int &jj, std::vector<int>& id)
{
	id.insert(id.end(), grid[ii*(N+1) + jj].particle_id.begin(), grid[ii * (N + 1) + jj].particle_id.end()); // частицы из самой €чейки

	//индексы частиц из соседних €чеек

	if (ii == 0 || jj == 0 || ii == N || jj == N) { 
		
		
		if (ii==N&&(jj!=N&&jj!=0)) { 

			//побокам(добавлено по 2 строки)
			id.insert(id.end(), grid[ii * (N + 1) + jj - 1].particle_id.begin(), grid[ii * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[ii * (N + 1) + jj + 1].particle_id.begin(), grid[ii * (N + 1) + jj + 1].particle_id.end());

			//плюс еще три

			id.insert(id.end(), grid[(ii-1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii-1) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.end());
	    }
		else if (ii == 0 && (jj != N && jj != 0)) { 

			id.insert(id.end(), grid[ii * (N + 1) + jj - 1].particle_id.begin(), grid[ii * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[ii * (N + 1) + jj + 1].particle_id.begin(), grid[ii * (N + 1) + jj + 1].particle_id.end());

			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.end());
			
		}
		else if (jj == 0 && (ii != N && ii != 0)) {//нижн€€ грань
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());

			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.end());
			id.insert(id.end(), grid[(ii) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii) * (N + 1) + jj + 1].particle_id.end());
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.end());
		}
		else if (jj == N && (ii != N && ii != 0)) {//верхн€€ нижн€€ грань
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());

			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.end());
		}

		// отдельно узловые точки
		if (ii==0&&jj==0) { //добавлено сверху и снизу
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());

			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.end());
			id.insert(id.end(), grid[(ii) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii) * (N + 1) + jj + 1].particle_id.end());
		
		}
		
		if (ii==0&&jj==N) {
			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());

			id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii) * (N + 1) + jj - 1].particle_id.end());
		}

		if (ii==N&&jj==0) {
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());

			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.end());
			id.insert(id.end(), grid[(ii) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii) * (N + 1) + jj + 1].particle_id.end());
		}

		if (ii == N && jj == N) {
			id.insert(id.end(), grid[(ii) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii) * (N + 1) + jj - 1].particle_id.end());

			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.end());
			id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());
		}
	}
	else {
		// вадрат
		id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj + 1].particle_id.end());
		id.insert(id.end(), grid[(ii) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii) * (N + 1) + jj + 1].particle_id.end());
		id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj + 1].particle_id.end());
	
		id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj].particle_id.end());
		id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj].particle_id.end());

		id.insert(id.end(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii + 1) * (N + 1) + jj - 1].particle_id.end());
		id.insert(id.end(), grid[(ii) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii) * (N + 1) + jj - 1].particle_id.end());
		id.insert(id.end(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.begin(), grid[(ii - 1) * (N + 1) + jj - 1].particle_id.end());
	}


}

void Solver::ClearCellsData()
{
	
	int index;
	for (int i = 0; i < cellsWithParticles.size(); i++) { 

		index = cellsWithParticles[i].i * (N + 1) + cellsWithParticles[i].j;
		grid[index].NumParticles = 0;
		grid[index].particle_id.clear();
		grid[index].particle_neigb_id.clear();//added

	}

	cellsWithParticles.clear();
	
}

