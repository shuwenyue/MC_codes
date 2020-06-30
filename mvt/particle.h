#ifndef PARTICLE_H
#define PARTICLE_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include "VectorMath.h"
#include "energy.h"
#include "box.h"
#include "histogram.h"

using namespace std;

class particle
{
	public:
		void setPosition(double BoxLength, int NParticles);
		void moveParticle (double a, double b, double c, double rmax, int randparticle);
		void restoreParticleMove (double3 OldCoord, int randparticle);
		double3 getPosition(int randparticle);
		vector<double3> getVector();
		void addParticle(double randx, double randy, double randz);
		void deleteParticle(int randparticle, int NParticles);
	protected:
			
	private:
		double3 position;
		vector<double3> positionVector;
};

#endif
