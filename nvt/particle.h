#ifndef PARTICLE_H
#define PARTICLE_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <vector>
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

class particle
{
	public:
		particle(int NParticles, double BoxLength);
		void setPosition();
		vector<double3> moveParticle (double a, double b, double c, double rmax, int randparticle);
		vector<double3> restoreOld (double3 OldCoord, int randparticle);
		double3 getPosition(int randparticle);
		vector<double3> getVector();
	protected:
			
	private:
		int NParticles_;
		double BoxLength_;
		double3 position;
		vector<double3> positionVector;
};

#endif
