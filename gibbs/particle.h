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

using namespace std;

class particle
{
	public:
		void setPosition(double BoxLength1, double BoxLength2, int NParticles1, int NParticles2);
		void moveParticle (double a, double b, double c, double rmax, int randparticle);
		void restoreParticleMove (double4 OldCoord, int randparticle);
		double4 getPosition(int randparticle);
		vector<double4> getVector();
		void rescaleBox(int whichbox, double boxChange, int NTotalParticles);
		void restoreBox(vector<double4> OldVector, int NTotalParticles);
		void swapParticle(int randparticle, double randx, double randy, double randz, int newBox);
	protected:
			
	private:
		double4 position;
		vector<double4> positionVector;
};

#endif
