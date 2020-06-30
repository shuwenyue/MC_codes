#ifndef ENERGY_H
#define ENERGY_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "VectorMath.h"
#include "particle.h"
#include "box.h"

using namespace std;

class energy
{
	public:
		energy(int NParticles, double BoxLength, double Esum, double sigma, double epsilon, double cutoff);
		double getTotalEnergy(vector<double3> coord);
		double getIndividualEnergy(vector<double3> coord, int randparticle);
	protected:
			
	private:
		double E_;
		double r_;
		double BoxLength_;
		int NParticles_;
		double Esum_;
		double sigma_;
		double epsilon_;
		double xlength_;
		double ylength_;
		double zlength_;
		double cutoff_;
};

#endif
