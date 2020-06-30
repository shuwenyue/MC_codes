#ifndef ENERGY_H
#define ENERGY_H

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
#include "particle.h"
#include "box.h"
#include "histogram.h"

using namespace std;

class energy
{
	public:
		energy(double Esum, double sigma, double epsilon, double cutoff);
		double getTotalEnergy(vector<double3> coord, double BoxLength, double NParticles);
		double getIndividualEnergy(vector<double3> coord, double BoxLength, int NParticles, int randparticle);
		double getECorrection(double BoxLength, double NParticles);
	protected:
			
	private:
		double E_;
		double r_;
		double Esum_;
		double sigma_;
		double epsilon_;
		double xlength_;
		double ylength_;
		double zlength_;
		double cutoff_;
};

#endif
