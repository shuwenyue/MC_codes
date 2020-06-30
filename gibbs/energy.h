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

using namespace std;

class energy
{
	public:
		energy(double Esum, double sigma, double epsilon, double cutoff);
		double getTotalEnergy(int whichbox, vector<double4> coord, double BoxLength1, double BoxLength2, double NParticles1, double NParticles2);
		double getIndividualEnergy(int whichbox, vector<double4> coord, double BoxLength, double NParticles, int randparticle);
		double getECorrection(double BoxLength, double NParticles);
    //extern int var=0;
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
