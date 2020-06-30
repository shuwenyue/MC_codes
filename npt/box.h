#ifndef BOX_H
#define BOX_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>

//using namespace std;

class box
{
	public:
		box (double BoxLength, double NParticles);
		double boundaries(double length);
		double rho(double BoxLength);
	protected:
			
	private:
		double BoxLength_;
		double NParticles_;

};

#endif
