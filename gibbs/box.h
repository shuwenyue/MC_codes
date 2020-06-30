#ifndef BOX_H
#define BOX_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>

using namespace std;

class box
{
	public:
		double boundaries(double length, double BoxLength, double NParticles);
	protected:
			
	private:
		double BoxLength_;
		double NParticles_;

};

#endif
