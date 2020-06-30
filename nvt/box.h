#ifndef BOUNDARIES_H
#define BOUNDARIES_H

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
		box (double BoxLength);
		double boundaries(double length);

	protected:
			
	private:
		double BoxLength_;

};

#endif
