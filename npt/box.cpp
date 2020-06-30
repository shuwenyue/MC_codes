#include <iostream>
#include <vector>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

box::box(double BoxLength, double NParticles)
{
	BoxLength_ = BoxLength;
	NParticles_ = NParticles;
} 

double box::boundaries(double length)
{
    	while (length > (BoxLength_/2))
        {
                length = length - BoxLength_;
        }

        while (length < (-BoxLength_/2))
        {
                length = length + BoxLength_;
        }

	return length;
}

double box::rho(double BoxLength)
{
        double rho = NParticles_/(BoxLength*BoxLength*BoxLength);
	//cout << "BoxLength from rho function: " << BoxLength << endl;
        return rho;
}


