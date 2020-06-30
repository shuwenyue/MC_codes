#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"
#include "histogram.h"

using namespace std;

double box::boundaries(double length, double BoxLength)
{
    	while (length > (BoxLength/2))
        {
                length = length - BoxLength;
        }

        while (length < (-BoxLength/2))
        {
                length = length + BoxLength;
        }

	return length;
}

double box::rho(double BoxLength, double NParticles)
{
        double rho = NParticles/(BoxLength*BoxLength*BoxLength);
        return rho;
}


