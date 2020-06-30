#include <iostream>
#include <vector>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

box::box(double BoxLength)
{
	BoxLength_ = BoxLength;
} 

double box::boundaries(double length)
{
	//particle particleObj;
    	while (length > (BoxLength_/2))
        {
                length = length - BoxLength_;
		//particleObj.backInBox();
        }

        while (length < (-BoxLength_/2))
        {
                length = length + BoxLength_;
		//particleObj.backInBox();
        }

	return length;
}


