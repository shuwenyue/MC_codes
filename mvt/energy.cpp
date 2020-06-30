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

energy::energy(double Esum, double sigma, double epsilon, double cutoff)
{
	Esum_ = Esum;
	sigma_ = sigma;
	epsilon_ = epsilon;
	cutoff_ = cutoff;
}

double energy::getTotalEnergy(vector<double3> coord, double BoxLength, double NParticles)
{
	Esum_ = 0;

	for (int particle1 = 0; particle1 < NParticles; particle1++)
    {
    	for (int particle2 = particle1+1; particle2 < NParticles; particle2++)
        {
        	xlength_ = coord[particle1].x - coord[particle2].x;
           	ylength_ = coord[particle1].y - coord[particle2].y;
           	zlength_ = coord[particle1].z - coord[particle2].z;

           	// periodic boundary conditions
           	box boxObj;
           	xlength_ = boxObj.boundaries(xlength_, BoxLength);
           	ylength_ = boxObj.boundaries(ylength_, BoxLength);
           	zlength_ = boxObj.boundaries(zlength_, BoxLength);

       		// calculate r
        	r_ = sqrt( (xlength_*xlength_) + (ylength_*ylength_) + (zlength_*zlength_) );
			particle particleObj;

       		// LJ formula
       		if (r_ < cutoff_ )
        	{
         		//cout << "r = " << r_ << endl;
            	E_ = 4 * epsilon_ * ( pow((sigma_/r_),12) - pow((sigma_/r_),6));
            	Esum_ = Esum_ + E_;
    		}

		}
    }
	return Esum_;
}



double energy::getIndividualEnergy(vector<double3> coord, double BoxLength, int NParticles, int randparticle)
{
	Esum_ = 0;
	for (int particle = 0; particle < NParticles; particle++)
    {
    	if (randparticle != particle)
        {
        	xlength_ = coord[particle].x - coord[randparticle].x;
            ylength_ = coord[particle].y - coord[randparticle].y;
            zlength_ = coord[particle].z - coord[randparticle].z;

            // periodic boundary conditions
            box boxObj;
            xlength_ = boxObj.boundaries(xlength_, BoxLength);
            ylength_ = boxObj.boundaries(ylength_, BoxLength);
            zlength_ = boxObj.boundaries(zlength_, BoxLength);
	
            // calculate r
            r_ = sqrt( (xlength_*xlength_) + (ylength_*ylength_) + (zlength_*zlength_) );

           	// LJ formula
            if (r_ < cutoff_ )
            {
            	E_ = 4 * epsilon_ * ( pow((sigma_/r_),12) - pow((sigma_/r_),6));
				Esum_ = Esum_ + E_;
                
                if ((Esum_*Esum_) < 0.00000000000001)
                {
                    cout << "E problem " << endl;
                    cout << "r: " << r_ << endl;
                    cout << "Esum: " << Esum_ << endl;
                    cout << "E_ : " << E_ << endl;
                    cout << "randparticle: " << randparticle << endl;
                    
                    exit(0);
                }
            }
            

            
        }
	}
    return Esum_;
}

double energy::getECorrection(double BoxLength, double NParticles)
{
	double Ecorr = ((double) 8/(double) 3)*3.14159265*(NParticles*NParticles/(BoxLength*BoxLength*BoxLength))*epsilon_*sigma_*sigma_*sigma_*(((double)1/(double)3)*pow((sigma_/cutoff_),9)-pow((sigma_/cutoff_),3));

	return Ecorr;
}
