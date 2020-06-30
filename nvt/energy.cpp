#include <iostream>
#include <vector>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

energy::energy(int NParticles, double BoxLength, double Esum, double sigma, double epsilon, double cutoff)
{
	NParticles_ = NParticles;
	BoxLength_ = BoxLength;
	Esum_ = Esum;
	sigma_ = sigma;
	epsilon_ = epsilon;
	cutoff_ = cutoff;
}

double energy::getTotalEnergy(vector<double3> coord)
{
	Esum_ = 0;
	for (int particle1 = 0; particle1 < NParticles_; particle1++)
    {
    	for (int particle2 = particle1+1; particle2 < NParticles_; particle2++)
       	{
            xlength_ = coord[particle1].x - coord[particle2].x;
            ylength_ = coord[particle1].y - coord[particle2].y;
            zlength_ = coord[particle1].z - coord[particle2].z;

			// periodic boundary conditions
			box boxObj(BoxLength_);
			xlength_ = boxObj.boundaries(xlength_);
			ylength_ = boxObj.boundaries(ylength_);
			zlength_ = boxObj.boundaries(zlength_);

			// calculate r
			r_ = sqrt( (xlength_*xlength_) + (ylength_*ylength_) + (zlength_*zlength_) );

			// LJ formula
		    if (r_ < cutoff_ )
    		{
				//	cout << "r = " << r_ << endl;
        		E_ = 4 * epsilon_ * ( pow((sigma_/r_),12) - pow((sigma_/r_),6));
				Esum_ = Esum_ + E_;
       		}			
		}
	}	
	return Esum_;
}



double energy::getIndividualEnergy(vector<double3> coord, int randparticle)
{
	Esum_ = 0;
	for (int particle2 = 0; particle2 < NParticles_; particle2++)
    {
    	if (randparticle != particle2)
        {
        	xlength_ = coord[particle2].x - coord[randparticle].x;
            ylength_ = coord[particle2].y - coord[randparticle].y;
            zlength_ = coord[particle2].z - coord[randparticle].z;

            // periodic boundary conditions
			box boxObj(BoxLength_);
            xlength_ = boxObj.boundaries(xlength_);
            ylength_ = boxObj.boundaries(ylength_);
            zlength_ = boxObj.boundaries(zlength_);

            // calculate r
            r_ = sqrt( (xlength_*xlength_) + (ylength_*ylength_) + (zlength_*zlength_) );

            // LJ formula
            if (r_ < cutoff_ )
            {
            	E_ = 4 * epsilon_ * ( pow((sigma_/r_),12) - pow((sigma_/r_),6));
				Esum_ = Esum_ + E_;
            }
        }
	}
        return Esum_;
}

