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

using namespace std;

energy::energy(double Esum, double sigma, double epsilon, double cutoff)
{
	Esum_ = Esum;
	sigma_ = sigma;
	epsilon_ = epsilon;
	cutoff_ = cutoff;
}

double energy::getTotalEnergy(int whichbox, vector<double4> coord, double BoxLength1, double BoxLength2, double NParticles1, double NParticles2)
{
	Esum_ = 0;
	particle particleObj;

	for (int particle1 = 0; particle1 < (NParticles1 + NParticles2); particle1++)
    {
    	for (int particle2 = particle1+1; particle2 < (NParticles1 + NParticles2); particle2++)
        {
			// check all coordinates to only calculate energies for particles in box x
			if (coord[particle1].n == whichbox && coord[particle2].n == whichbox)
			{
        		xlength_ = coord[particle1].x - coord[particle2].x;
           		ylength_ = coord[particle1].y - coord[particle2].y;
            	zlength_ = coord[particle1].z - coord[particle2].z;

            	// periodic boundary conditions
            	box boxObj;
				if (whichbox == 1)
				{
            		xlength_ = boxObj.boundaries(xlength_, BoxLength1, NParticles1);
            		ylength_ = boxObj.boundaries(ylength_, BoxLength1, NParticles1);
            		zlength_ = boxObj.boundaries(zlength_, BoxLength1, NParticles1);
				}
				else
				{
					xlength_ = boxObj.boundaries(xlength_, BoxLength2, NParticles2);
            		ylength_ = boxObj.boundaries(ylength_, BoxLength2, NParticles2);
            		zlength_ = boxObj.boundaries(zlength_, BoxLength2, NParticles2);
				}

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
    }
	
	return Esum_;
}



double energy::getIndividualEnergy(int whichbox, vector<double4> coord, double BoxLength, double NTotalParticles, int randparticle)
{
	Esum_ = 0;
	particle particleObj;



	//cout << "calculating individual energy" << endl;
	//cout << "whichbox: " << whichbox << endl;

	for (int particle = 0; particle < NTotalParticles; particle++)
    {
		if (coord[particle].n == whichbox)
		{
        	if (randparticle != particle)
            {
				//cout << "randparticle number: " << randparticle << "-" << particleObj.getPosition(randparticle).n << endl;
            	xlength_ = coord[particle].x - coord[randparticle].x;
                ylength_ = coord[particle].y - coord[randparticle].y;
                zlength_ = coord[particle].z - coord[randparticle].z;

                // periodic boundary conditions
                box boxObj;

            	xlength_ = boxObj.boundaries(xlength_, BoxLength, NTotalParticles);
            	ylength_ = boxObj.boundaries(ylength_, BoxLength, NTotalParticles);
            	zlength_ = boxObj.boundaries(zlength_, BoxLength, NTotalParticles);

                // calculate r
                r_ = sqrt( (xlength_*xlength_) + (ylength_*ylength_) + (zlength_*zlength_) );

                // LJ formula
                if (r_ < cutoff_ )
                {
                    E_ = 4 * epsilon_ * ( pow((sigma_/r_),12) - pow((sigma_/r_),6));
                    /*if (E_ != E_)
                    {
                        cout << "r: " << r_ << endl;
                        cout << "sigma: " << sigma_ << endl;
                        cout << "ep: " << epsilon_ << endl;
                        cout << "randparticl: " << randparticle << endl;
                        cout << "particle: " << particle << endl;
                        cout << "E: " << E_ << endl;
                        cout << "coord1x: " <<coord[particle].x << endl;
                        cout << "coord1y: " <<coord[particle].y << endl;
                        cout << "coord1z: " <<coord[particle].z << endl;
                        cout << "randx: " <<coord[randparticle].x << endl;
                        cout << "randy: " <<coord[randparticle].y << endl;
                        cout << "randz: " <<coord[randparticle].z << endl;

                        exit(0);
                    }*/

					Esum_ = Esum_ + E_;
                }
        	}
        }
	}
    
    /*cout << "whichbox: " << whichbox << endl;
    cout << "randparticle " << randparticle << endl;
    cout << "nparticlestotal: " << NTotalParticles << endl;
    cout << "rand-x: " << coord[randparticle].x << endl;
    cout << "rand-y: " << coord[randparticle].y << endl;
    cout << "rand-z: " << coord[randparticle].z << endl;
    cout << "esum: " << Esum_ << endl;
    cout << "-----------" << endl;*/
    return Esum_;
}

double energy::getECorrection(double BoxLength, double NParticles)
{
	double Ecorr = ((double) 8/(double) 3)*3.14159265*(NParticles*NParticles/(BoxLength*BoxLength*BoxLength))*epsilon_*sigma_*sigma_*sigma_*(((double)1/(double)3)*pow((sigma_/cutoff_),9)-pow((sigma_/cutoff_),3));

	return Ecorr;
}
