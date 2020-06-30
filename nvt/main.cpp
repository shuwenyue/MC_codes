#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <ctime>
#include "time.h"
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

int main () 
{
	//double Kb = 0.0013806; // (angstrom^2)*kg/(s^-2*K)
	//double epsilon = 0.165; // kg*angstroms/s^2
	//double sigma = 3.4; // angstroms
	//double T = 0.9*epsilon/Kb; // K
	
	double sigma = 1;
	double epsilon = 1;
	double T = 0.85;
    	double rho = 0.776;
    	int NParticles = 512;
	//double BoxLength = 7;
    	double BoxLength = pow(((double)NParticles/rho), 1.0/3.0);
	cout << "boxlength: " << BoxLength << endl;
	//int NParticles = (int)(rho*BoxLength*BoxLength*BoxLength);
    	//exit(0);
	double cutoff = 3;
	double rmax = 0.1;

	int TrialNumber = 1000000;
	double ETotal, E1, E2, deltaE, Ecorr;
	int acceptance = 0;
	int rejection = 0; 

	clock_t start, end;

	// Initialize positions
	particle particleObj(NParticles, BoxLength);
	particleObj.setPosition();

	// calculate inital total energy
	double Esum = 0;
	energy energyObj(NParticles, BoxLength, Esum, sigma, epsilon, cutoff);

	ETotal = energyObj.getTotalEnergy(particleObj.getVector());

	Ecorr = ((double) 8/(double) 3)*3.14159265*(NParticles/(BoxLength*BoxLength*BoxLength))*epsilon*sigma*sigma*sigma*(((double)1/(double)3)*pow((sigma/cutoff),9)-pow((sigma/cutoff),3));
	
	cout << "Ecorr: " << Ecorr << endl;
        //exit(0);
	// seed random number generator
	srand(time(0));

	// print to output file
    ofstream outFile;
    outFile.open("MDcompare.txt");


	start = clock();

	// begin trial loop
    for (int trial = 1; trial < TrialNumber; trial++)
    {
		//cout << "Trial: " << trial << endl;
		//cout << "Initial Total Energy: " << ETotal+Ecorr << endl;

        // choose a random particle     
        int randparticle = rand()%((NParticles-1));
		//cout << "randparticle: " << randparticle << endl;

		// save old coordinates 
        double3 OldCoord = particleObj.getPosition(randparticle);

        // calculate energy of random particle
        E1 = 0;
		E1 = energyObj.getIndividualEnergy(particleObj.getVector(),randparticle);
		//cout << "E1= " << E1 << endl;

        // random number between 0 and 1
        double a = ((double) rand() / (RAND_MAX)) ;
        double b = ((double) rand() / (RAND_MAX)) ;
        double c = ((double) rand() / (RAND_MAX)) ;

		// adjust position to new coordinates
		particleObj.moveParticle(a,b,c,rmax,randparticle);

		// calculate new energy
		E2 = 0;
		E2 = energyObj.getIndividualEnergy(particleObj.getVector(),randparticle);
		//cout << "E2= " << E2 << endl;
		
		// calculate energy difference
		double deltaE = E2 - E1;
		//cout << "deltaE: " << deltaE << endl;

		// calculate acceptance criteria	
		double p = exp(-deltaE/T);
		
		//cout << "acceptance criteria: " << p << endl;

		// accept new values if criteria is met	
		if (((double) rand() / (RAND_MAX)) < p)
		{
			ETotal = ETotal + deltaE;
			acceptance = acceptance + 1;	
			//cout << "accepted: " << acceptance << endl;
		}
		else
		{
			//cout << "rejected: " << rejection << endl;
			particleObj.restoreOld(OldCoord, randparticle);
			rejection = rejection + 1;
		}

		if ((trial)%1000 == 0)
		{
			outFile << trial << "   " << ETotal << endl; 
			cout << trial << "   " << ETotal/NParticles << endl; 
		}

		//cout << "Etotal: " << ETotal+Ecorr << endl;

		// check
		//Esum = 0;
		//double ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector());
		//cout << "Etotal check: " << ETotalcheck+Ecorr << endl;

		//cout << "rmax = " << rmax << endl;
		//cout << "--------------------------" << endl;

		// optimze rmax every 100 loops
		if ((trial)%1000 == 0)
		{
            	if(trial>0) rmax = 2*(double(acceptance)/trial)*rmax;
		}

	}

	end = clock();

	//cout << "CPU time: " << ((float)(end - start)/CLOCKS_PER_SEC) << endl;

	outFile.close();
	
	return 0;
}





