#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

int main () 
{
	//double Kb = 0.0013806; // (angstrom^2)*kg/(s^2*K)
	//double epsilon = 0.1654; // kg*angstroms/s^2
	//double sigma = 3.4; // angstroms
	//double T = (0.9)*epsilon/Kb; // K
	//double P = (7.36626e-3)*epsilon/(sigma*sigma*sigma); // kg/(angstrom*s^2)

	// reduced units
	double sigma = 1;
	double epsilon = 1;
	double T = 1.5; 
	double P = 0.1; 
	double rho = 0.5;
	int NParticles = 216; 
	double BoxLength = pow(((double)NParticles/rho), 1.0/3.0);
	double cutoff = 3;

	int TrialNumber = 500000;
	int acceptance = 0;
	int rejection = 0; 

	clock_t start, end;

	// particle displacement variables
	double ETotal, E1, E2, deltaE, Ecorr;
	double rmax = 1;
	int Rchange = 1;

	// volume move variables
	double ETotal1, ETotal2, deltaETotal, EcorrV, EcorrV2;
	double deltavmax = 0.1*BoxLength*BoxLength*BoxLength; 
	double rho2, oldV, newV, newBoxLength;
	int Vchange = 1;


	// Initialize positions
	particle particleObj(NParticles, BoxLength);
	particleObj.setPosition();

	// calculate inital total energy
	double Esum = 0;
	energy energyObj(NParticles, Esum, sigma, epsilon, cutoff);
	ETotal = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength);
	//&cout << "boxlength: " << BoxLength << endl;

	// Ecorr for particle movement
	rho = NParticles/(BoxLength*BoxLength*BoxLength);
	Ecorr = energyObj.getECorrection(BoxLength);  // only used for particle displacement

	//cout << "ecorr: " << Ecorr << endl;

	// seed random number generator
	srand(time(0));

	// print to output file
        ofstream outFile;
        outFile.open("assignment1.txt");

	start = clock();

	// begin trial loop
        for (int trial = 1; trial < TrialNumber; trial++)
        {
			//cout << "Trial: " << trial << endl;

		// chose particle displacement or volume change
		double whichmove = ((double) rand() / (RAND_MAX));

		//whichmove = 0; 

		box boxObj(BoxLength, NParticles);

        if ((trial)%1000 == 0)
		{
        	outFile << trial << "   " << (boxObj.rho(BoxLength)) << endl; 
        	//cout << "density: " << boxObj.rho(BoxLength) << endl;
		}

		if (whichmove < ((double) 1/((double) 1+NParticles)))
		{
			// perform volume change
			//cout << "Volume change: " << Vchange++ << endl;

			

			// initial total energy
			//cout << "old rho: " << boxObj.rho(BoxLength) << endl;
			ETotal1 = ETotal;
			double Ecorr1 = energyObj.getECorrection(BoxLength);
			//cout << "ETotal1 = " << ETotal + energyObj.getECorrection(BoxLength) << endl;
			//cout << "ecorr 1: " << Ecorr1 << endl;

			// save old volume and coordinates
       		oldV = BoxLength*BoxLength*BoxLength;
			vector<double3> OldPosition = particleObj.getVector();

			// random number to change V
			double d = ((double) rand() / (RAND_MAX)) ;
			
			// adjust volume
			//cout << "oldV=" << oldV << endl;
			//double lognewV = log(oldV) + (2*d-1) * deltavmax; 
			//newV = exp((double)lognewV);
			newV = oldV + (d-0.5)*deltavmax;
			//cout << "newV=" << newV << endl;
			newBoxLength = pow((double)newV, 1.0/3.0);

			// rescale particles in box
			particleObj.rescaleBox(newBoxLength/BoxLength);

			// calculate new total energy
			//cout << "new rho: " << boxObj.rho(newBoxLength) << endl;
			ETotal2 = energyObj.getTotalEnergy(particleObj.getVector(), newBoxLength);
			//cout << "ETotal2 = " << ETotal2 + energyObj.getECorrection(newBoxLength)  << endl;
			
			double Ecorr2 = energyObj.getECorrection(newBoxLength);
			//cout << "ecorr 2: " << Ecorr1 << endl;

			// acceptance criteria
			deltaETotal = (ETotal2 + Ecorr2) - (ETotal1 + Ecorr1);
			//cout << "deltaETotal = " << deltaETotal << endl;
			double beta = 1/T;
			double accV = exp(-beta*(deltaETotal)-beta*P*(newV-oldV)+(NParticles)*log(newV/oldV));
			//cout << "acceptance criteria: " << accV << endl;

			// accept new values if criteria is met 
            if (((double) rand() / (RAND_MAX)) < accV)
            {
                acceptance = acceptance + 1;
               // cout << "Vchange accepted: " << acceptance << endl;
				BoxLength = newBoxLength;
				//cout << "E from volume move: " << ETotal2 + energyObj.getECorrection(BoxLength) << endl;
				ETotal = ETotal2;
            }
            else
            {
				particleObj.restoreBox(OldPosition);
                rejection = rejection + 1;
                //cout << "Vchange rejected: " << rejection << endl;
				//cout << "E from volume move (rej): " << ETotal1 + energyObj.getECorrection(BoxLength) << endl;

				// ETotal check2
				//double ETotalcheck1 = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength);
                //cout << "ETotalcheck2 = " << ETotalcheck1 + energyObj.getECorrection(BoxLength)  << endl;

             }

			//cout << "deltavmax= " << deltavmax << endl;
            //cout << "--------------------------" << endl;

			// optimze deltavmax every 100 loops
            if ((trial)%1000 == 0)
            {
				//cout << "density: " << boxObj.rho(BoxLength) << endl;
				
            	if(trial>0) deltavmax = 2*(double(acceptance)/trial)*deltavmax;
            }
			
			//outFile << trial << "   " << ETotal + energyObj.getECorrection(BoxLength) << endl;
			//outFile << trial << "   " << (boxObj.rho(BoxLength)) << endl;

			
		}
		else
		{

			// perform particle displacement
			//&cout << "Particle displacement: " << Rchange++ << endl; 

			box boxObj(BoxLength,NParticles);

			//&cout << "Initial Total energy: " << ETotal + energyObj.getECorrection(BoxLength) << endl;

			//&cout << "boxlength: " << BoxLength << endl;
			// choose a random particle     
            int randparticle = rand()%NParticles;
            //&cout << "randparticle: " << randparticle << endl;

            // save old coordinates 
            double3 OldCoord = particleObj.getPosition(randparticle);

            // calculate energy of random particle
            E1 = 0;
            E1 = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, randparticle);
            //&cout << "E1= " << E1 << endl;

            // random number between 0 and 1
            double a = ((double) rand() / (RAND_MAX)) ;
            double b = ((double) rand() / (RAND_MAX)) ;
            double c = ((double) rand() / (RAND_MAX)) ;

            // adjust position to new coordinates
            particleObj.moveParticle(a,b,c,rmax,randparticle);

            // calculate new energy
            E2 = 0;
            E2 = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, randparticle);
            //&cout << "E2= " << E2 << endl;

            // calculate energy difference
            double deltaE = E2 - E1;
            //&cout << "deltaE: " << deltaE << endl;

            // calculate acceptance criteria        
            double acc = exp(-deltaE/T);
            //&cout << "acceptance criteria: " << acc << endl;

            // accept new values if criteria is met 
            if (((double) rand() / (RAND_MAX)) < acc)
            { 
				ETotal = ETotal + deltaE;
                acceptance = acceptance + 1;
                //&cout << "accepted: " << acceptance << endl;
            }
            else
            {
                particleObj.restoreOld(OldCoord, randparticle);
                rejection = rejection + 1;
                //&cout << "rejected: " << rejection << endl;
            }

			// total energy printout	
            //&cout << "Etotal from particle move: " << ETotal+Ecorr << endl;

            // check
            Esum = 0;
            //double ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength);
            //&cout << "Etotal check: " << ETotalcheck+Ecorr << endl;
	
			//&cout << "rmax = " << rmax << endl;
            //&cout << "--------------------------" << endl;

            // optimze rmax every 1000 loops
            if ((trial)%1000 == 0)
            {
			//	cout << "density: " << boxObj.rho(BoxLength) << endl;
			//	outFile << trial << "   " << (boxObj.rho(BoxLength)) << endl;
            	if(trial>0) rmax = 2*(double(acceptance)/trial)*rmax;
            }
			
			//outFile << trial << "   " << ETotal+Ecorr <<  endl;
			//outFile << trial << "   " << boxObj.rho(BoxLength) <<  endl;
        	}


	}

	end = clock();

	cout << "CPU time: " << ((float)(end - start)/CLOCKS_PER_SEC) << endl;

	outFile.close();
	return 0;
} 



