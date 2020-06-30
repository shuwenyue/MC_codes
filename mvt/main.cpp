#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include <cmath>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"
#include "histogram.h"

using namespace std;

int main () 
{
	int TrialNumber = 1000000;

	double sigma = 1;
	double epsilon = 1;

	double T = 1.3;
	double beta = 1/T;
	double mu = 1;
    double lambda = 1;
	double rho = 0.3;
	int NParticles = 216;
	//double BoxLength = 7.56;
	//int NParticles = rho*BoxLength*BoxLength*BoxLength; 
	double BoxLength = pow(((double)NParticles/(double)rho), 1.0/3.0);
    //cout << "NParticles: " << NParticles << endl;
	double cutoff = 3; 
	cout << "BoxLength: " << BoxLength << endl;

    double whichmove;

	// particle displacement variables
	double ETotal, E1, E2, deltaE, Ecorr;
	double rmax = 0.5; 		
	int particlemove = 1;
	int partmoveAcc = 0;
    int partmoveRej = 0;
 
    // add/delete variables
    double EDeletedParticle, EAddedParticle;
	int exchangeAcc = 0;
	int exchangeRej = 0;
	int exchangemove = 1;

	// Initialize positions
	particle particleObj;
	particleObj.setPosition(BoxLength, NParticles);
    
	// calculate inital total energy
	double Esum = 0;
	energy energyObj(Esum, sigma, epsilon, cutoff);
	ETotal = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength, NParticles)+energyObj.getECorrection(BoxLength, NParticles);
	//cout << "initial ETotal: " << ETotal << endl;

	// seed random number generator
	srand(time(0));

	// print to output file
   	ofstream outFile;
    outFile.open("Energies.txt");

	// begin trial loop
   	for (int trial = 1; trial < TrialNumber; trial++)
    {
   
    
		//cout << "Trial: " << trial << endl;

        // chose particle displacement or exchange move
        double randmove = ((double) rand() / (RAND_MAX));

        int numberOfParticleMoves = 60;
        int numberOfExchangeMoves = 40;

       	whichmove = randmove * (numberOfParticleMoves+numberOfExchangeMoves);
        
        // ETotal check
        double ETotalCheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength, NParticles)+energyObj.getECorrection(BoxLength, NParticles);
        
        if ((ETotalCheck-ETotal)*(ETotalCheck-ETotal) > 0.00000001)
        {
            cout << "different ETotal:  " << ETotalCheck - ETotal << endl;
            cout << "trial: " << trial << endl;
            cout << "whichmove: " << whichmove << endl;
            
            exit(0);
        }
        
		//whichmove = 100; //CHECK

		if (whichmove < numberOfParticleMoves && NParticles > 0)
        {

          	//cout << "particle move: " << particlemove << endl;
          	//cout << "Initial Total energy: " << ETotal << endl;
     
            // choose a random particle from reservoir to move
			int randparticle = rand()%(NParticles);

			//cout << "randparticle: " << randparticle << endl;

            // save old coordinates
            double3 OldCoord = particleObj.getPosition(randparticle);

            // calculate energy of random particle
            E1 = 0;
            E1 = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, NParticles, randparticle);
            //cout << "E1= " << E1 << endl;

            // random number between 0 and 1
            double a = ((double) rand() / (RAND_MAX));
            double b = ((double) rand() / (RAND_MAX));
            double c = ((double) rand() / (RAND_MAX));

            //cout << "old position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

            // adjust position to new coordinates
            particleObj.moveParticle(a,b,c,rmax,randparticle);
			
            //cout << "new position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

            // calculate new energy
            E2 = 0;
            E2 = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, NParticles, randparticle);
            //cout << "E2= " << E2 << endl;

            // calculate energy difference
            double deltaE = E2 - E1;
            //cout << "deltaE: " << deltaE << endl;

            // calculate acceptance criteria
            double acc = exp(-deltaE/T);
            //cout << "acceptance criteria: " << acc << endl;
 
            // accept new values if criteria is met
            if (((double) rand() / (RAND_MAX)) < acc)
            {
            	ETotal = ETotal + deltaE;
                partmoveAcc = partmoveAcc + 1;
                //cout << "particle move in box 1 accepted: " << partmoveAcc << " out of " << particlemove << " moves" << endl;

				//cout << "position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
            }
            else
            {
            	particleObj.restoreParticleMove(OldCoord, randparticle);
                partmoveRej = partmoveRej + 1;
                //cout << "particle move in box 1 rejected: " << partmoveRej << " out of " << particlemove << " moves" << endl;

				//cout << "position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
            }

            // total energy printout
            //cout << "Etotal from particle move: " << ETotal << endl;
           
            // check
			double ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength, NParticles)+energyObj.getECorrection(BoxLength, NParticles);
			//cout << "ETotalCheck: " << ETotalcheck << endl;

            // optimze rmax every 100 loops
            if ((particlemove)%100 == 0)
            {
            	if(particlemove>0) rmax = 2*((double)partmoveAcc/(double)particlemove)*rmax;
            }
				
            particlemove = particlemove + 1;
			//cout << "rmax: " << rmax << endl;
        }
        else 
        {

  			// Exchange move - add or delete?
         	double whichmove = ((double) rand() / (RAND_MAX));
            
			if (whichmove < 0.5)
			{

            	//cout << "exchange move: " << exchangemove << " - add" << endl;
				//cout << "number of particles (before): " << NParticles << endl;
				double EcorrOld = energyObj.getECorrection(BoxLength, NParticles);
				//cout << "Total Energy (before): " << ETotal << endl;
				double ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength, NParticles)+energyObj.getECorrection(BoxLength, NParticles);
				//cout << "ETotalCheck: " << ETotalcheck << endl;

				// random position for new particle
                double randx = BoxLength * ((double) rand() / (RAND_MAX));
                double randy = BoxLength * ((double) rand() / (RAND_MAX));
                double randz = BoxLength * ((double) rand() / (RAND_MAX));
                
  				// add random particle to box
                particleObj.addParticle(randx, randy, randz);
				NParticles = NParticles + 1;
               
                
				// individual energy calculation of newly added particle
              	EAddedParticle = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, NParticles, NParticles-1);
                double EcorrNew = energyObj.getECorrection(BoxLength, NParticles);
				double deltaE = EAddedParticle + (EcorrNew-EcorrOld);
                //cout << "E added Particle: " << EAddedParticle << endl;
				//cout << "deltaE: " << deltaE << endl;
			
                
                //double lambda = sqrt((6.626e-34)*(6.626e-34)/(2*3.14159265359*6.6335209e-26*T));
                
                //double muPrime = mu - T*log(lambda*lambda*lambda*T);
                
				// acceptance criteria
				double accAdd = (BoxLength*BoxLength*BoxLength/(lambda*lambda*lambda*NParticles))*exp(beta*mu-beta*deltaE);
                //double accAdd = -beta*deltaE+beta*muPrime+log((BoxLength*BoxLength*BoxLength)/(NParticles));
                
                //cout << "acceptance criteria: " << accAdd << endl;
                
                //accAdd = 1; // CHECK
				if (((double) rand() / (RAND_MAX)) < accAdd)
            	{
                	exchangeAcc = exchangeAcc + 1;
                	//cout << "Add accepted: " << exchangeAcc << endl;
                	ETotal = ETotal + deltaE;
                    //cout << "vector size: " << particleObj.getVector().size() << endl;
           	 	}
            	else
            	{
					exchangeRej = exchangeRej + 1;
                    particleObj.getVector().pop_back();
					//cout << "Add rejected: " << exchangeRej << endl;
					NParticles = NParticles - 1;
                    //cout << "vector size: " << particleObj.getVector().size() << endl;
           		}
                
                //cout << "Etotal: " << ETotal << endl;

                
				// calculate density
				double rho = exp(mu/T)/(lambda*lambda*lambda);

				// ETotal check
				ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength, NParticles)+energyObj.getECorrection(BoxLength, NParticles);
				//cout << "ETotalCheck: " << ETotalcheck << endl;

                

			}
			else if (NParticles > 0)
			{
                
                
				//cout << "exchange move: " << exchangemove << " - delete" << endl;
				//cout << "number of particles (before): " << NParticles << endl;
				double EcorrOld = energyObj.getECorrection(BoxLength, NParticles);
				//cout << "Total Energy (before): " << ETotal << endl;
				double ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength,NParticles)+energyObj.getECorrection(BoxLength, NParticles);
				//cout << "ETotalCheck: " << ETotalcheck << endl;

 				// choose a random particle from box to remove
				int randparticle = rand()%(NParticles);
                //cout << "randparticle: " << randparticle << endl;

				// save coordinates of old position
				double3 OldCoord = particleObj.getPosition(randparticle);
                
				// individual energy calculation of deleted particle prior to deletion
              	EDeletedParticle = energyObj.getIndividualEnergy(particleObj.getVector(), BoxLength, NParticles, randparticle);
                //cout << "EDeleted Particle: " << EDeletedParticle << endl;
                
                // CHECKies
                /*for (int t=0; t<NParticles; t++)
                {
                    cout << t << ": " << particleObj.getPosition(t).x << ", " << particleObj.getPosition(t).y << ", " << particleObj.getPosition(t).z << endl;
                } // CHECK */
                
                
				// delete random particle
                particleObj.deleteParticle(randparticle, NParticles);
				NParticles = NParticles - 1;
                

				double EcorrNew = energyObj.getECorrection(BoxLength, NParticles);

				double deltaE = -EDeletedParticle + (EcorrNew-EcorrOld);
				//cout << "deltaE: " << deltaE << endl;
                
				//double lambda = sqrt((6.626e-34)*(6.626e-34)/(2*3.14159265359*6.6335209e-26*T));

                //double muPrime = mu - T*log(lambda*lambda*lambda*T);
                
				// acceptance criteria
				double accDel = ((lambda*lambda*lambda*(NParticles+1))/(BoxLength*BoxLength*BoxLength))*exp(-beta*mu-beta*deltaE);
                //double accDel = -beta*muPrime - beta*deltaE + log((NParticles+1)/(BoxLength*BoxLength*BoxLength));
 
				//cout << "acceptance criteria: " << accDel << endl;
                
				if (((double) rand() / (RAND_MAX)) < accDel)
            	{
                	exchangeAcc = exchangeAcc + 1;
                	//cout << "Delete accepted: " << exchangeAcc << endl;
                	ETotal = ETotal + deltaE;
                    //cout << "vector size: " << particleObj.getVector().size() << endl;
           	 	}
            	else
            	{
					exchangeRej = exchangeRej + 1;
					//cout << "Delete rejected: " << exchangeRej << endl;
                    NParticles = NParticles + 1;
                    particleObj.addParticle(OldCoord.x, OldCoord.y, OldCoord.z);
           		}
                

				//ETotalcheck = energyObj.getTotalEnergy(particleObj.getVector(), BoxLength,NParticles)+energyObj.getECorrection(BoxLength, NParticles);
				/*if (abs(ETotalcheck -ETotal) < 0.000001)
				{
					cout << "ETotals don't match" << endl;
					exit(0);
				}*/
				
   
            }

            exchangemove = exchangemove + 1;

        }

		//cout << "NParticles in box: " << NParticles << endl;
		//cout << "Final Box Energy: " << ETotal << endl;


		//cout << "--------------------------" << endl;

		if ((trial)%100 == 0)
        {
        	outFile << trial << "   " << NParticles/(BoxLength*BoxLength*BoxLength) << endl;
			cout << "energy: " << ETotal << endl;
            //cout << trial << "----" << NParticles << "--- Energy: " << ETotal << endl;

        }
    }
  
	outFile.close();
    
    histogram histogramObj;
    histogramObj.plot(T, mu, BoxLength);
    
	return 0;
}



