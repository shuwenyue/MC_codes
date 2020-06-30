#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h> 
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
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
	//double epsilon = 0.1654; // kg*angstroms/s^2
	//double sigma = 3.4; // angstroms
	//double T = (0.9)*epsilon/Kb; // K
    //double P = (7.36626e-3)*epsilon/(sigma*sigma*sigma); // kg/(angstrom*s^2)

	int TrialNumber = 1000000;

	double sigma = 1;
	double epsilon = 1;

	double T = 1.0;
	double beta = 1/T;
	double rho1 = 0.3;
	double rho2 = 0.3;
	int NParticles1 = 216;
	int NParticles2 = 216;
	int NTotalParticles = NParticles1 + NParticles2;
	double BoxLength1 = pow((NParticles1/rho1), 1.0/3.0);
	double BoxLength2 = pow((NParticles2/rho2), 1.0/3.0);

	double cutoff = 3; 

	// particle displacement variables
	double ETotal1, ETotal2, E1, E2, deltaE, Ecorr;
	double rmax1 = 0.5; 	
	double rmax2 = 0.5; 	
	int particlebox1move = 1;
	int particlebox2move = 1;
	int partMoveBox1acc = 0;
    int partMoveBox2acc = 0;
    int partMoveBox1rej = 0;
    int partMoveBox2rej = 0;
	int randparticle;

	// volume move variables
	double deltaETotal1, deltaETotal2, ETotal1new, ETotal2new, EcorrV1, EcorrV2;
	double deltavmax = 1*BoxLength1*BoxLength1*BoxLength1;
	double oldV1, oldV2, newV1, newV2, totalV, newBoxLength1, newBoxLength2;
    int volumemove = 1;
    int Vacceptance = 0;
    int Vrejection = 0;
    
    // swap variables
    double EDeletedParticle, EAddedParticle;
	int SwapAcceptance = 0;
	int SwapRejection = 0;
	int swapMove = 1;

	clock_t p1, p2, p3, p4, p5, p6, v1, v2, s1, s2;

	// Initialize positions
	particle particleObj;
	particleObj.setPosition(BoxLength1, BoxLength2, NParticles1, NParticles2);

	// calculate inital total energy
	double Esum = 0;
	energy energyObj(Esum, sigma, epsilon, cutoff);
	ETotal1 = energyObj.getTotalEnergy(1, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength1, NParticles1);
	ETotal2 = energyObj.getTotalEnergy(2, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength2, NParticles2);

	//cout << "initial ETotal 1: " << ETotal1<< endl;
	//cout << "initial ETotal 2: " << ETotal2<< endl;

	// seed random number generator
	srand(time(0));

	// print to output file
   	ofstream outFile;
    outFile.open("test.txt");
	//FILE * outFile;
	//outFile = fopen ("Energies.txt","w+");

	/*if (outFile.is_open())
	{
		cout << "file is open" << endl;
	}*/

	// begin trial loop
   	for (int trial = 1; trial < TrialNumber; trial++)
    {
		 //cout << "Trial: " << trial << endl;
		// outFile << trial  << "  " << NParticles1/(BoxLength1*BoxLength1*BoxLength1) << "   " << NParticles2/(BoxLength2*BoxLength2*BoxLength2)  << endl; //<< "  " << ((float)t/CLOCKS_PER_SEC) << endl;
        // chose particle displacement or volume change
        double whichmove;
       /* double test1 = ETotal1 -(energyObj.getTotalEnergy(1, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength1, NParticles1));
        double test2 = ETotal2 -(energyObj.getTotalEnergy(2, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength2, NParticles2));
        if ((test1*test1) > 1e-12 || (test2*test2) > 1e-12)
        {
            cout << "test1" << test1 << ", test2" << test2 << endl;
            cout << "number of trials: " << trial << endl;
            cout << "whichmove: " << whichmove << endl;
            
            exit(0);
        }*/
        
        
        
        
        double randmove = ((double) rand() / (RAND_MAX));

        int particlemoves = 100;
        int volumemoves = 1;
        int swaps = 200;

       	whichmove = randmove * (particlemoves+volumemoves+swaps);
	
                //cout << "Etotal1 check: " << energyObj.getTotalEnergy(1, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength1, NParticles1) << endl;		//whichmove = particlemoves+volumemoves+swaps;  //CHECK
        
		if (whichmove < particlemoves)
        {
            p1 = clock();
            
            // perform particle displacement
	
           	//chose whichbox to move particle
			int randparticle = rand()%((NTotalParticles));
			
			//whichbox = 0; //CHECK
           	if (particleObj.getPosition(randparticle).n == 1 && NParticles1 > 0)
			{

          		//cout << "particle displacement in box 1: " << particlebox1move << endl;
          		//cout << "Initial Total energy box 1: " << ETotal1<< endl;
				
                p2 = clock();
                
				// save old coordinates
                double4 OldCoord = particleObj.getPosition(randparticle);
                
				//cout << "rand particle in box1 : " << randparticle << " - box: " << particleObj.getPosition(randparticle).n << endl;
		
                // calculate energy of random particle
                E1 = 0;
                E1 = energyObj.getIndividualEnergy(1, particleObj.getVector(), BoxLength1, NTotalParticles, randparticle);
                //cout << "E1= " << E1 << endl;
                
                p3 = clock();
                
                // random number between 0 and 1
                double a = ((double) rand() / (RAND_MAX)) ;
                double b = ((double) rand() / (RAND_MAX)) ;
                double c = ((double) rand() / (RAND_MAX)) ;
                

				//cout << "old position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;

                // adjust position to new coordinates
                particleObj.moveParticle(a,b,c,rmax1,randparticle);
                
                			//cout << "new position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;

                p4 = clock();
                
                // calculate new energy
                E2 = 0;
                E2 = energyObj.getIndividualEnergy(1, particleObj.getVector(), BoxLength1, NTotalParticles, randparticle);
                //cout << "E2= " << E2 << endl;

                p5 = clock();
                
                // calculate energy difference
                double deltaE = E2 - E1;
                //cout << "deltaE: " << deltaE << endl;

                // calculate acceptance criteria
                double acc = exp(-deltaE/T);
                //cout << "acceptance criteria: " << acc << endl;

                // accept new values if criteria is met
                if (((double) rand() / (RAND_MAX)) < acc)
                {
                    ETotal1 = ETotal1 + deltaE;
                    partMoveBox1acc = partMoveBox1acc + 1;
                    //cout << "particle move accepted: " << partMoveBox1acc << " out of " << particlebox1move << " moves in box 1" << endl;
                    //cout << "accepted position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;
                }
                else
                {
                    particleObj.restoreParticleMove(OldCoord, randparticle);
                    partMoveBox1rej = partMoveBox1rej + 1;
                    //cout << "particle move rejected: " << partMoveBox1rej << " out of " << particlebox1move << " moves in box 1" << endl;
                    //cout << "rejected position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;
                }

                p6 = clock();
                
                // total energy printout
               	//cout << "Etotal1 from particle move: " << ETotal1<< endl;

				// rmax printout
                //cout << "rmax1 = " << rmax1 << endl;

				// Total Energy check
                //cout << "Etotal1 check: " << energyObj.getTotalEnergy(1, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength1, NParticles1) << endl;

                // optimze rmax1 every 100 loops
                if ((particlebox1move)%100 == 0)
                {
                    if(particlebox1move>0) rmax1 = 2*((double)partMoveBox1acc/(double)particlebox1move)*rmax1;
                }
				
                particlebox1move = particlebox1move + 1;

            }
        
            else if (particleObj.getPosition(randparticle).n == 2 && NParticles2 > 0)
            {
          		//cout << "particle displacement in box 2: " << particlebox2move << endl;
          		//cout << "Initial Total energy box 2: " << ETotal2<< endl;
         
                //cout << "randparticle: " << randparticle << endl;

				// save old coordinates
                double4 OldCoord = particleObj.getPosition(randparticle);

                // calculate energy of random particle
                E1 = 0;
                E1 = energyObj.getIndividualEnergy(2, particleObj.getVector(), BoxLength2, NTotalParticles, randparticle);
                //cout << "E1= " << E1 << endl;

                // random number between 0 and 1
                double a = ((double) rand() / (RAND_MAX)) ;
                double b = ((double) rand() / (RAND_MAX)) ;
                double c = ((double) rand() / (RAND_MAX)) ;

				//cout << "old position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;

                // adjust position to new coordinates
                particleObj.moveParticle(a,b,c,rmax2,randparticle);
			
				//cout << "new position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;


                // calculate new energy
                E2 = 0;
                E2 = energyObj.getIndividualEnergy(2, particleObj.getVector(), BoxLength2, NTotalParticles, randparticle);
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
                    ETotal2 = ETotal2 + deltaE;
                    partMoveBox2acc = partMoveBox2acc + 1;
                    //cout << "particle move accepted: " << partMoveBox2acc << " out of " << particlebox2move << " moves in box 2" << endl;
                    //cout << "accepted position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;
                }
                else
                {
                    particleObj.restoreParticleMove(OldCoord, randparticle);
                    partMoveBox2rej = partMoveBox2rej + 1;
                    //cout << "particle move rejected: " << partMoveBox2rej << " out of " << particlebox2move << " moves in box 2" << endl;
                    //cout << "rejected position check: " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << " in box " << particleObj.getPosition(randparticle).n << endl;
                }

                // total energy printout
               	//cout << "Etotal2 from particle move: " << ETotal2<< endl;

				// rmax printout
                //cout << "rmax2 = " << rmax2 << endl;

				// Total Energy check
                //cout << "Etotal2 check: " << energyObj.getTotalEnergy(2, particleObj.getVector(), BoxLength1, BoxLength2, NParticles1, NParticles2)+energyObj.getECorrection(BoxLength2, NParticles2) << endl;

                // optimze rmax2 every 100 loops
                if ((particlebox2move)%100 == 0)
                {
                    if(particlebox2move>0) rmax2 = 2*((double)partMoveBox2acc/(double)particlebox2move)*rmax2;
                }
				
                particlebox2move = particlebox2move + 1;
			}
            p2 = clock();
        }
       	else if (whichmove < (particlemoves+volumemoves))
        {
            v1 = clock();
        
            //cout << "volume move: " << volumemove << endl;

			double Ecorr1Old = energyObj.getECorrection(BoxLength1, NParticles1);
			double Ecorr2Old = energyObj.getECorrection(BoxLength2, NParticles2);
            //cout << "Initial ETotal1 = " << ETotal1<< endl;
            //cout << "Initial ETotal2 = " << ETotal2<< endl;
        
            // save old volume and coordinates
            oldV1 = BoxLength1*BoxLength1*BoxLength1;
            oldV2 = BoxLength2*BoxLength2*BoxLength2;
            totalV = oldV1 + oldV2;
            vector<double4> OldVector = particleObj.getVector();

            // random number to change V1
            double d = ((double) rand() / (RAND_MAX)) ;
        
            // adjust volume of V1
            //cout << "oldV1=" << oldV1 << endl;
            //double lognewV1 = log(oldV1) + (2*d-1) * deltavmax; 
			//newV1 = exp((double)lognewV1);
			newV1 = oldV1 + (d-0.5)*deltavmax;
            //cout << "newV1=" << newV1 << endl;
			//cout << "old BoxLength1: " << BoxLength1 << endl;
            newBoxLength1 = pow((double)newV1, 1.0/3.0);
            //cout << "newBoxLength1: " << newBoxLength1  << endl;
        
            // adjust volume of V2
            //cout << "oldV2=" << oldV2 << endl;
            newV2 = totalV - newV1;
            //cout << "newV2=" << newV2 << endl;
			//cout << "old BoxLength2: " << BoxLength2 << endl;
            newBoxLength2 = pow((double)newV2, 1.0/3.0);
            //cout << "newBoxLength2: " << newBoxLength2  << endl;

            // rescale particles
            particleObj.rescaleBox(1, (newBoxLength1/BoxLength1), NTotalParticles);
    		particleObj.rescaleBox(2, (newBoxLength2/BoxLength2), NTotalParticles);

            // calculate new total energy
			ETotal1new = energyObj.getTotalEnergy(1, particleObj.getVector(), newBoxLength1, newBoxLength2, NParticles1, NParticles2);
			ETotal2new = energyObj.getTotalEnergy(2, particleObj.getVector(), newBoxLength1, newBoxLength2, NParticles1, NParticles2);
			double Ecorr1New = energyObj.getECorrection(newBoxLength1, NParticles1);
			double Ecorr2New = energyObj.getECorrection(newBoxLength2, NParticles2);
	
            //cout << "ETotal1 new= " << ETotal1new + Ecorr1New  << endl;
            //cout << "ETotal2 new= " << ETotal2new + Ecorr2New  << endl;
            
            // acceptance criteria
            deltaETotal1 = (ETotal1new + Ecorr1New) - (ETotal1);
            deltaETotal2 = (ETotal2new + Ecorr2New) - (ETotal2);
            //cout << "deltaETotal1 = " << deltaETotal1 << endl;
            //cout << "deltaETotal2 = " << deltaETotal2 << endl;
            double accV = exp(-beta*(deltaETotal1)+(NParticles1)*log(newV1/oldV1)-beta*(deltaETotal2)+(NParticles2)*log(newV2/oldV2));
            //cout << "acceptance criteria: " << accV << endl;
            
            // accept new values if criteria is met
            if (((double) rand() / (RAND_MAX)) < accV)
            {
				if ((newBoxLength1/2) < cutoff || (newBoxLength2/2) < cutoff)
				{
					cout << "box length too small!" << endl;
					exit(0);
				}

                Vacceptance = Vacceptance + 1;
                //cout << "Vchange accepted: " << Vacceptance << endl;
                BoxLength1 = newBoxLength1;
                BoxLength2 = newBoxLength2;
                ETotal1 = ETotal1new+Ecorr1New;
                ETotal2 = ETotal2new+Ecorr2New;
                //cout << "E box1 from volume move: " << ETotal1<< endl;
                //cout << "E box2 from volume move: " << ETotal2<< endl;
            }
            else
            {
                particleObj.restoreBox(OldVector, NTotalParticles);
                Vrejection = Vrejection + 1;
                //cout << "Vchange rejected: " << Vrejection << endl;
				//cout << "E box1 from volume move (rej): " << ETotal1 << endl;
                //cout << "E box2 from volume move (rej): " << ETotal2 << endl;
            }

            //cout << "deltavmax= " << deltavmax << endl;

            // optimze deltavmax1 every 100 loops
            if ((volumemove)%100 == 0)
            {
                if(volumemove>0) deltavmax = 2*(double(Vacceptance)/volumemove)*deltavmax;
            }
            
            volumemove = volumemove + 1;
            v2 = clock();
        }
        else //if (whichmove > (particlemoves+volumemoves))
     	{
            s1 = clock();
            
            
            // perform swap
            //cout << "swap move: " << swapMove << endl;

			//chose whichbox to move particle
			double whichbox = ((double) rand() / (RAND_MAX));
            //whichbox = 1; // CHECK
            
           	if (whichbox < 0.5 && NParticles1 > 0)
			{
				//cout << "particle move from box 1 to box 2" << endl;
				//cout << "initial NParticles1 : " << NParticles1 << endl;
				//cout << "initial NParticles2 : " << NParticles2 << endl;

				// check initial energy
				double Ecorr1Old = energyObj.getECorrection(BoxLength1, NParticles1);
				double Ecorr2Old = energyObj.getECorrection(BoxLength2, NParticles2);
				//cout << "ETotal1 before swap: " << ETotal1 << endl;
				//cout << "ETotal2 before swap: " << ETotal2 << endl;


				// choose random particle
				int id = int(((double) rand() / (RAND_MAX)) * NParticles1);
				int k=0;
				int i;

				for (i=0; i < NTotalParticles; i++)
				{
                    
					if (particleObj.getPosition(i).n == 1)
					{
                        if (k == id) break;
						k++;
					}
										
				}
				double randparticle = i;

                //cout << "before swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

            	// save particle coordinates
            	double4 OldCoord = particleObj.getPosition(randparticle);

				// individual energy calculation of deleted particle in box 1 prior to deletion 
              	EDeletedParticle = energyObj.getIndividualEnergy(1, particleObj.getVector(), BoxLength1, NTotalParticles, randparticle);
				
				//cout << "EDeletedParticle: " << EDeletedParticle << endl;

				// random position for new particle in box2
                double randx = BoxLength2 * ((double) rand() / (RAND_MAX));
                double randy = BoxLength2 * ((double) rand() / (RAND_MAX));
                double randz = BoxLength2 * ((double) rand() / (RAND_MAX));

				// move particle from box 1 to new position in box 2
				particleObj.swapParticle(randparticle, randx, randy, randz, 2);
				NParticles1 = NParticles1 - 1;
				NParticles2 = NParticles2 + 1;

				double EcorrDeleted = energyObj.getECorrection(BoxLength1, NParticles1);
                //cout << "after swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

				// individual energy calculation of newly added particle in box 2 
              	EAddedParticle = energyObj.getIndividualEnergy(2, particleObj.getVector(), BoxLength2, NTotalParticles, randparticle);
				double EcorrAdded = energyObj.getECorrection(BoxLength2, NParticles2);

				double deltaE1Swap = -EDeletedParticle+(EcorrDeleted - Ecorr1Old);
				double deltaE2Swap = EAddedParticle+(EcorrAdded - Ecorr2Old);
	
				//cout << "EAddedParticle: " << EAddedParticle << endl;
				//cout << "deltaE1Swap: " << deltaE1Swap << endl;
				//cout << "deltaE2Swap: " << deltaE2Swap << endl;

				// acceptance criteria
				double accSwap = ((NParticles1+1)*BoxLength2*BoxLength2*BoxLength2/((NParticles2)*BoxLength1*BoxLength1*BoxLength1))*exp((-beta*(deltaE1Swap))-(beta*(deltaE2Swap)));
				//double accSwap = ((NParticles1*BoxLength2*BoxLength2*BoxLength2/((NParticles2+1)*BoxLength1*BoxLength1*BoxLength1))*exp((-beta*(deltaE1Swap))-(beta*(deltaE2Swap)));
				//double accSwap = exp(log(BoxLength2*BoxLength2*BoxLength2/BoxLength1*BoxLength1*BoxLength1)-beta*deltaE1Swap-beta*deltaE2Swap);

				//cout << "acc: " << accSwap << endl;

				if (((double) rand() / (RAND_MAX)) < accSwap)
            	{
                	SwapAcceptance = SwapAcceptance + 1;
                	//cout << "Swap accepted: " << SwapAcceptance << endl;
                	ETotal1 = ETotal1 + deltaE1Swap;  
                	ETotal2 = ETotal2 + deltaE2Swap;

                    //cout << "accepted swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
               
           	 	}
				else
				{
					SwapRejection = SwapRejection + 1;
					//cout << "Swap rejected: " << SwapRejection << endl;

					// restore particle back in box 1
					particleObj.swapParticle(randparticle, OldCoord.x, OldCoord.y, OldCoord.z, 1);
					NParticles1 = NParticles1 + 1;
					NParticles2 = NParticles2 - 1;
                    //cout << "rejected swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
                    
				}

                swapMove = swapMove + 1;

				//cout << "swaps accepted: " << SwapAcceptance << " out of " << swapMove << endl;

			}
			else if (whichbox>0.5 && NParticles2 > 0)
			{
                
                
      			//cout << "particle move from box 2 to box 1" << endl;
				//cout << "initial NParticles1 : " << NParticles1 << endl;
				//cout << "initial NParticles2 : " << NParticles2 << endl;

				// check initial energy
				double Ecorr1Old = energyObj.getECorrection(BoxLength1, NParticles1);
				double Ecorr2Old = energyObj.getECorrection(BoxLength2, NParticles2);
				//cout << "ETotal1 before swap: " << ETotal1 << endl;
				//cout << "ETotal2 before swap: " << ETotal2 << endl;
		

                
                
				// choose random particle
				int id = int(((double) rand() / (RAND_MAX)) * NParticles2);
                
                
				int k=0;
				int i;

                if (trial == 659)
                {

                    cout << " id: " << id << endl;
                    //exit(0);
                }
				for (i=0; i < NTotalParticles; i++)
				{
                    
					if (particleObj.getPosition(i).n == 2)
					{
                        if (k == id) break;
						k++;
					}
											
				}
				double randparticle = i;


                //cout << "before swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

            	// save particle coordinates
            	double4 OldCoord = particleObj.getPosition(randparticle);


                
				// individual energy calculation of deleted particle in box 2 prior to deletion
              	EDeletedParticle = energyObj.getIndividualEnergy(2, particleObj.getVector(), BoxLength2, NTotalParticles, randparticle);
              
                /*if (trial == 659)
                {
                    cout << "Npart1 : " << NParticles1 << ", Npart2 : " << NParticles2 << endl;
                    //cout << "Etotal1: " << ETotal1 << endl;
                    //cout << "Etotal2: " << ETotal2 << endl;
                    //cout << "whichbox: " << whichbox << endl;
                    cout << "randparticle" << randparticle << endl;
                    //cout << "Etotal2: " << ETotal2 << endl;
                    //cout << "deltaEswap: " << deltaE2Swap << endl;
                    cout << "Edeleteparticle: " << EDeletedParticle << endl;
                    //cout << "EcorrDeleted: " << EcorrDeleted << endl;
                    //cout << "ecorr2Old: " << Ecorr2Old << endl;
                    //exit(0);
                 }*/
            
				//cout << "EDeletedParticle: " << EDeletedParticle << endl;

				// random position for new particle in box1
                double randx = BoxLength1 * ((double) rand() / (RAND_MAX));
                double randy = BoxLength1 * ((double) rand() / (RAND_MAX));
                double randz = BoxLength1 * ((double) rand() / (RAND_MAX));

				// move particle from box 2 to new position in box 1
				particleObj.swapParticle(randparticle, randx, randy, randz, 1);
				NParticles1 = NParticles1 + 1;
				NParticles2 = NParticles2 - 1;

				double EcorrDeleted = energyObj.getECorrection(BoxLength2, NParticles2);
                //cout << "after swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;

				// individual energy calculation of newly added particle in box 1 
              	EAddedParticle = energyObj.getIndividualEnergy(1, particleObj.getVector(), BoxLength1, NTotalParticles, randparticle);
				double EcorrAdded = energyObj.getECorrection(BoxLength1, NParticles1);

				double deltaE1Swap = EAddedParticle+(EcorrAdded - Ecorr1Old);
				double deltaE2Swap = -EDeletedParticle+(EcorrDeleted - Ecorr2Old);

                /*if (trial == 746)
                {
                    cout << "Npart1 : " << NParticles1 << ", Npart2 : " << NParticles2 << endl;
                    cout << "Etotal1: " << ETotal1 << endl;
                    cout << "Etotal2: " << ETotal2 << endl;
                    cout << "whichbox: " << whichbox << endl;
            
                    cout << "Etotal2: " << ETotal2 << endl;
                    cout << "deltaEswap: " << deltaE2Swap << endl;
                    cout << "Edeleteparticle: " << EDeletedParticle << endl;
                    cout << "EcorrDeleted: " << EcorrDeleted << endl;
                    cout << "ecorr2Old: " << Ecorr2Old << endl;
                    var = 1;
                }*/
                
				//cout << "EAddedParticle: " << EAddedParticle << endl;
				//cout << "deltaE1Swap: " << deltaE1Swap << endl;
				//cout << "deltaE2Swap: " << deltaE2Swap << endl;

				// acceptance criteria
				double accSwap = ((NParticles2+1)*BoxLength1*BoxLength1*BoxLength1/((NParticles1)*BoxLength2*BoxLength2*BoxLength2))*exp((-beta*(deltaE1Swap))-(beta*(deltaE2Swap)));
                //double accSwap = (NParticles2*BoxLength1*BoxLength1*BoxLength1/((NParticles1+1)*BoxLength2*BoxLength2*BoxLength2))*exp((-beta*(deltaE1Swap))-(beta*(deltaE2Swap)));
				//double accSwap = exp(log(BoxLength1*BoxLength1*BoxLength1/BoxLength2*BoxLength2*BoxLength2)-beta*deltaE2Swap-beta*deltaE1Swap);
				//cout << "acc: " << accSwap << endl;
                
                

                
                
				if (((double) rand() / (RAND_MAX)) < accSwap)
            	{
                	SwapAcceptance = SwapAcceptance + 1;
                	//cout << "Swap accepted: " << SwapAcceptance << endl;
                	ETotal1 = ETotal1 + deltaE1Swap;  
                	ETotal2 = ETotal2 + deltaE2Swap;

//cout << "accepted swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
           	 	}
            	else
            	{
					SwapRejection = SwapRejection + 1;
					//cout << "Swap rejected: " << SwapRejection << endl;

					// restore particle back in box 2
					particleObj.swapParticle(randparticle, OldCoord.x, OldCoord.y, OldCoord.z, 2);
					NParticles1 = NParticles1 - 1;
					NParticles2 = NParticles2 + 1;

//cout << "rejected swap - randparticle: " << randparticle << " in box " << particleObj.getPosition(randparticle).n << " at " << particleObj.getPosition(randparticle).x << "," << particleObj.getPosition(randparticle).y << "," << particleObj.getPosition(randparticle).z << endl;
           		} 

				//cout << "swaps accepted: " << SwapAcceptance << " out of " << swapMove << endl;
                swapMove = swapMove + 1;
            	
            }
            s2 = clock();
		}
        
        //cout << "Final NParticles1: " << NParticles1 << endl;
        //cout << "Final NParticles2: " << NParticles2 << endl;
        //cout << "Final Box 1 Energy: " << ETotal1 << endl;
        //cout << "Final Box 2 Energy: " << ETotal2 << endl;
   
		//t = clock() - t; //end timer
		if ((trial)%10000 == 0)
        {
            //cout << "NParticles total: " << NParticles1 + NParticles2 << endl;
            //cout << trial << "  " << NParticles1/(BoxLength1*BoxLength1*BoxLength1) << "   " << NParticles2/(BoxLength2*BoxLength2*BoxLength2) << endl; //"  " << p2-p1<< "," << p3-p2 << "," << p4-p3 << "," << p5-p4 << "," << p6-p5 << endl;//" / " << v2-v1 << " / " << s2-s1 << endl;
            //cout << "NParticles1: " << NParticles1 << " / " << "NParticles2: " << NParticles2 << endl;

            outFile << trial  << "  " << NParticles1/(BoxLength1*BoxLength1*BoxLength1) << "   " << NParticles2/(BoxLength2*BoxLength2*BoxLength2) << endl; //<< ((float)t/CLOCKS_PER_SEC) << endl;
		cout << "energy: " << ETotal1 << "," << ETotal2 << endl;

		//	fprintf (outFile, "%d  %f  %f  \n", trial, NParticles1/(BoxLength1*BoxLength1*BoxLength1), NParticles2/(BoxLength2*BoxLength2*BoxLength2)); //, (float)t/CLOCKS_PER_SEC);
			//fflush (outFile);
			//outFile.flush();
		}
		//cout << "--------------------------" << endl;

    }
	
	//fclose (outFile);
	outFile.close(); 

	return 0;
}



