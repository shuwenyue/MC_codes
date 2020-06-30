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

void particle::setPosition(double BoxLength1, double BoxLength2, int NParticles1, int NParticles2)
{
	int particles = 0; 
	
	positionVector.resize(NParticles1 + NParticles2);

	double ParticleLength1 = BoxLength1/pow(NParticles1,1.0/3.0);
	double ParticleLength2 = BoxLength2/pow(NParticles2,1.0/3.0);


	double ParticlesPerSide1 = pow(NParticles1, 1.0/3.0);	
	double ParticlesPerSide2 = pow(NParticles2, 1.0/3.0);


	for (int x1 = 0; x1 < ParticlesPerSide1; x1++)
	{	
		for (int y1 = 0; y1 < ParticlesPerSide1; y1++)
		{
			for (int z1 = 0; z1 < ParticlesPerSide1; z1++) 
			{
					position.x = x1 * ParticleLength1;
					position.y = y1 * ParticleLength1;
					position.z = z1 * ParticleLength1;
					position.n = 1;
			
					positionVector[particles] = position;
			
					//cout << particles << ": " << position.x << "," << position.y << "," << position.z << "," << "in box " << position.n << endl;
					particles = particles + 1;
			}
		}
	}
				
	for (int x2 = 0; x2 < ParticlesPerSide2; x2++)
	{	
		for (int y2 = 0; y2 < ParticlesPerSide2; y2++)
		{
			for (int z2 = 0; z2 < ParticlesPerSide2; z2++) 
			{
					position.x = x2 * ParticleLength2;
					position.y = y2 * ParticleLength2;
					position.z = z2 * ParticleLength2;
					position.n = 2;
			
					positionVector[particles] = position;
					
					//cout << particles << ": " << position.x << "," << position.y << "," << position.z << "," << "in box " << position.n << endl;
					particles = particles + 1;
			}
		}
	}			

}



void particle::moveParticle (double a, double b, double c, double rmax, int randparticle)
{
		positionVector[randparticle].x = position.x + (2*a - 1)*rmax;
		positionVector[randparticle].y = position.y + (2*b - 1)*rmax;
		positionVector[randparticle].z = position.z + (2*c - 1)*rmax;
}


void particle::restoreParticleMove (double4 OldCoord, int randparticle)
{
		positionVector[randparticle].x = OldCoord.x;
		positionVector[randparticle].y = OldCoord.y;
		positionVector[randparticle].z = OldCoord.z;
		positionVector[randparticle].n = OldCoord.n;
}


double4 particle::getPosition(int randparticle)
{

		position.x = positionVector[randparticle].x;
		position.y = positionVector[randparticle].y;
		position.z = positionVector[randparticle].z;
		position.n = positionVector[randparticle].n;
	
		return position;
} 


vector<double4> particle::getVector()
{
		return positionVector;
}

void particle::rescaleBox(int whichbox, double boxChange, int NTotalParticles)
{
	for (int i = 0; i < NTotalParticles; i++) 
	{
		if (positionVector[i].n == whichbox)
		{
            //cout << "old position x: " << positionVector1[i].x << endl;
            //cout << "old position y: " << positionVector1[i].y << endl;
            //cout << "old position z: " << positionVector1[i].z << endl;
            //cout << "newBL/oldBL: " << boxChange << endl;
            positionVector[i].x = positionVector[i].x*(boxChange);
            positionVector[i].y = positionVector[i].y*(boxChange);
            positionVector[i].z = positionVector[i].z*(boxChange);
            //cout << "new position x: " << positionVector1[i].x << endl;
            //cout << "new position y: " << positionVector1[i].y << endl;
            //cout << "new position z: " << positionVector1[i].z << endl;
    	}
	}
}


void particle::restoreBox(vector<double4> OldVector, int NTotalParticles)
{
	for (int i = 0; i < NTotalParticles; i++)
    {		
        positionVector[i].x =  OldVector[i].x;
        positionVector[i].y =  OldVector[i].y;
        positionVector[i].z =  OldVector[i].z;
	}
}

void particle::swapParticle(int randparticle, double randx, double randy, double randz, int newBox)
{
	positionVector[randparticle].x = randx;
	positionVector[randparticle].y = randy;
	positionVector[randparticle].z = randz;
	positionVector[randparticle].n = newBox;
}






