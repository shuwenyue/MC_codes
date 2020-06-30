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

void particle::setPosition(double BoxLength, int NParticles)
{
	int particles = 0;

	positionVector.resize(NParticles);

	double ParticleLength = BoxLength/pow(NParticles,1.0/3.0);

	double ParticlesPerSide = pow(NParticles, 1.0/3.0);	

	for (int x = 0; x < ParticlesPerSide; x++)
	{
		for (int y = 0; y < ParticlesPerSide; y++)
		{
			for (int z = 0; z < ParticlesPerSide; z++) 
			{
				position.x = x * ParticleLength;
				position.y = y * ParticleLength;
				position.z = z * ParticleLength;

				positionVector[particles] = position;

				// print out coordinates
				//cout << particles << ": " << position.x << "," << position.y << "," << position.z << endl;

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

void particle::restoreParticleMove (double3 OldCoord, int randparticle)
{
	positionVector[randparticle].x = OldCoord.x;
	positionVector[randparticle].y = OldCoord.y;
	positionVector[randparticle].z = OldCoord.z;
}

double3 particle::getPosition(int randparticle)
{
	position.x = positionVector[randparticle].x;
    position.y = positionVector[randparticle].y;
    position.z = positionVector[randparticle].z;	

	return position;
} 

vector<double3> particle::getVector()
{
	return positionVector;
}

void particle::addParticle(double randx, double randy, double randz)
{
	position.x = randx;
	position.y = randy;
	position.z = randz;

	positionVector.push_back(position);
}

void particle::deleteParticle(int randparticle, int NParticles)
{
	positionVector[randparticle].x = positionVector[NParticles-1].x;
	positionVector[randparticle].y = positionVector[NParticles-1].y;
	positionVector[randparticle].z = positionVector[NParticles-1].z;
	
	positionVector.pop_back();
}






