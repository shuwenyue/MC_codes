#include <iostream>
#include <vector>
#include "particle.h"
#include "VectorMath.h"
#include "energy.h"
#include "box.h"

using namespace std;

particle::particle(int NParticles, double BoxLength)
{
	NParticles_ = NParticles;
	BoxLength_ = BoxLength;
}

void particle::setPosition()
{
	int particles = 0; 
	
	double ParticleLength = BoxLength_/pow(NParticles_,1.0/3.0);

	positionVector.resize(NParticles_);

	double ParticlesPerSide = pow(NParticles_, 1.0/3.0);	

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

				particles = particles + 1;

			}
		}
	}
	
}

vector<double3> particle::moveParticle (double a, double b, double c, double rmax, int randparticle)
{

	positionVector[randparticle].x = position.x + (2*a - 1)*rmax;
	positionVector[randparticle].y = position.y + (2*b - 1)*rmax;
	positionVector[randparticle].z = position.z + (2*c - 1)*rmax;

	return positionVector;
}

vector<double3> particle::restoreOld (double3 OldCoord, int randparticle)
{
	positionVector[randparticle].x = OldCoord.x;
	positionVector[randparticle].y = OldCoord.y;
	positionVector[randparticle].z = OldCoord.z;

	return positionVector;
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

vector<double3> particle::rescaleBox(double boxChange)
{
	for (int i = 0; i < NParticles_; i++)
        {
		//cout << "old position x: " << positionVector[i].x << endl;
		//cout << "old position y: " << positionVector[i].y << endl;
		//cout << "old position z: " << positionVector[i].z << endl;
		//cout << "newBL/oldBL: " << boxChange << endl;
       		positionVector[i].x = positionVector[i].x*(boxChange);
                positionVector[i].y = positionVector[i].y*(boxChange);
                positionVector[i].z = positionVector[i].z*(boxChange);
		//cout << "new position x: " << positionVector[i].x << endl;
		//cout << "new position y: " << positionVector[i].y << endl;
		//cout << "new position z: " << positionVector[i].z << endl;
        }
	return positionVector;
}

vector<double3> particle::restoreBox(vector<double3> OldPosition)
{
	for (int i = 0; i < NParticles_; i++)
        {
                positionVector[i].x = OldPosition[i].x;
                positionVector[i].y = OldPosition[i].y;
                positionVector[i].z = OldPosition[i].z;
	}
	return positionVector;
}



