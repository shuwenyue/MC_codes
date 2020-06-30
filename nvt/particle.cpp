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

	//BoxLength_ = BoxLength_ + (BoxLength_/((pow(NParticles_, 1.0/3.0) - 1)));
	
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
	
	/*for (int i = 0; i < particles; i++)
	{
		cout << "particle: " << i << endl;
		cout << "x: " << positionVector[i].x << endl;
		cout << "y: " << positionVector[i].y << endl;
		cout << "z: " << positionVector[i].z << endl;
	} */
}

vector<double3> particle::moveParticle (double a, double b, double c, double rmax, int randparticle)
{

        //cout << "old x: " << positionVector[randparticle].x << endl;
        //cout << "old y: " << positionVector[randparticle].y << endl;
        //cout << "old z: " << positionVector[randparticle].z << endl;

	positionVector[randparticle].x = position.x + (2*a - 1)*rmax;
	positionVector[randparticle].y = position.y + (2*b - 1)*rmax;
	positionVector[randparticle].z = position.z + (2*c - 1)*rmax;

	//cout << "new x: " << positionVector[randparticle].x << endl;	
	//cout << "new y: " << positionVector[randparticle].y << endl;	
	//cout << "new z: " << positionVector[randparticle].z << endl;	

	return positionVector;
}

//double3 particle::backInBox()

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


