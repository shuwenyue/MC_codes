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


void histogram::plot(double T, double mu, double BoxLength)
{
   	ifstream inFile;
    inFile.open("Energies.txt");

	int totalbinCount = 0;
	int bin1Count = 0;
	int bin2Count = 0;
	int bin3Count = 0;
	int bin4Count = 0;
	int bin5Count = 0;

	double a,b;

	while (!inFile.eof())
	{
		inFile >> a >> b;
        totalbin.push_back(b);
		totalbinCount = totalbinCount + 1;
	}
	totalbinCount = totalbinCount - 1;
    
	for (int i=0; i < totalbinCount; i++)
	{
		if (totalbin[i] >= 0 && totalbin[i] <= 0.2)
		{
			bin1.push_back(totalbin[i]);
			bin1Count = bin1Count + 1;
		}
		if (totalbin[i] > 0.2 && totalbin[i] <= 0.4)
		{
			bin2.push_back(totalbin[i]);
			bin2Count = bin2Count + 1;
		}
		if (totalbin[i] > 0.4 && totalbin[i] <= 0.6)
		{
			bin3.push_back(totalbin[i]);
			bin3Count = bin3Count + 1;
		}
		if (totalbin[i] > 0.6 && totalbin[i] <= 0.8)
		{
			bin4.push_back(totalbin[i]);
			bin4Count = bin4Count + 1;
		}
		if (totalbin[i] > 0.8 && totalbin[i] <= 1.0)
		{
			bin5.push_back(totalbin[i]);
			bin5Count = bin5Count + 1;
		}
	}
	inFile.close(); 


	ofstream outFile;
	outFile.open("histogram.txt");
	
    double densitywidth = 0.2;
    
    outFile << "   T       mu        width    x- y- zdim" << endl;
    outFile << "   " << T << "       " << mu << "       " << densitywidth << "     " << BoxLength << "  " << BoxLength << "  " << BoxLength << endl;

    int width = 8;
    
	int column1 = 0;
    outFile << setw(width) << "1" << setw(width) << bin1Count << endl;  // ENERGY ??
	for (int a=0; a<bin1Count; a++)
	{
		if (column1 > 8) 
		{
			column1 = 0;
			outFile << endl;
		}

		outFile << setw(width) << bin1[a] << ".";
		column1 = column1 + 1;
	}
	outFile << endl;

    outFile << setw(width) << "2" << setw(width) << bin2Count << endl;  // ENERGY ??
	int column2 = 0;
	for (int a=0; a<bin2Count; a++)
	{
		if (column2 > 8) 
		{
			column2 = 0;
			outFile << endl;
		}
		outFile << setw(width) << bin2[a] << ".";
		column2 = column2 + 1;
	}
	outFile << endl;

    outFile << setw(width) << "3" << setw(width) << bin3Count << endl;  // ENERGY ??
	int column3 = 0;
	for (int a=0; a<bin3Count; a++)
	{
		if (column3 > 8) 
		{
			column3 = 0;
			outFile << endl;
		}
		outFile << setw(width) << bin3[a] << ".";
		column3 = column3 + 1;
	}
	outFile << endl;

    outFile << setw(width) << "4" << setw(width) << bin4Count << endl;  // ENERGY ??
	int column4 = 0;
	for (int a=0; a<bin4Count; a++)
	{
		if (column4 > 8) 
		{
			column4 = 0;
			outFile << endl;
		}
		outFile << setw(width) << bin4[a] << ".";
		column4 = column4 + 1;
	}
	outFile << endl;

    outFile << setw(width) << "5" << setw(width) << bin5Count << endl;  // ENERGY ??
	int column5 = 0;
	for (int a=0; a<bin5Count; a++)
	{
		if (column5 > 8) 
		{
			column5 = 0;
			outFile<< endl;
		}
		outFile << setw(width) << bin5[a] << ".";
		column5 = column5 + 1;
	}
	outFile << endl;

	outFile.close();

}



