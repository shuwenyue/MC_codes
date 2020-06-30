#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <fstream>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>

using namespace std;

class histogram
{
    
    public:
        void plot(double T, double mu, double BoxLength);

    protected:
    
    private:
        vector<double> totalbin;
        vector<double> bin1;
        vector<double> bin2;
        vector<double> bin3;
        vector<double> bin4;
        vector<double> bin5;
    
};

#endif
