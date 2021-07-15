#include<vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "Triangulation.h"

using namespace std;

double PARAM_SA_C_0 = -0.08;
double PARAM_SA_C_CS;   // initialized at the beginning of each step
int PARAM_SA_NB_IT;     // depends on the number of points
int PARAM_SA_NB_CS = 8;
double PARAM_SA_C_EF = 7.0;
double PARAM_SA_C_IT = 7.0;
int PARAM_SA_NB_TR = 24;
double PARAM_SA_C_TR = 0.99;
double PARAM_SA_C_SAR = -3.0;
int PARAM_SA_NB_SC = 100;


int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        cerr << "Usage: " << argv[0] << " instanceFile solutionFile maxOrMin" << endl;
        cerr << "maxOrMin: 0 for minimization, 1 for maximization" << endl;
        return 1;
    }

    string instanceFile=argv[1];
    string solutionFile=argv[2];
    bool minimization=(strtol(argv[3], nullptr, 10)==0);

    Triangulation T;
    int nbPoints = T.initLP(instanceFile);

    if (nbPoints==-1)
    {
        cout << "File not found: " << instanceFile << endl;
        return 1;
    }

    PARAM_SA_NB_IT = static_cast<int>(round(152.6*sqrt(nbPoints+10550)));

    T.areaOptimization(minimization);

    cout << "Best area found: " << T.getArea() << endl;

    T.savePolygon(solutionFile);
}
