#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include "polygon.h"
#include<set>

void printCpuInfo();


using namespace std;

inline bool exists (const std::string& filename) {
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0);
}


int main(int argc, char** argv)
{

    printCpuInfo();
    time_t startingTime=time(NULL);
    polygon poly;

    srand(getpid());

    if(argc<3){
        cout<<"Argc="<<argc<<endl;
        cout<<"Use: polygon filename maximisation[0:minimisation 1:maximisation] nbNodePerCell[default=10]] startingTemp endingTemp incrTemp"<<endl;
        exit(1);
    }
    string filename=argv[1];
    bool maximisation=atoi(argv[2]);

    int nbNodePerCell = (argc>=4) ? atoi(argv[3]) : 10;
    double startingTemp = (argc>=5) ? atof(argv[4]) : 1;
    double endingTemp = (argc>=6) ? atof(argv[5]) : 200;
    double incrTemp = (argc>=7) ? atof(argv[6]) : 1;
    int nbIterPerTemp = (argc>=8) ? atoi(argv[7]) : 100000;



    // chargement du polygone
    poly.load(filename, nbNodePerCell);
    poly.init(maximisation);
    cout<<"Polygon loaded"<<endl;

    // chagement d'une éventuelle solution précédente
    cout<<"Recherche de : "<<"./solution/"+poly.filenameSolution<<endl;
    bool exist = exists("./solution/"+poly.filenameSolution);
    if(exist){
        poly.loadSolution("./solution/"+poly.filenameSolution);
    }
    long long val=poly.getArea();
    cout<<"Aire du polygon au lancement : "<<abs(val)<<endl;

    vector<int> vBestSol=poly.vNext;
    long bestArea = abs(val);
    time_t tsaveConfig=time(NULL);


        poly.recuitWheelSolver(maximisation, startingTemp, endingTemp, incrTemp, nbIterPerTemp);
        val=abs(poly.getArea());
        printCpuInfo();
        cout<<difftime(time(NULL), startingTime)/3600<<" heures => Aire après optimisation (";
        cout<<startingTemp<<" , "<<endingTemp<<" , "<<incrTemp<<") : "<<val<<"     BestArea="<<bestArea<<endl;
        if((val<bestArea && maximisation) || (val>bestArea && !maximisation)){
            poly.setSequence(vBestSol);
        }
        else if(val != bestArea){
            bestArea=val;
            vBestSol=poly.vNext;
            cout<<poly.graphName<<" new bestVal => "<<bestArea<<endl;
        }
        if(difftime(time(NULL), tsaveConfig)>=600){
            tsaveConfig=time(NULL);
            poly.save();
        }


    val=poly.getArea();
    cout<<"Aire après optimisation : "<<abs(val)<<endl;
    poly.save();

    return 0;
}


void printCpuInfo(){
    string line;
    ifstream finfo("/proc/cpuinfo");
    while(getline(finfo,line)) {
        stringstream str(line);
        string itype;
        string info;
        if ( getline( str, itype, ':' ) && getline(str,info) && itype.substr(0,10) == "model name" ) {
            cout << info << endl;
            break;
        }
    }
}
