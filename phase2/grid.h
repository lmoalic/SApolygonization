#ifndef GRID_H
#define GRID_H

#include<list>
#include<vector>

#include"polygon.h"

using namespace std;

class polygon;

class mesh{
    public:
        mesh(){};
        virtual ~mesh(){};


        list<int> lPoints;
        list<int> lCote;
};


class grid
{
    int pdmBrowseX, pdmBrowseY;
    long nextBrowseX, nextBrowseY;
    double alphaBrowseX, alphaBrowseY;
    point p0Browse, p1Browse;

    public:
        grid(){nbMicroSecInitGrid=nbMicroSecUpdateGrid=0;};
        virtual ~grid(){};

        void init(polygon& p, int nbNodePerCell);
        void resetCote();
        void addCote(point& p0, point& p1);
        void updateCote(point& p);
        void initBrowse(point& p0, point& p1);
        void nextMeshes(int& meshX, int& meshY);

        long long xMin, yMin, xMax, yMax;
        int nbMeshX, nbMeshY;
        int pdm;

        double nbMicroSecInitGrid, nbMicroSecUpdateGrid; // pour analyse

        vector<vector<mesh>> vvMeshes;
        vector<vector<pair<list<int>*, list<int>::iterator>>> vvAssoCoteMesh;

    private:
        polygon* poly;
};

#endif // GRID_H
