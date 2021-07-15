#ifndef POLYGON_H
#define POLYGON_H

#include<string>
#include<vector>

#include "point.h"
#include "grid.h"

using namespace std;
class grid;
class polygon
{
    vector<bool> vVisitedEdges;
public:
    string filename, graphName, dirName, filenameSolution;
    vector<point> vPoints;
    vector<int> vNext; // autant d'élément que de points, avec l'id du noeud suivant
    vector<int> vPrev; // autant d'élément que de points, avec l'id du noeud précédent

    grid* gr;

    point xMinMin, xMinMax, xMaxMin, xMaxMax;
    point yMinMin, yMinMax, yMaxMin, yMaxMax;


    // pour analyse
//    vector<unsigned long long> vFreqVisibleEdges;
//    vector<unsigned long long> vFreqTestedEdges;
    long nbTotalIter;
    long nbPrevNextVisible;
    long nbSelectedEdgeNew;
    long nbNewBest;

    int nbCutTested;
    int nbRayTestForAllVisible;



    polygon();
    void load(string filename, int nbNodePerCell=10);
    void loadSolution(string _filename);
    void setSequence(vector<int>& vSeq);
    void save();
    void init(bool maximisation);
    void buildInitSolXincrease();
    bool buildInitSolRand();
    bool buildInitSolRandTry();
    int getNotUsedNode();
    void insertPointAfter(int p0Id, int pPredId);
    bool demidroiteCoupe(point& p, point& p1, point& p2);
    bool estDansPolygone(point& p, point& startingPoly);
    void randomizeSol();
    void determineBounds();
    void determineBounds(vector<point>& vPt, point& ptxMinMin, point& ptxMinMax, point& ptxMaxMin, point& ptxMaxMax, point& ptyMinMin, point& ptyMinMax, point& ptyMaxMin, point& ptyMaxMax);
    long long getArea();
    long long getArea(point& pt0, point& pt1, point& pt2);
    void prepareMin();
    void prepareMax();
    bool isInTriangle(point& p1, point& p2, point& p3, point& pTest);
    bool isTriangleEmpty(point& p1, point& p2, point& p3);
    bool cuttingEdges(point& p0, point& p1,point& p2, point& p3);
    bool cuttingLine(point& p0, point& p1,point& p2, point& p3);
    point getIntersection(point& p0, point& p1, point& p2, point& p3);
    double getDistToIntersection(point& p0, point& p1, point& p2, point& p3);
    pair<int, double> getFirstCrossedEdge(point& p0, point& p1);
    int getFirstCrossedEdgeAfterP1(point& p0, point& p1, int orientation);
    void getVisibleSides(point& p0, vector<int>& vVisible);
    bool getSectorVisibleSides(point& p0, vector<int>& vVisible, point& startPoint, vector<int>& vNextUsed, int orientation, bool startIsVisible);

    void recuitWheelSolver(bool maximization, double firstTemp, double lastTemp, double increTemp, int nbIterPerTemp);
    pair<int, long> getCote(int nodeId, vector<int> vCote, double power);
};

#endif // POLYGON_H
