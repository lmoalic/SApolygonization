#include<math.h>
#include<limits>
#include<iostream>

#include "grid.h"

void grid::init(polygon& p, int nbNodePerCell)
{
    clock_t startTime = clock ();
    poly=&p;
    p.determineBounds();
    xMin = p.xMinMin.x;
    yMin = p.yMinMin.y;
    xMax = p.xMaxMin.x;
    yMax = p.yMaxMin.y;

    // calcul du PDM
    double area = (xMax-xMin)*(yMax-yMin);
    int nbParCase = nbNodePerCell;
    double pdm2 = area/p.vPoints.size() * nbParCase;
    pdm = 1.2*sqrt(pdm2);
    nbMeshX = (xMax - xMin) / pdm + 1;
    nbMeshY = (yMax - yMin) / pdm + 1;

    vvMeshes.resize(nbMeshX);
    for(int i=0; i<nbMeshX; i++)
        vvMeshes[i].resize(nbMeshY);

    vvAssoCoteMesh.resize(p.vPoints.size());

    // ajout des points aux bonnes mailles
    for(auto it=p.vPoints.begin(); it!=p.vPoints.end(); it++){
        vvMeshes[(it->x-xMin)/pdm][(it->y-yMin)/pdm].lPoints.push_back(it->id);
    }
    nbMicroSecInitGrid += (clock()-startTime);
}


void grid::resetCote(){
    for(int i=0; i<nbMeshX; i++){
        for(int j=0; j<nbMeshY; j++){
            vvMeshes[i][j].lCote.clear();
        }
    }

    vvAssoCoteMesh.clear();
    vvAssoCoteMesh.resize(poly->vPoints.size());


//    cout<<"Association des cotes aux mailles..."<<endl;
    for(int i=0; i<poly->vPoints.size(); i++){
        if(poly->vNext[i]>-1)
            addCote(poly->vPoints[i], poly->vPoints[poly->vNext[i]]);
    }
}

void grid::addCote(point& p0, point& p1){

    long nextX = xMin+((long)(p0.x-xMin)/pdm)*pdm+1;
    long nextY = yMin+((long)(p0.y-yMin)/pdm)*pdm+1;
    int pdmX=p1.x>p0.x ? pdm : -pdm;
    int pdmY=p1.y>p0.y ? pdm : -pdm;

    if(pdmX>0)
        nextX+=(pdm-2);
    if(pdmY>0)
        nextY+=(pdm-2);

    double alphaX = (fabs(p1.x-p0.x)>0.00001) ? (nextX-p0.x)/(p1.x-p0.x) : numeric_limits<double>::max();
    double alphaY = (fabs(p1.y-p0.y)>0.00001) ? (nextY-p0.y)/(p1.y-p0.y) : numeric_limits<double>::max();

//    cout<<"       alpha = "<<alphaX<<" / "<<alphaY<<endl;
    list<int>& liste = vvMeshes[(nextX-xMin)/pdm][(nextY-yMin)/pdm].lCote;
    liste.push_front(p0.id);
    vvAssoCoteMesh[p0.id].push_back( pair<list<int>*, list<int>::iterator>(&(liste), liste.begin()) );
    while(alphaX<1 || alphaY<1)
    {
        if(alphaX<alphaY){
            nextX += pdmX;
            alphaX = (nextX-p0.x)/(p1.x-p0.x);
        }
        else{
            nextY += pdmY;
            alphaY = (nextY-p0.y)/(p1.y-p0.y);
        }
//        cout<<"   ajout a la maille "<<(nextX-xMin)/pdm<<" / "<<(nextY-yMin)/pdm<<endl;
//        cout<<"       alpha = "<<alphaX<<" / "<<alphaY<<endl;
        list<int>& liste = vvMeshes[(nextX-xMin)/pdm][(nextY-yMin)/pdm].lCote;
        liste.push_front(p0.id);
        vvAssoCoteMesh[p0.id].push_back( pair<list<int>*, list<int>::iterator>(&(liste), liste.begin()) );
    }

}

void grid::updateCote(point& p){
    clock_t startTime = clock ();
    // suppression des liens précédents
    for(auto it=vvAssoCoteMesh[p.id].begin(); it!=vvAssoCoteMesh[p.id].end(); it++)
        it->first->erase( it->second );
    vvAssoCoteMesh[p.id].clear();

    // ajout des nouveaux liens
    addCote(p, poly->vPoints[poly->vNext[p.id]]);

    nbMicroSecUpdateGrid += (clock()-startTime);
}



void grid::initBrowse(point& p0, point& p1){

    nextBrowseX = xMin+((long)(p0.x-xMin)/pdm)*pdm+1;
    nextBrowseY = yMin+((long)(p0.y-yMin)/pdm)*pdm+1;
    pdmBrowseX=p1.x>p0.x ? pdm : -pdm;
    pdmBrowseY=p1.y>p0.y ? pdm : -pdm;

    if(pdmBrowseX>0)
        nextBrowseX+=(pdm-2);
    if(pdmBrowseY>0)
        nextBrowseY+=(pdm-2);

    alphaBrowseX = (fabs(p1.x-p0.x)>0.0000001) ? (nextBrowseX-p0.x)/(p1.x-p0.x) : numeric_limits<double>::max();
    alphaBrowseY = (fabs(p1.y-p0.y)>0.0000001) ? (nextBrowseY-p0.y)/(p1.y-p0.y) : numeric_limits<double>::max();

    p0Browse = p0;
    p1Browse = p1;
}


void grid::nextMeshes(int& meshX, int& meshY){
    if(alphaBrowseX<alphaBrowseY){
        nextBrowseX += pdmBrowseX;
        alphaBrowseX = (nextBrowseX-p0Browse.x)/(p1Browse.x-p0Browse.x);
        meshX = (nextBrowseX-xMin)/pdm;
    }
    else{
        nextBrowseY += pdmBrowseY;
        alphaBrowseY = (nextBrowseY-p0Browse.y)/(p1Browse.y-p0Browse.y);
        meshY = (nextBrowseY-yMin)/pdm;
    }
}
