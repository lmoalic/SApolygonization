#include<iostream>
#include <sstream>
#include <iomanip>
#include<limits>
#include<algorithm>
#include <unistd.h>
#include "polygon.h"
#include "util/gfile.h"


using namespace std;


bool maxPairIntLong(const pair<int,long>& p1, const pair<int,long>& p2) { return p1.second>p2.second; }
bool maxPairIntDouble(const pair<int,double>& p1, const pair<int,double>& p2) { return p1.second>p2.second; }

polygon::polygon()
{
}

void polygon::load(string _filename, int nbNodePerCell){
    filename = _filename;
    int pos=filename.find_last_of("/\\");
    dirName=filename.substr(0,pos);
    graphName=filename.substr(pos+1);

    cout<<"Repertoire : "<<dirName<<endl;
    cout<<"Graph : "<<graphName<<endl;

    ostringstream oss;
    oss<<"./bestSolutions/"<<graphName<<".solution";
    filenameSolution=oss.str();

    GInputFile infile(filename) ;
    char* buf ;

    infile.open();

    if(! infile.isOpen()){
        cout<<"Fichier introuvable: "<<filename<<endl;
        exit(10);
    }

    vPoints.clear();
    vPrev.clear();
    vNext.clear();
    while ((buf=infile.readUncommentedLine())!=NULL) {
        int id = infile.getNextIntToken();
        int x = infile.getNextIntToken();
        int y = infile.getNextIntToken();

        vPoints.push_back(point(id, x, y));
    }

    vPrev.resize(vPoints.size(), -1);
    vNext.resize(vPoints.size(), -1);

    gr=new grid();
    gr->init(*this, nbNodePerCell);

    init(true);


    cout<<"Nombre de points chargés : "<<vPoints.size()<<endl;
}


// chargement de l'ordre des sommets du ploygone a partir d'un fichier
void polygon::loadSolution(string _filename)
{
    cout<<"Chargement de la solution : "<<_filename<<endl;
    ifstream f(_filename);
    int nodeId;
    int first, prev;

    // suprimme les commentaires
    int position=0;
    char c;
    f>>c;
    while(c=='#'){
        string str;
        getline(f, str);
        position = f.tellg();
        f>>c;
    }
    f.seekg(position);
    //

    f >> first >> ws;
    prev = first;


    while (f >> nodeId >> ws){
        vNext[prev]=nodeId;
        vPrev[nodeId]=prev;
        prev=nodeId;
    }
    vNext[prev]=first;
    vPrev[first]=prev;

    gr->resetCote();
}


void polygon::setSequence(vector<int>& vSeq){
    for(int i=0; i<vPoints.size(); i++){
        vNext[i] = vSeq[i];
        vPrev[vSeq[i]]=i;
    }

    gr->resetCote();
}

void polygon::save(){
    cout<<"Sauvegarde de la solution"<<endl;
    ostringstream oss;

    oss<<"bestSolutions/"<<filenameSolution<<getpid();
    FILE *f;
    f = fopen(oss.str().c_str(), "w");

    int current=0;
    for(int i=0; i<vNext.size(); i++){
        fprintf(f, "%d\n", current);
        current = vNext[current];
    }
    fclose(f);
}

void polygon::init(bool maximisation){
    determineBounds();
    buildInitSolXincrease();

    if(maximisation){
        cout<<"Préparation pour maximisation"<<endl;
        prepareMax();
    }
    else{
        cout<<"Préparation pour minimisation"<<endl;
        prepareMin();
    }
}


void polygon::determineBounds(){
    determineBounds(vPoints, xMinMin, xMinMax, xMaxMin, xMaxMax, yMinMin, yMinMax, yMaxMin, yMaxMax);
}

void polygon::determineBounds(vector<point>& vPt, point& ptxMinMin, point& ptxMinMax, point& ptxMaxMin, point& ptxMaxMax, point& ptyMinMin, point& ptyMinMax, point& ptyMaxMin, point& ptyMaxMax){
    // détermination des points extrêmes
    ptxMinMin=ptxMinMax=ptxMaxMin=ptxMaxMax=vPt[0];
    ptyMinMin=ptyMinMax=ptyMaxMin=ptyMaxMax=vPt[0];
    for(auto it=vPt.begin(); it!=vPt.end(); it++){
        if(it->x <= ptxMinMin.x){
            if(it->x < ptxMinMin.x)
                ptxMinMin = ptxMinMax = *it;
            else if(it->y < ptxMinMin.y)
                ptxMinMin = *it;
            else if(it->y > ptxMinMax.y)
                ptxMinMax = *it;
        }
        if(it->x >= ptxMaxMin.x){
            if(it->x > ptxMaxMin.x)
                ptxMaxMin = ptxMaxMax = *it;
            else if(it->y < ptxMaxMin.y)
                ptxMaxMin = *it;
            else if(it->y > ptxMaxMax.y)
                ptxMaxMax = *it;
        }
        if(it->y <= ptyMinMin.y){
            if(it->y < ptyMinMin.y)
                ptyMinMin = ptyMinMax = *it;
            else if(it->x < ptyMinMin.x)
                ptyMinMin = *it;
            else if(it->x > ptyMinMax.x)
                ptyMinMax = *it;
        }
        if(it->y >= ptyMaxMin.y){
            if(it->y > ptyMaxMin.y)
                ptyMaxMin = ptyMaxMax = *it;
            else if(it->x < ptyMaxMin.x)
                ptyMaxMin = *it;
            else if(it->x > ptyMaxMax.x)
                ptyMaxMax = *it;
        }
    }

}



void polygon::buildInitSolXincrease(){
    int currentUp, currentDown;
    currentUp=currentDown=vPoints[0].id;

    // droite coupant les points en 2 moitiés : (ptxMinMin, ptxMaxMax)
    // y=ax+b
    double a = (xMaxMax.y - xMinMin.y)/((double) xMaxMax.x - xMinMin.x);
    double b = -a * xMinMin.x + xMinMin.y;

    auto it=vPoints.begin();
    it++;
    for(it; it!=vPoints.end(); it++){
        if(it->y > (a*it->x+b)){
            vPrev[currentUp] = it->id;
            vNext[it->id] = currentUp;
            currentUp = it->id;
        }
        else{
            vNext[currentDown] = it->id;
            vPrev[it->id] = currentDown;
            currentDown = it->id;
        }
    }
    vPrev[currentUp] = currentDown;
    vNext[currentDown] = currentUp;

    gr->resetCote();
}

bool polygon::buildInitSolRand(){
    while(!buildInitSolRandTry()){

    }
}
bool polygon::buildInitSolRandTry(){
    vPrev.clear();
    vNext.clear();
    vPrev.resize(vPoints.size(), -1);
    vNext.resize(vPoints.size(), -1);
    vector<bool> vUsedNode(vPoints.size(), false);

    int id1=getNotUsedNode();
    int id2=getNotUsedNode();
    int id3=getNotUsedNode();

    while(id2==id1)
        id2=getNotUsedNode();

    while(id3==id2 || id3==id1 || vPoints[id1].left(vPoints[id2],vPoints[id3])==0)
        id3=getNotUsedNode();

    if(vPoints[id1].left(vPoints[id2],vPoints[id3])<0){
        int tmp=id1;
        id1=id2;
        id2=tmp;
    }

    vPrev[id1]=id3;
    vNext[id1]=id2;
    vPrev[id2]=id1;
    vNext[id2]=id3;
    vPrev[id3]=id2;
    vNext[id3]=id1;

    gr->resetCote();

    int nbAdded=3;
    int idToAdd=-1;
    int nbEchec=0;
    while(nbAdded!=vPoints.size() && nbEchec<100){
        vector<int> vPointsVisible;

        idToAdd=getNotUsedNode();

        vVisitedEdges.clear();
        vVisitedEdges.resize(vPoints.size(), false);
        int orientation = estDansPolygone(vPoints[idToAdd], vPoints[id1]) ? 1 : -1;
        bool complete=getSectorVisibleSides(vPoints[idToAdd], vPointsVisible, vPoints[id1], vNext, orientation, false);
        if(!complete)
            getSectorVisibleSides(vPoints[idToAdd], vPointsVisible, vPoints[id1], vPrev, -1*orientation, false);
        if(vPointsVisible.size()>0){
            insertPointAfter(idToAdd, vPointsVisible[0]);
            //printSolution();
            nbAdded++;
            nbEchec=0;
        }
        else
            nbEchec++;
    }
    return nbEchec<100;
}


int polygon::getNotUsedNode(){
    int id=rand()/((double)RAND_MAX+1) * vPrev.size();

    while(vPrev[id]>-1){
        id++;
        if(id==vPrev.size())
            id=0;
    }
    return id;
}

void polygon::insertPointAfter(int p0Id, int pPredId){
    int next=vNext[pPredId];
    vNext[pPredId]=p0Id;
    vPrev[next]=p0Id;
    vNext[p0Id]=next;
    vPrev[p0Id]=pPredId;

    gr->updateCote(vPoints[pPredId]);
    gr->updateCote(vPoints[p0Id]);
}

//demie droite horizontale issue de p (vers la droite)
//coupe strictement segment [p1,p2]
//p à gauche de (p1,p2) ssi determinant(p1,p2,p)>0

bool polygon::demidroiteCoupe(point& p, point& p1, point& p2)
{
  return ( (p1.y<p.y && p2.y>p.y) && p.left(p1,p2)==1 )
     ||  ( (p1.y>p.y && p2.y<p.y) && p.left(p1,p2)==-1 );
}

//point dans polygone si demie-droite issue coupe nombre impair de cotes
bool polygon::estDansPolygone(point& p, point& startingPoly)
{
  int nbcoupe=0;
  int courant = startingPoly.id;
  do
  {
    if (demidroiteCoupe(p, vPoints[courant], vPoints[vNext[courant]]))
        nbcoupe++;
    courant = vNext[courant];
  }while (courant != startingPoly.id);
  return nbcoupe%2==1;
}


void polygon::randomizeSol(){
    for(int i=0; i<10000; i++){
        int nodeId=rand()/((double)RAND_MAX+1) * vPoints.size();
        if(isTriangleEmpty(vPoints[vPrev[nodeId]], vPoints[nodeId], vPoints[vNext[nodeId]])){
            vector<int> vCote;
            getVisibleSides(vPoints[nodeId], vCote);

            if(vCote.size()>0){
                int bestCote=vCote[rand()/((double)RAND_MAX+1) * vCote.size()];
                vPrev[vNext[nodeId]] = vPrev[nodeId];
                vNext[vPrev[nodeId]] = vNext[nodeId];
                gr->updateCote(vPoints[vPrev[nodeId]]);

                vNext[nodeId] = vNext[bestCote];
                vPrev[nodeId] = bestCote;
                gr->updateCote(vPoints[nodeId]);

                vPrev[vNext[bestCote]] = nodeId;
                vNext[bestCote] = nodeId;
                gr->updateCote(vPoints[bestCote]);
            }
        }
    }
}

long long polygon::getArea(){
    int xCentre = (xMinMin.x + xMaxMax.x)/2;
    int yCentre = (yMinMin.y + yMaxMax.y)/2;

    // Ajoute le côté du dernier point au premier
    long long x1 = (vPoints[vPrev[0]].x-xCentre);
    long long y1 = (vPoints[vPrev[0]].y-yCentre);
    long long x2, y2;

    long long area=0;
    int current=0;
    for(int i=0; i<vNext.size(); i++){
        x2 = (vPoints[current].x-xCentre);
        y2 = (vPoints[current].y-yCentre);
        area+=( x1*y2  -  y1*x2);
        x1 = x2;
        y1 = y2;
        current=vNext[current];
    }

    return area/2;
}


long long polygon::getArea(point& pt0, point& pt1, point& pt2){
    long long x1 = (pt1.x-pt0.x);
    long long y1 = (pt1.y-pt0.y);
    long long x2 = (pt2.x-pt0.x);
    long long y2 = (pt2.y-pt0.y);

    long long area = x1*y2  -  y1*x2;

    return area/2;
}

void polygon::prepareMax(){
    if(vPoints[vPrev[0]].left(vPoints[0], vPoints[vNext[0]]) < 0)
        vPrev.swap(vNext);
    gr->resetCote();

    int pos=graphName.find(".instance");
    graphName=graphName.substr(0,pos);

    ostringstream oss;
    oss<<"./"<<graphName<<".max.solution" ;
    filenameSolution=oss.str();
}


void polygon::prepareMin(){
    if(vPoints[vPrev[0]].left(vPoints[0], vPoints[vNext[0]]) > 0)
        vPrev.swap(vNext);
    gr->resetCote();

    int pos=graphName.find(".instance");
    graphName=graphName.substr(0,pos);

    ostringstream oss;
    oss<<"./"<<graphName<<".min.solution" ;
    filenameSolution=oss.str();
}



bool polygon::isInTriangle(point& p1, point& p2, point& p3, point& pTest){
    bool res=false;
    int xMin=p1.x, xMax=p1.x, yMin=p1.y, yMax=p1.y;

    if(p2.x<xMin) xMin=p2.x;
    else if(p2.x>xMax) xMax=p2.x;
    if(p3.x<xMin) xMin=p3.x;
    else if(p3.x>xMax) xMax=p3.x;

    if(p2.y<yMin) yMin=p2.y;
    else if(p2.y>yMax) yMax=p2.y;
    if(p3.y<yMin) yMin=p3.y;
    else if(p3.y>yMax) yMax=p3.y;


    if(pTest.x>=xMin && pTest.x<=xMax && pTest.y>=yMin && pTest.y<=yMax){
        long long ux=p1.x-pTest.x;
        long long uy=p1.y-pTest.y;
        long long vx=p2.x-pTest.x;
        long long vy=p2.y-pTest.y;
        long long wx=p3.x-pTest.x;
        long long wy=p3.y-pTest.y;

        long long uv = ux*vy-uy*vx;
        long long vw = vx*wy-vy*wx;
        long long wu = wx*uy-wy*ux;

        if(uv==0)
            res = (vw>0 == wu>0);
        else if(vw==0)
            res = (uv>0 == wu>0);
        else if(wu==0)
            res = (uv>0 == vw>0);
        else
            res = (uv>0 == vw>0 && vw>0 == wu>0);
    }
    return res;
}



bool polygon::isTriangleEmpty(point& p1, point& p2, point& p3){
    int pdm=gr->pdm;
    int mXmin = (p1.x-gr->xMin)/pdm;
    int mYmin = (p1.y-gr->yMin)/pdm;
    int mXmax=mXmin;
    int mYmax=mYmin;
    int tmp;
    if((tmp=(p2.x-gr->xMin)/pdm)<mXmin) mXmin=tmp;
    else if(tmp>mXmax)mXmax=tmp;
    if((tmp=(p2.y-gr->yMin)/pdm)<mYmin) mYmin=tmp;
    else if(tmp>mYmax)mYmax=tmp;
    if((tmp=(p3.x-gr->xMin)/pdm)<mXmin) mXmin=tmp;
    else if(tmp>mXmax)mXmax=tmp;
    if((tmp=(p3.y-gr->yMin)/pdm)<mYmin) mYmin=tmp;
    else if(tmp>mYmax)mYmax=tmp;

    for(int i=mXmin; i<=mXmax; i++){
        for(int j=mYmin; j<=mYmax; j++){
            for(auto it=gr->vvMeshes[i][j].lPoints.begin(); it!=gr->vvMeshes[i][j].lPoints.end(); it++){
                if(*it != p1.id && *it != p2.id && *it != p3.id){
                    if(isInTriangle(p1, p2, p3, vPoints[*it]))
                        return false;
                }
            }
        }
    }

    return true;
}

// determine if p2p3 cut the edge p0p1.
bool polygon::cuttingEdges(point& p0, point& p1,point& p2, point& p3){
    nbCutTested++;
    int prod1=p0.left(p2,p3)*p1.left(p2,p3);
    int prod2=p2.left(p0,p1)*p3.left(p0,p1);
    bool res=prod1<0 && prod2<0;
    if(prod1 == 0 || prod2 == 0){ // if at least 3 points are aligned
        if(prod1 == 0 && prod2 == 0) // the 4 points are aligned
            res = (p0.prodScal(p3,p1,p3)<0); // true if p2p3 is included in p0p1
        else
            res = (prod1+prod2 < 0); // true if p2 or p3 is aligned with and between p0 and p1
    }
    return res;
}

// determine if p2p3 cut the half line p0p1.
bool polygon::cuttingLine(point& p0, point& p1,point& p2, point& p3){
    nbCutTested++;
//    cout<<"\t\tintersection de p"<<p0<<p1<<" avec "<<p2<<p3<<" : "<<p2.left(p0,p1)<<"*"<<p3.left(p0,p1)<<"="<<p2.left(p0,p1)*p3.left(p0,p1)<<endl;
    if(p2.left(p0,p1)*p3.left(p0,p1)<=0 && p0.prodScal(p1, p0, getIntersection(p0,p1,p2,p3))>0 )
        return true;
    return false;
}

// determine the intersection between p0p1 and p2p3
point polygon::getIntersection(point& p0, point& p1, point& p2, point& p3){
    point I(-1, p1.x - p0.x, p1.y - p0.y);
    point J(-1, p3.x - p2.x, p3.y - p2.y);
    double m=0;
    double diviseur = (I.x * J.y - I.y * J.x);
    point res;
    if(diviseur != 0)
    {
        m = (I.x * p0.y
             - I.x * p2.y
             - I.y * p0.x
             + I.y * p2.x
            ) / diviseur;
        res.x = p2.x + m * J.x;
        res.y = p2.y + m * J.y;
    }
    else{
        double dist = p0.dist(p2);
        double d2 = p0.dist(p3);
        if(d2<dist)
            dist=d2;
        double alpha = dist/p0.dist(p1);
        res.x = p0.x + alpha * (p1.x-p0.x);
        res.y = p0.y + alpha * (p1.y-p0.y);
    }
    return res;
}

// determine the distance between p0 and the intersection between p0p1 and p2p3
double polygon::getDistToIntersection(point& p0, point& p1, point& p2, point& p3){
    return p0.dist(getIntersection(p0, p1,p2,p3));
}

// renvoie la paire <id, distance> du point source du côté le plus proche de p0 dans la direction (p0,p1)
pair<int, double> polygon::getFirstCrossedEdge(point& p0, point& p1){
    int bestId=p1.id;
    double bestDist=p0.dist(p1);
    int mailleXbest=-1, mailleYbest=-1;

    vector<bool> vVisitedSide(vPoints.size(), false);
    vVisitedSide[p0.id]=true;

    int mailleX = (p0.x-gr->xMin)/gr->pdm;
    int mailleY = (p0.y-gr->yMin)/gr->pdm;

    gr->initBrowse(p0, p1);

    bool finished=false;

    // check if [p0.pred ; p0.next] is cut
    if(cuttingEdges(p0, p1, vPoints[vPrev[p0.id]], vPoints[vNext[p0.id]])){
        bestDist = getDistToIntersection(p0, p1, vPoints[vPrev[p0.id]], vPoints[vNext[p0.id]]);
        bestId = vPrev[p0.id];
    }

    nbCutTested=0;
    while( !finished ){
        // >>recherche du côtés le plus proche (non encore traités)
        list<int>& liste = gr->vvMeshes[mailleX][mailleY].lCote;
        int nbElt = liste.size();

        for(auto& pt : liste){
            if(!vVisitedSide[pt]){
                if(cuttingEdges(p0, p1, vPoints[pt], vPoints[vNext[pt]])){
                    double dist = getDistToIntersection(p0, p1, vPoints[pt], vPoints[vNext[pt]]);
                    if(dist<bestDist){
                        bestDist = dist;
                        bestId = pt;
                    }
                }
                vVisitedSide[pt]=true;
            }
        }
        // <<recherche du côtés le plus proche (non encore traités)


        // >> détermination de la maille contenant l'intersection entre la direction et l'actuel bestId
        if(bestId>-1){
            point bestPoint = getIntersection(p0, p1, vPoints[bestId], vPoints[vNext[bestId]]);
            mailleXbest = (bestPoint.x-gr->xMin)/gr->pdm;
            mailleYbest = (bestPoint.y-gr->yMin)/gr->pdm;
        }
        // << détermination de la maille contenant l'intersection entre la direction et l'actuel bestId

        if(bestId>-1 && mailleXbest==mailleX && mailleYbest==mailleY)
            finished=true;
        else{ // mise à jour de la maille à explorer
            gr->nextMeshes(mailleX, mailleY);
        }
        if(mailleX<0 || mailleX>=gr->nbMeshX || mailleY<0 || mailleY>=gr->nbMeshY){
            //cout<<"ERREUR : pas d'intersection avec "<<p0<<" et "<<p1<<endl;
            finished=true;
        }
    }
//    vFreqTestedEdges[nbCutTested]++;
    return pair<int, double>(bestId, bestDist);
}



// return the point id of the cutted edge after p1 with one side in the orientation
int polygon::getFirstCrossedEdgeAfterP1(point& p0, point& p1, int orientation){
    int bestId=-1;
    double bestDist=numeric_limits<double>::max();
    int mailleXbest=-1, mailleYbest=-1;

    vector<bool> vVisitedSide(vPoints.size(), false);
    vVisitedSide[p0.id]=true;
    if(vPrev[p0.id]>-1){
        vVisitedSide[vPrev[p0.id]]=true;
        vVisitedSide[vNext[p0.id]]=true;
    }

    int mailleX = (p0.x-gr->xMin)/gr->pdm;
    int mailleY = (p0.y-gr->yMin)/gr->pdm;

    gr->initBrowse(p0, p1);

    bool finished=false;

    nbCutTested=0;
    while( !finished ){
        // >>recherche du côtés le plus proche (non encore traités)
        list<int>& liste = gr->vvMeshes[mailleX][mailleY].lCote;

        for(auto& pt : liste){
            if(!vVisitedSide[pt] && pt!=p1.id && vNext[pt]!=p1.id){
//                cout<<"\t\t\t essai de "<<vPoints[pt]<<"  :  "<<cuttingLine(p0, p1, vPoints[pt], vPoints[vNext[pt]])<<endl;
                if(cuttingLine(p0, p1, vPoints[pt], vPoints[vNext[pt]]) && (vPoints[pt].left(p0,p1) == orientation || vPoints[vNext[pt]].left(p0,p1) == orientation)){
                    double dist = getDistToIntersection(p0, p1, vPoints[pt], vPoints[vNext[pt]]);
                    if(dist<bestDist){
                        bestDist = dist;
                        bestId = pt;
                    }
                }
                vVisitedSide[pt]=true;
            }
        }
        // <<recherche du côtés le plus proche (non encore traités)


        // >> détermination de la maille contenant l'intersection entre la direction et l'actuel bestId
        if(bestId>-1){
            point bestPoint = getIntersection(p0, p1, vPoints[bestId], vPoints[vNext[bestId]]);
            mailleXbest = (bestPoint.x-gr->xMin)/gr->pdm;
            mailleYbest = (bestPoint.y-gr->yMin)/gr->pdm;
        }
        // << détermination de la maille contenant l'intersection entre la direction et l'actuel bestId

        if(bestId>-1 && mailleXbest==mailleX && mailleYbest==mailleY)
            finished=true;
        else{ // mise à jour de la maille à explorer
            gr->nextMeshes(mailleX, mailleY);
            if(mailleX<0 || mailleY<0 || mailleX>=gr->nbMeshX || mailleY>=gr->nbMeshY)
                finished=true;
        }
    }
    return bestId;
}

void polygon::getVisibleSides(point& p0, vector<int>& vVisible){
    vVisible.clear();
    int orientation = p0.left(vPoints[vPrev[p0.id]], vPoints[vNext[p0.id]]);
    if(orientation == 0) // p0 ne peut prendre une autre place (il est aligné, entre son pred et son suivant)
        return;
    point& next = vPoints[vNext[p0.id]];
    point& prev = vPoints[vPrev[p0.id]];
    vNext[prev.id]=next.id;
    vPrev[next.id]=prev.id;




    vVisitedEdges.clear();
    vVisitedEdges.resize(vPoints.size(), false);
    bool complete=getSectorVisibleSides(p0, vVisible, next, vNext, orientation, true);
    if(!complete)
        getSectorVisibleSides(p0, vVisible, prev, vPrev, -1*orientation, true);

    vNext[prev.id]=p0.id;
    vPrev[next.id]=p0.id;
}


bool polygon::getSectorVisibleSides(point& p0, vector<int>& vVisible, point& startPoint, vector<int>& vNextUsed, int orientation, bool startIsVisible){
    point* current=&vPoints[vNextUsed[startPoint.id]];
    point* pred=&startPoint;

    double angleParcouru=0;
    bool finished = false;
    point* currentOld = pred; // a supprimer => pour angle

//    cout<<"Recherche des cotes visibles de "<<p0<<endl;
    while(!vVisitedEdges[current->id] && !finished){
        vVisitedEdges[current->id]=true;
//        cout<<"\t current est "<<*current<<endl;
        if(current->left(p0, *pred) != orientation){//current go back
//            cout<<"\t\t going back => recherche du suivant allant dans le bon sens"<<endl;
            int edgeId = getFirstCrossedEdgeAfterP1(p0, *pred, orientation);
//            cout<<"\t\t trouve p"<<vPoints[edgeId]<<endl;
            if(edgeId<0)
                finished=true;
            else{
                pred = (&vNextUsed == &vNext) ? &vPoints[edgeId] : &vPoints[vNext[edgeId]];
                current = &vPoints[vNextUsed[pred->id]];
                startIsVisible = false;
            }
        }
        else{
            angleParcouru+=p0.getAngle(*currentOld, *current);// a supprimer =>  pour angle
            currentOld=current;
            pair<int, double> plusProche = getFirstCrossedEdge(p0, *current);

            if( p0.dist(*current)<=plusProche.second){ // current is visible
//                cout<<"\t\t"<<*current<<" is visible from "<<p0<<endl;
                if(startIsVisible){
//                    cout<<"Ajout de "<<*pred<<endl;
                    if(&vNextUsed==&vNext){
                        if(pred->id!=vPrev[p0.id])
                            vVisible.push_back(pred->id);
                    }
                    else if(current->id!=vPrev[p0.id])
                        vVisible.push_back(current->id);
                }
                pred = current;
                current = &vPoints[vNextUsed[current->id]];
                startIsVisible = true;
            }
            else{ //current is not visible
                // go in the other direction until lastVisible
                pred = (&vNextUsed == &vNext) ? &vPoints[plusProche.first] : &vPoints[vNext[plusProche.first]];
                pair<int, double> plusProche = getFirstCrossedEdge(p0, *pred);
                if(p0.dist(*pred)<=plusProche.second)
                    startIsVisible = true;
                else
                    startIsVisible = false;
//                cout<<"***>> check other direction"<<endl;
                getSectorVisibleSides(p0, vVisible, *pred, (&vNextUsed == &vNext)?vPrev:vNext, -1*orientation, startIsVisible);
//                cout<<"***<< check other direction"<<endl;
                current = &vPoints[vNextUsed[pred->id]]; // à affiner, repartir en arrière ??
            }
        }
    }
//    cout<<"\t\tangle parcouru:"<<angleParcouru<<endl;
    return !finished;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void polygon::recuitWheelSolver(bool maximization, double firstTemp, double lastTemp, double increTemp, int nbIterPerTemp){
    //cout<<"maxim="<<maximization<<" firstTemp="<<firstTemp<<" lastTemp="<<lastTemp<<" incr="<<increTemp<<" nbIter="<<nbIterPerTemp<<endl;
    nbTotalIter=0;
    nbPrevNextVisible=0;
    nbSelectedEdgeNew=0;
    nbNewBest=0;



    if(maximization)
        prepareMax();
    else
        prepareMin();

    long startingArea = getArea();
    long bestArea = startingArea;
    vector<int> vBestNext = vNext;
    long improve=0;

    time_t startingTime=time(NULL);
    bool heure1=false;
    bool heure10=false;

        for(double power = firstTemp; power<=lastTemp; power+=increTemp){
            setSequence(vBestNext);
            improve=0;

            //////// pour afficher l'avancement
            double hourpassed = difftime(time(NULL), startingTime)/3600;
            if((hourpassed>1 && !heure1) || (hourpassed>10 && !heure10)){
                cout<<"Record"<<(int)hourpassed<<"h ("<<getpid()<<") & "<<graphName<<" & "<<bestArea<<" & "<<nbTotalIter<<" & "<<nbPrevNextVisible<<" & "<<nbSelectedEdgeNew<<" & "<<nbNewBest<<" \\tabularnewline"<< endl;
                if(!heure1) heure1=true; else heure10=true;
                save();
            }
            ////////

            for(int i=0; i<nbIterPerTemp; i++){
                nbTotalIter++;
                int nodeId=rand()/((double)RAND_MAX+1) * vPoints.size();

                if(isTriangleEmpty(vPoints[vPrev[nodeId]], vPoints[nodeId], vPoints[vNext[nodeId]])){
                    nbPrevNextVisible++;
                    vector<int> vCote;
                    getVisibleSides(vPoints[nodeId], vCote);
//                    vFreqVisibleEdges[vCote.size()]++;
                    pair<int,long> selectedEdge=getCote(nodeId, vCote, power);
                    int bestCote=selectedEdge.first;
                    improve+=(selectedEdge.second);

                    if(bestCote!=vPrev[nodeId]){
                        nbSelectedEdgeNew++;
    //                        cout<<"Prev="<<vPrev[nodeId]<<" et meilleur="<<bestCote<<endl;
                        vPrev[vNext[nodeId]] = vPrev[nodeId];
                        vNext[vPrev[nodeId]] = vNext[nodeId];
                        gr->updateCote(vPoints[vPrev[nodeId]]);

                        vNext[nodeId] = vNext[bestCote];
                        vPrev[nodeId] = bestCote;
                        gr->updateCote(vPoints[nodeId]);

                        vPrev[vNext[bestCote]] = nodeId;
                        vNext[bestCote] = nodeId;
                        gr->updateCote(vPoints[bestCote]);

                        if(improve>0){
                            nbNewBest++;
                            bestArea+=improve;
                            vBestNext=vNext;
                            improve=0;
                            cout<<difftime(time(NULL), startingTime)/3600<<"h ("<<getpid()<<") : \t\tnouvelle meilleur aire "<<bestArea<<" power="<<power<< endl;
                        }
                    }
                }
            }

//        cout<<"Power at the end : "<<power<<"\t => "<<bestArea<<endl;
    }


    double hourpassed = difftime(time(NULL), startingTime)/3600;
    cout<<"RecordFinal "<<hourpassed<<"h ("<<getpid()<<") & "<<graphName<<" & "<<bestArea<<" & "<<nbTotalIter<<" & "<<nbPrevNextVisible<<" & "<<nbSelectedEdgeNew<<" & "<<nbNewBest<<" \\tabularnewline"<< endl;


    long finalArea=getArea();
    cout<<"\t *** (t0="<<firstTemp<<" tfin="<<lastTemp<<") LANCEMENT="<<startingArea<<"  FIN="<<finalArea<<"  BEST="<<bestArea<<" ("<<(finalArea-startingArea)*100.0/startingArea<<"%)"<<endl;

    cout<<"TotalIter="<<setw(10) << nbTotalIter;
    cout<<"nbPrevNextVisible="<<setw(10)<<nbPrevNextVisible;
    cout<<"nbSelectedEdgeNew="<<setw(10)<<nbSelectedEdgeNew;
    cout<<"nbNewBest="<<setw(10)<<nbNewBest;


    setSequence(vBestNext);

//    cout<<"Nombre d'itérations : "<<nbIter<<endl;
}


pair<int, long> polygon::getCote(int nodeId, vector<int> vCote, double power){
    vector<long> vArea(vCote.size());
    long currentArea = getArea(vPoints[vPrev[nodeId]], vPoints[nodeId], vPoints[vNext[nodeId]]);
    double total= currentArea;

    for(int i=0; i<vCote.size(); i++){
        long area = getArea(vPoints[vCote[i]], vPoints[nodeId], vPoints[vNext[vCote[i]]]);
        vArea[i] = area;
        total+=area;
    }

    vector<pair<int, double>> vCoteArea(vCote.size()+1);
    for(int i=0; i<vCote.size(); i++){
        vCoteArea[i] = pair<int,double>(i, vArea[i]/total);
    }
    vCoteArea[vCote.size()] = pair<int, double>(vCote.size(), currentArea/total);

    // pour la minimisation on inverse les aires
    if(currentArea<0){
        double total2=0;
        for(int i=0; i<vCoteArea.size(); i++)
            total2+=(1/vCoteArea[i].second);
        for(int i=0; i<vCoteArea.size(); i++)
            vCoteArea[i].second = (1/vCoteArea[i].second) / total2;
    }

    // intégration de power
    double total3=0;
    for(int i=0; i<vCoteArea.size(); i++){
        vCoteArea[i].second = pow(vCoteArea[i].second, power);
        total3+=vCoteArea[i].second;
    }

    sort(vCoteArea.begin(), vCoteArea.end(), maxPairIntDouble);

    int res=0;
    double parcouru=vCoteArea[res].second;

    double randValue = rand()/((double)RAND_MAX+1)*total3;

    while(randValue>parcouru && res!=vCoteArea.size()-1){
        res++;
        parcouru+=vCoteArea[res].second;
    }
    int cote = (vCoteArea[res].first == vCote.size()) ? vPrev[nodeId] : vCote[vCoteArea[res].first];
    long area = (vCoteArea[res].first == vCote.size()) ? currentArea : vArea[vCoteArea[res].first];
    return pair<int, long> (cote, area-currentArea);
}
