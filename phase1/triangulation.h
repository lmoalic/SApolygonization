#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <string>
#include <vector>


using namespace std;


extern double PARAM_SA_C_0;
extern double PARAM_SA_C_CS;
extern int PARAM_SA_NB_CS;
extern double PARAM_SA_C_EF;
extern int PARAM_SA_NB_IT;
extern double PARAM_SA_C_IT;
extern int PARAM_SA_NB_TR;
extern double PARAM_SA_C_TR;
extern double PARAM_SA_C_SAR;
extern int PARAM_SA_NB_SC;

class Triangle;

class Point
{
	friend class Triangulation;
	friend class Triangle;
public:
	// Constructors
	Point() {}
	Point(long long abs, long long ord, int ind=0):x(abs),y(ord),num(ind),t{nullptr},e{-1} {}

    /* FAIT*/
    // Returns 1, -1, or 0 depending on whether the current point is on the left of, on the right of, or on the oriented straight-line (ab)
	int onLeft(Point *a, Point *b) const;


public:
	long long x,y;
	int num;	    // the index of the point in LP
	Triangle *t;    // the triangle on the left side of pp', where p is the current point and p' its neighbor in counterclockwise direction on the polygon
	int e;          // the N° of the edge pp' in t
};


class Triangle
{
	friend class Triangulation;
public:
	// Constructor
	Triangle(Point *s0, Point *s1, Point *s2, int col=-1, int ind=0):vertex(3),neighbor(3,nullptr),num(ind),color(col)
	{
		vertex[0]=s0;vertex[1]=s1;vertex[2]=s2;
	}

    // Returns the N° of the edge that the current triangle shares with triangle a
    // or -1 if the two triangles share no edge
	int edge(Triangle *a) const;

	// Returns twice the area of the triangle
    long long area2() const
    {
        return    (vertex[1]->x-vertex[0]->x)*(vertex[2]->y-vertex[0]->y)
                - (vertex[1]->y-vertex[0]->y)*(vertex[2]->x-vertex[0]->x);
    }


private:
	vector<Point*> vertex;          // the three vertices of the triangle in counterclockwise order
	vector<Triangle*> neighbor;     // the three neighbors of the triangle: neighbor[i] is opposite to vertex[i]
	int num;	// the index of the triangle in LT
	int color;  // the color of the triangle: 1 black, 0 white, -1 unbounded
};


class Triangulation
{
public:

    // Default constructor
    Triangulation():area2{0}
    {}

    // Copy constructor
    Triangulation(const Triangulation &T);

    // Overwrites the current triangulation by T.
    // The triangulations must have the same numbers of vertices and triangles.
    void overwriteBy(const Triangulation &T);

    // Destructor
    ~Triangulation();

    // Deletes the vertices of the triangulation
    void deleteVertices();

    // Deletes the triangles of the triangulation
    void deleteTriangles();

    // Returns the area of the current polygonization
    long long getArea()
    {
        return area2/2;
    }

    // Returns twice the standard deviation of the area of the colored triangles of the triangulation
    double areaStdDeviation2() const;

    // Stores in the triangulation twice the area of the Polygon
    void initArea(const vector<int>& Polygon);

    // Adds a new point (x,y) to the list of points of the triangulation.
    // Returns a pointer on this point.
	Point* addPoint(long long x, long long y);

    // Adds a new triangle to the list of triangles of the triangulation.
    // The points s0, s1, and s2 are supposed to already belong to the list of points
    // of the triangulation and become respectively the vertices N°0, 1, and 2
    // of the new triangle.
    // Returns a pointer to the new triangle.
	Triangle* addTriangle(Point *s0, Point *s1, Point *s2, int color=-1);

	// Pastes edge N°i of triangle a on edge N°j of triangle b.
	// a and b are supposed to be already present in the triangle list
	// and their edges i and j are supposed to be "pastable"
	void pasteTriangles(Triangle *a, int i, Triangle *b, int j);

    // Initializes the vector of points from a file with the format given for the competition.
    // The vector is sorted in lexicographic order with respect to (x,y).
    // The point at index 0 is the point at infinity.
    // Returns the number of (finite) points, or -1 if the file could not be opened.
    long long initLP(const std::string &fileName);

    // Saves the ordering of the vertices of the triangulated polygon in a file
    void savePolygon(const string &fichier);

    // Stores in Polygon the ordering of the vertices of the triangulated polygon
    void extractPolygon(std::vector<int>& Polygon);

    // Constructs in Polygon an initial polygonization of the point set LP
    void initialPolygonization(std::vector<int> &Polygon) const;



    // Flips the edge N° c of triangle t in the triangulation. The flip is supposed to be possible.
    void flip(Triangle *t, int c);

    // Swaps the colors of the triangles t1 and t2 which are incident in the vertex s.
    void triangleColorSwap(Triangle *t1, Triangle *t2, Point *s);

    // Lists the color swaps that are possible at vertex s. A color swap is possible only if there exists either exactly
    // one black or exactly one white triangle with vertex s. t1 stores this triangle (if it exists).
    // tri stores the triangles that can be swapped with t1 (if some exist).
    void possibleColorSwaps(Point *s, Triangle *&t1, vector<Triangle *> &tri);

    // Performs a swap color sequence in the given direction: 1 left to right, -1 right to left.
    // Returns true if the best polygonization TBest has been improved.
    // cmpOp is the operator < for minimization and > for maximization.
    template <typename C>
    bool colorSwapSequence(double temp, Triangulation &TBest, int direction, C cmpOp);

    // Performs a sequence of edge flips in the triangulation
    void edgeFlipSequence();

    // Optimizes the area of the current triangulation.
    // The area is minimized if minimization=true and maximized otherwise.
    void areaOptimization(bool minimization);



///////////////////////////////////////////////////////////////////////////////////////////////////
////          Sweep algorithm to construct a colored triangulation of a given polygon          ////
///////////////////////////////////////////////////////////////////////////////////////////////////

    void triangulatePolygon(const std::vector<int> &Polygon);

    void allVisibleFront(Triangle* frontCurrent, Point* rightmostVertex, Point* currentVertex, int order, Triangle*& firstFront, Triangle*& frontNew);

    void allFront(Triangle* frontCurrent, Point* currentVertex, int order, Triangle*& frontLast);

    void beginApexEvent(Triangle* frontSucc, Triangle* frontPred, Point* rightmostVertex, Point* currentVertex, Triangle*& tUp, Triangle*& tDown);

    void endApexEvent(Triangle* frontUpperSucc, Triangle* frontUpperPred, Point* rightmostVertexUpper,
        Triangle* frontLowerSucc, Triangle* frontLowerPred, Point* rightmostVertexLower,
        Triangle* frontMiddleSucc, Triangle* frontMiddlePred, Point* rightmostVertexMiddle,
        Point* currentVertex, Triangle*& tUp, Triangle*& tDown, bool previousDown);

    void intermediateVertexEvent(Triangle* frontUpperSucc, Triangle* frontUpperPred, Point* rightmostVertexUpper,
        Triangle* frontLowerSucc, Triangle* frontLowerPred, Point* rightmostVertexLower,
        Point* currentVertex, Triangle*& tUp, Triangle*& tDown, bool previousBefore);



private:
	vector<Point*> LP;      // list of vertices of the colored triangulation
	vector<Triangle*> LT;   // list of triangles of the colored triangulation

    long long area2;        // twice the area of the triangulated polygon
};


// Data structure and functions for the decomposition of the sweep-line

class SweeplineSegment              // segments that partition the sweep-line
{
public:

    SweeplineSegment(Point* f, Point* l, int c=0):first{f},last{l},rightmost{f},pred{nullptr},succ{nullptr},color{c} {}

    Point *first,*last;     // the first and last endpoints in lexicographic order with respect to (x,y)
    Point *rightmost;       // the rightmost vertex of the sub-front just above the segment
    Triangle *pred, *succ;  // the sub-front triangles out of rightmost (if exist) in CW and CCW order
    int color;              // 0 if the part of the sweep-line just above the segment is outside of the polygon, 1 otherwise
};

// Returns 1, -1, or 0 depending on whether p is on the left of, on the right of, or on the oriented straight-line (s.first,s.last)
int left(const SweeplineSegment &s, const Point &p);

// Returns true if s1 < s2 on the sweepline
bool lessSweepline(const SweeplineSegment &s1, const SweeplineSegment &s2);


#endif // TRIANGULATION_H
