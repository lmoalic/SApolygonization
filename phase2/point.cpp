#include<cmath>
#include "point.h"

point::point()
{
    id=x=y=-1;
}

point::point(int _id, int _x, int _y)
{
    id=_id;
    x=_x;
    y=_y;
}


double point::dist(const point& pt){
    return sqrt(pow(pt.x-x, 2) + pow(pt.y-y, 2));
}


// returns 1, -1, or 0 depending on whether the current point (this) is on the left of, on the right of, or on the oriented straight-line (p0,p1)
int point::left(point& p0, point& p1)
{
    long long det=((long long)p1.x-p0.x)*(y-p0.y) - ((long long)x-p0.x)*(p1.y-p0.y);

    if (det > 0)
        return 1;
    if (det < 0)
        return -1;
    return 0;
}

// return the scalar between pp0 and p1p2
double point::prodScal(const point& p0, const point& p1, const point& p2){
    double ux=p0.x-x;
    double uy=p0.y-y;
    double vx=p2.x-p1.x;
    double vy=p2.y-p1.y;

    return ux*vx + uy*vy;
}



double point::getAngle(point& p1, point& p2){
    double ux=p1.x-x;
    double uy=p1.y-y;
    double vx=p2.x-x;
    double vy=p2.y-y;
    return acos((ux*vx+uy*vy)/(sqrt(ux*ux+uy*uy)*sqrt(vx*vx+vy*vy)));
}


std::ostream & operator<<(std::ostream &out, const point& pt){
    out<<'p'<<pt.id<<"(x="<<pt.x<<", y="<<pt.y<<')';
    return out;
}
