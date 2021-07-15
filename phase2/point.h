#ifndef POINT_H
#define POINT_H

#include<ostream>

class point
{
public:
    point();
    point(int _id, int _x, int _y);
    double dist(const point& pt);
    int left(point& p0, point& p1);
    double prodScal(const point& p0, const point& p1, const point& p2);
    double getAngle(point& p1, point& p2);
    friend std::ostream& operator<<(std::ostream &out, const point& pt);

    int id;
    double x, y;
};

#endif // POINT_H
