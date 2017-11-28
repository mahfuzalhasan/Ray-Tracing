#ifndef INCLUDE_POINT_HPP
#define INCLUDE_POINT_HPP

#include <cmath>
#include <fstream>
#include <algorithm>

using namespace std;
#define EPS 1e-10
#define FL float
class Point{
 public:	
	/*double*/ float x,y,z;
	Point(){x=0;y=0;z=0;}
	Point(double tx,double ty,double tz){
		this->x=tx;
		this->y=ty;
		this->z=tz;
	}
	void setPoint(double tx,double ty,double tz){
		this->x=tx;
		this->y=ty;
		this->z=tz;
	}
	void setPoint(Point p){
		this->x=p.x;
		this->y=p.y;
		this->z=p.z;
	}
	FL getVectX(){return x;}
	FL getVectY(){return y;}
	FL getVectZ(){return z;}
	Point VectAdd(Point v){
		return Point(getVectX()+v.getVectX(),getVectY()+v.getVectY(), getVectZ()+v.getVectZ());
	}
	double VectDot(Point v){
		return (getVectX()*v.getVectX()+getVectY()*v.getVectY()+getVectZ()*v.getVectZ());
	}
	Point VectCross(Point v){
		double a = getVectY()*v.getVectZ() - getVectZ()*v.getVectY();
		double b = getVectZ()*v.getVectX() - getVectX()*v.getVectZ();
		double c= getVectX()*v.getVectY() - getVectY()*v.getVectX();
		return Point(a, b, c);
	}
	Point UnitVect(){
		double a = getVectX() * getVectX();
		double b = getVectY() * getVectY();
		double c=  getVectZ() * getVectZ();
		double res= sqrt((a+b+c));
		double X,Y,Z;
		X=getVectX()/res;
		Y=getVectY()/res;
		Z=getVectZ()/res;
		if(X < EPS && X > -EPS) X = 0.;
		if(Y < EPS && Y > -EPS) Y = 0.;
		if(Z < EPS && Z > -EPS) Z = 0.;
		return Point(X,Y,Z);
	}
	
	double magnitude(){
		return (x*x+y*y+z*z);
	}
	double magnitude2(){
		return sqrt(this->magnitude());
	}
	double projection(Point b){
		return (this->VectDot(b)/b.magnitude2());
	}
	friend Point operator+(Point a,Point b){
		return Point(a.x+b.x, a.y+b.y, a.z+b.z);	
	}

	Point operator- (){
		return Point(-this->x, -this->y, -this->z);	
	}
	bool operator==(Point a){
		return (this->x==a.x && this->y==a.y && this->z==a.z);
	}
	friend Point operator-(Point a,Point b){
		return Point(a.x-b.x, a.y-b.y, a.z-b.z);	
	}

	friend Point operator*(Point a,Point b){
		return a.VectCross(b);	
	}
	
	friend Point operator*(Point a, double b){		// A*2.
		return Point(a.x*b, a.y*b, a.z*b);	
	}
	friend Point operator*(double b,Point a){		// A*2.
		return Point(a.x*b, a.y*b, a.z*b);	
	}
	friend Point operator/(Point a, double b){		// A*2.
		return Point(a.x/b, a.y/b, a.z/b);	
	}
	
};


#endif