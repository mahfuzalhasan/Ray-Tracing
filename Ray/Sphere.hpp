#ifndef INCLUDE_SPHERE_HPP
#define INCLUDE_SPHERE_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"
#include "Shape.hpp"
#include "Ray.hpp"
#include<windows.h>

#include "glut.h"

#define EPS 1e-10
#define INF 5000
class Sphere:public Shape{
public:
	Point center;
	float radius;
	int sz;
	Sphere(){sz=0;}
	void setvalues(Point c, float r){
			center.x=c.x;
			center.y=c.y;
			center.z=c.z;
			radius=r;
	}
	void setcenter(Point p){
		this->center.x=p.x;
		this->center.y=p.y;
		this->center.z=p.z;
	}
	void setradius(float r){
		this->radius=r;
	}
	void drawSphere(){
		glPushMatrix();{
		glColor3f(color[0],color[1],color[2]);
		glTranslatef(center.x,center.y,center.z);
		glutSolidSphere(radius,100,100);
		}glPopMatrix();
	}
	Point getNormal(Point P){
		Point n=P - this->center;
		return n.UnitVect();
	}
	Point findintersection(Ray Rh){
		double d,A=Rh.direction.VectDot(Rh.direction);
		Point tmp=Rh.origin-this->center;
		double B=2*tmp.VectDot(Rh.direction);
		double C=tmp.VectDot(tmp) - (this->radius*this->radius);
		double img= (B*B) -(4.0*A*C);
		if(img<0) return Point(INF,INF,INF);
		d=sqrt(img);
		//cout<<"A:"<<A<<" B:"<<B<<"C:"<<C<<endl;
		double t,t2,t1= (-B + d)/2*A;
		t2= (-B - d)/2*A;
		//cout<<"sphre t1:"<<t1<<" t2:"<<t2<<"img "<<"B:"<<B<<endl;
		if(t1<0 && t2<0)     return Point(INF,INF,INF);
		else if(t1 < 0)  t = t2;
		else if(t2 < 0) t = t1;
		else if(t1 < t2) t = t1;
		else t = t2;
		if(t <= EPS)    return Point(INF,INF,INF);

		return (Rh.origin + Rh.direction*t);
	}
	bool doesIntersect(Point P){
		if(P.x==INF && P.y==INF && P.z==INF) return false;
		else return true;

	}
};

#endif
