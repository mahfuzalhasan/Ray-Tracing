#ifndef INCLUDE_CYLINDER_HPP
#define INCLUDE_CYLINDER_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"
#include "Shape.hpp"
#include "Ray.hpp"
#include<windows.h>

#include "glut.h"
#define PI 3.141592653589793238462643383
#define EPS 1e-10
#define INF 5000

class Cylinder:public Shape{
public:
	float xCenter, yCenter, radius, height,zMin,zMax;
	GLUquadricObj *quadratic;
	void setvalues(float xCenter, float yCenter, float radius, float yMin,float yMax){
		this->xCenter=xCenter;
		this->yCenter=yCenter;
		this->radius=radius;
		this->zMin=yMin;
		this->zMax=yMax;
		this->height=yMax-yMin;
	}
	void setxCenter(float tmp){this->xCenter=tmp;}
	void setyCenter(float tmp){this->yCenter=tmp;}
	void setradius(float tmp){this->radius=tmp;}
	void setzMin(float tmp){this->zMin=tmp;}
	void setzMax(float tmp){this->zMax=tmp;}
	void drawCircle(float x, float y, float z,float radius) {
		int circle_points=100;
		//glClear(GL_COLOR_BUFFER_BIT);
		glPushMatrix();
		//glLoadIdentity();
		glTranslatef(x, y, z);
		double angle = 2*  PI/circle_points ;
		glPolygonMode( GL_FRONT, GL_FILL );
		//glColor3f(1.0, 1.0, 1.0 );
		glBegin(GL_POLYGON);
			double angle1=0.0;
			glVertex2d( radius * cos(0.0) , radius * sin(0.0));
			int i;
			for ( i=0 ; i< circle_points ;i++)
			{

				glVertex2d(radius * cos(angle1), radius *sin(angle1));
				angle1 += angle ;
			}
		glEnd();
		glPopMatrix();
	}

	void drawCylinder(){
		quadratic = gluNewQuadric();
		glPushMatrix();{
				glColor3f(color[0],color[1],color[2]);
				glTranslatef(xCenter,yCenter,zMin);
				//glRotatef(-90,1,0,0);
				//glBegin;
				gluCylinder(quadratic, radius, radius, height, 20, 20);
				//glEnd();
				drawCircle(0,0,0,radius);
				drawCircle(0,0,height,radius);
		}glPopMatrix();
	}
	float findlowestPositive(double t1,double t2,double t3,double t4){
		float mini=INF;
		if(t1>=0 && t1<mini) mini=t1;
		if(t2>=0 && t2<mini) mini=t2;
		if(t3>=0 && t3<mini) mini=t3;
		if(t4>=0 && t4<mini) mini=t4;
		return mini;
	}
	Point getNormal(Point P){
		if(P.getVectZ()==zMin) return Point(0,0,1);
		else if(P.getVectZ()==zMax) return Point(0,0,1);
		else {
			Point C(xCenter,yCenter,P.getVectZ());
			Point n=P-C;
			return n.UnitVect();
		}
	}
	Point findintersection(Ray Rh){
		double a,b,d,B,A,C;
		a=Rh.origin.x-this->xCenter;
		b=Rh.origin.y-this->yCenter;
		A= Rh.direction.x*Rh.direction.x + Rh.direction.y*Rh.direction.y;
		B= 2*(Rh.direction.x*a + Rh.direction.y*b);
		C= (a*a) + (b*b) - (this->radius*this->radius);
		d= (B*B) - (4.0*A*C);
		if(d<0) return Point(INF,INF,INF);
		d=sqrt(d);
		double t,t3=-1,t4=-1,t2,t1;
		t1= (-B + d)/(2*A);
		t2= (-B - d)/(2*A);
		double zcomp2,zcomp1=Rh.origin.z+ t1* Rh.direction.z;
		zcomp2=Rh.origin.z+ t2* Rh.direction.z;
		if(zcomp1<this->zMin || zcomp1>this->zMax ) t1=-1e10;
		if(zcomp2<this->zMin || zcomp2>this->zMax ) t2=-1e10;
		Point n(0,0,1);
		double D=-this->zMax;
		t3= -(D+n.VectDot(Rh.origin))/(n.VectDot(Rh.direction));
		n.setPoint(0,0,-1);
		D=-this->zMin;
		t4= -(D+n.VectDot(Rh.origin))/(n.VectDot(Rh.direction));
		double d1,d2;
		Point P1=Rh.origin + t3 * Rh.direction;
		Point c(this->xCenter,this->yCenter,this->zMax);
		d1= (P1-c).magnitude();
		P1= Rh.origin + t4 * Rh.direction;
		c.setPoint(this->xCenter,this->yCenter,this->zMin);
		d2=(P1-c).magnitude();
		if(d1>(this->radius * this->radius)) t3=-1e10;
		if(d2>this->radius) t4=-1e10;
		//cout<<"cylinder t1:"<<t1<<" t2:"<<t2<<" t3:"<<t3<<" t4:"<<t4<<endl;
		t=this->findlowestPositive(t1,t2,t3,t4);
		if(t>=INF || t<EPS) return Point(INF,INF,INF);
		else return (Rh.origin + t* Rh.direction);

	}
	bool doesIntersect(Point P){
		if(P.x==INF /*&& P.y==INF && P.z==INF*/) return false;
		else return true;

	}
};

#endif
