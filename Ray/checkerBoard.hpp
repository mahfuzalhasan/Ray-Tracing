#ifndef INCLUDE_CHECKERBOARD_HPP
#define INCLUDE_CHECKERBOARD_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"
#include "Shape.hpp"
#include "Ray.hpp"
#include<windows.h>

#include "glut.h"
#define EPS 1e-10
#define INF 5000

class checkerBoard:public Shape{
public:
	float noX,noY,blocksize;
	float color2[3];
	checkerBoard(){
		noX=1000;
		noY=1000;
		blocksize=8;
	}
	void setColor2(float r, float g, float b){
		color2[0]=r;
		color2[1]=g;
		color2[2]=b;
	}
	void drawCheckerBoard(){
		int i,j,SizeX=this->blocksize,SizeY=this->blocksize;
		float ref=-500;
		glBegin(GL_QUADS);
		for(i=0;i<noX;i++){
			for(j=0;j<noY;j++){
				if((i+j)%2==0) {
					glColor3f( color[0], color[1], color[2]);
				}
				else {
					glColor3f( color2[0], color2[1], color2[2]);
				}
				glVertex2f(ref+i*SizeX,ref+j*SizeY);
				glVertex2f(ref+(i+1)*SizeX,ref+j*SizeY);
				glVertex2f(ref+(i+1)*SizeX,ref+(j+1)*SizeY);
				glVertex2f(ref+i*SizeX,ref+(j+1)*SizeY);

			}
		}
		glEnd();
	}
	Point getNormal(){
		return Point(0,0,1);
	}
	int getColor(Point P){
		int i,j;
		Point X(1,0,0);
		Point Y(0,1,0);
		i=(int) (P.projection(X)/this->blocksize);
		j=(int) (P.projection(Y)/this->blocksize);
		if(((i+j)&1)==0) return 1;
		else return 2;
	}
	Point findintersection(Ray Rh){
		float D=0,t;
		Point n=this->getNormal();
		t=-(D+ n.VectDot(Rh.origin))/(n.VectDot(Rh.direction));
		if(t<0) return Point(INF,INF,INF);
		else return Rh.origin + t*Rh.direction;

	}
	bool doesIntersect(Point P){
		if(P.x==INF && P.y==INF && P.z==INF) return false;
		else return true;
	}
};

#endif
