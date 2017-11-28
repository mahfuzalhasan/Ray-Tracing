#ifndef INCLUDE_SHAPE_HPP
#define INCLUDE_SHAPE_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"
#include "Ray.hpp"
#define INF 5000
class Shape{
public:
	float ambCoeff,difCoeff,refCoeff,specCoeff,specExp;
	float color[3];
	Point IP;
public:
	void setIP(Point P){
		IP.x=P.x;
		IP.y=P.y;
		IP.z=P.z;
	}
	Point getIP(){return IP;}
	void setCoeff(float amb,float dif, float ref, float spec, float exp){
		ambCoeff=amb;
		difCoeff=dif;
		refCoeff=ref;
		specCoeff=spec;
		specExp=exp;
	}
	void setambCoeff(float amb){ambCoeff=amb;}
	void setdifCoeff(float dif){difCoeff=dif;}
	void setrefCoeff(float ref){refCoeff=ref;}
	void setspecCoeff(float spec){specCoeff=spec;}
	void setspecExp(float exp){specExp=exp;}
	void setColor(float r,float g,float b){
		color[0]=r;
		color[1]=g;
		color[2]=b;
	}
	
};

#endif