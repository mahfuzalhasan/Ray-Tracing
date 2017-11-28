#ifndef INCLUDE_RAY_HPP
#define INCLUDE_RAY_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"

#define EPS 1e-10
class Ray{
 public:
	Point origin, direction;
	Ray(){
		origin= Point(0,0,0);
		direction= Point(0,0,0);
	}
	Ray(Point o, Point d){
		this->origin=o;
		this->direction=d.UnitVect();

	}
	void setRay(Point o, Point d){
		this->origin.setPoint(o);
		this->direction.setPoint(d.UnitVect());
	}
	Ray getReflectedRay(Point hitP, Point normal){
		Point r=this->direction - 2*(this->direction.VectDot(normal)) * normal;
		hitP= hitP + 0.01 * r;
		return Ray(hitP,r);
	}
	Ray getRefractedRay(Point hitP, Point normal, float refIdx){
		double Cosi= this->direction.VectDot(normal);
		double Sini= this->direction.VectCross(normal).magnitude2();
		Sini= Sini*Sini;
		Sini=sqrt((1-Sini));
		float r=1/refIdx;
		double scalar= r*Cosi - Sini;
		Point d= (r*this->direction) + (scalar * normal);
		hitP= hitP + 0.01 * d;
		return Ray(hitP,d);

	}
};

#endif