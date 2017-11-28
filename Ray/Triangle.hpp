#ifndef INCLUDE_TRIANGLE_HPP
#define INCLUDE_TRIANGLE_HPP

#include <cmath>
#include <fstream>
#include "Point.hpp"
#include "Shape.hpp"
#include "Ray.hpp"
#include<windows.h>

#include "glut.h"
#define EPS 1e-10
#define INF 5000
class Triangle:public Shape{
public:
	//float a[3],b[3],c[3];
	Point a,b,c;
	float refIdx;
	void setvalues(Point a,Point b,Point c){
		this->a.setPoint(a);
		this->b.setPoint(b);
		this->c.setPoint(c);
	}
	void seta(Point p){
		this->a.setPoint(p);
	}
	void setb(Point p){
		this->b.setPoint(p);
	}
	void setc(Point p){
		this->c.setPoint(p);
	}
	void setrefIdx(float ref){
		this->refIdx=ref;
	}
	Point getNormal(){
		Point u=this->b-this->a;
		Point v=this->c-this->a;
		Point n= u*v;
		return n.UnitVect();
	}
	void drawTriangle(){
		glColor3f(color[0],color[1],color[2]);
		glBegin(GL_POLYGON);
			glVertex3f(a.x,a.y,a.z);
			glVertex3f(b.x,b.y,b.z);
			glVertex3f(c.x,c.y,c.z);
		glEnd();
	}
	Point findintersection(Ray Rh){
		Point    u, v, n;              // triangle vectors
		Point    dir, w0, w;           // ray vectors
		float     r, c, d;              // params to calc ray-plane intersect

		// get triangle edge vectors and plane normal
		u = this->b - this->a;
		v = this->c - this->a;
		n = u * v;              // cross product
		if (n == Point(0,0,0))             // triangle is degenerate
			return Point(INF,INF,INF);                  // do not deal with this case

		//dir = R.P1 - R.P0;              // ray direction vector
		w0 = Rh.origin - this->a;
		c = (float) - n.VectDot(w0);
		d = (float) n.VectDot(Rh.direction);
		if (fabs(d) < EPS) {     // ray is  parallel to triangle plane
			if (c == 0) return Point(INF,INF,INF);      // ray lies in triangle plane
			else return Point(INF,INF,INF);              // ray disjoint from plane
		}

		// get intersect point of ray with triangle plane
		r = c / d;
		if (r < 0.0)                    // ray goes away from triangle
			return Point(INF,INF,INF);                   // => no intersect
		// for a segment, also test if (r > 1.0) => no intersect

		Point I = Rh.origin + r * Rh.direction;            // intersect point of ray and plane

		// is I inside T?
		float    uu, uv, vv, wu, wv, D;
		uu = u.VectDot(u);
		uv = u.VectDot(v);
		vv = v.VectDot(v);
		w = I - this->a;
		wu = w.VectDot(u);
		wv = w.VectDot(v);
		D = uv * uv - uu * vv;

		// get and test parametric coords
		float s, t;
		s = (uv * wv - vv * wu) / D;
		if (s < 0.0 || s > 1.0)         // I is outside T
			return Point(INF,INF,INF);
		t = (uv * wu - uu * wv) / D;
		if (t < 0.0 || (s + t) > 1.0)  // I is outside T
			return Point(INF,INF,INF);

		return I;
	}

	bool doesIntersect(Point P){
		if(P.x==INF) return false;
		else return true;

	}


};

#endif
