#ifndef INCLUDE_LIGHTSOURCE_HPP
#define INCLUDE_LIGHTSOURCE_HPP

#include <cmath>
#include <fstream>
#include <algorithm>
#include "Point.hpp"
#include "Shape.hpp"
#include "Ray.hpp"
#include<windows.h>

#include "glut.h"

using namespace std;
#define N 10
#define EPS 1e-10
#define INF 5000
class LightSource{
public:
	Point Light[N];
	int sz;
	float Ia,Ip;
	LightSource(){
		sz=0.;
		Ia=1.;
		Ip=1;
	}
	void insertLight(double x,double y,double z){
		Light[sz].setPoint(x,y,z);
		sz++;
	}
	double findLambert(Point P,int i,Point n){
		Point L= (this->Light[i]- P).UnitVect();
		return L.VectDot(n);

	}
	double findSpecular(Point P,int i,Point n,Point v){
		Point l= (this->Light[i]-P).UnitVect();
		Point h= l + (P-v).UnitVect();
		return h.VectDot(n);
	}
	Point findIntensity(Point P,Point n,Point v,float amb,float dif,float spec,float exp){
		int i;
		Point I(0,0,0);
		float lam,pp,p,totalI=0;
		//this->Ip=(1.0/this->sz);
		for(i=0;i<this->sz;i++){
			float tp=findLambert(P,i,n);
			float tmp=max(tp,0.0f);
			lam=dif*tmp;
			p=findSpecular(P,i,n,v);
			tp=pow(p,exp); //(r.v)^k
			tmp=max(tp,0.0f);
			pp=spec * tmp; //phong
			totalI= totalI + Ip*(lam +pp);
		}
		I.x = Ia*amb +  totalI;
		I.y = Ia*amb +  totalI;
		I.z = Ia*amb +  totalI;
		return I;
	}
};
#endif
