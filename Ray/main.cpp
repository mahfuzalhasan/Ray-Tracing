
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string>
#include<vector>
#include<utility>
#include<algorithm>
#include<fstream>
#include<iostream>
#include <sstream>
#define NOMINMAX
#include<windows.h>
#include<bits/stdc++.h>

#include "glut.h"
#include "bitmap_image.hpp"



#define BLACK 0, 0, 0
#define N 10
#define VIEW 548.35
#define PI 3.141592653589793238462643383
#define Cos(th) cos(PI/180.0*th)
#define Sin(th) sin(PI/180.0*th)
#define INF 5000

#define N 10
#define EPS 1e-10
#define FL float

using namespace std;



int recursion_level;
int window_size;
int no_of_obj;

int no_of_lights;

//make a global variable -- for tracking the anglular position of camera //in radian
//GLfloat unset[]={0,0,0,1};
int recDepth,pixels;



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



class TriangleList{
public:
	Triangle TR[N];
	int sz;
	TriangleList(){sz=0;}
	void insertTriangle(Triangle t){
		TR[sz].setCoeff(t.ambCoeff,t.difCoeff,t.refCoeff,t.specCoeff,t.specExp);
		TR[sz].setColor(t.color[0],t.color[1],t.color[2]);
		TR[sz].setvalues(t.a,t.b,t.c);
		TR[sz].setrefIdx(t.refIdx);
		sz++;
		//cout<<"Size:"<<sz<<endl;
	}
};




queue<Triangle>sT;


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



queue<Sphere>spheres;

void copyQueue(queue<Sphere>s1,queue<Sphere>s2)
{
    int s = s1.size();
    for(int i=0;i<s1.size();i++)
    {
        Sphere s = s1.front();
        s1.pop();
        s2.push(s);
        s1.push(s);
    }
}

class SphereList{
public:
	Sphere SP[N];
	int sz;
	SphereList(){sz=0;}
	void insertSphere(Sphere s){
		SP[sz].setCoeff(s.ambCoeff,s.difCoeff,s.refCoeff,s.specCoeff,s.specExp);
		SP[sz].setColor(s.color[0],s.color[1],s.color[2]);
		SP[sz].setvalues(s.center,s.radius);
		sz++;
	}
};
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

/*class SphereList{
public:
	Sphere SP[N];
	int sz;
	SphereList(){sz=0;}
	void insertSphere(Sphere s){
		SP[sz].setCoeff(s.ambCoeff,s.difCoeff,s.refCoeff,s.specCoeff,s.specExp);
		SP[sz].setColor(s.color[0],s.color[1],s.color[2]);
		SP[sz].setvalues(s.center,s.radius);
		sz++;
	}
};*/





Point Front,Side,Up,CameraPos,CameraLook;
float mat[16];
double Len(double x,double y,double z){
		return sqrt(x*x+y*y+z*z);
}
void drawGrid(){
	int i;
		glColor3f(0.3, 0.3, 0.3);	//grey
		glBegin(GL_LINES);{
			for(i=-10;i<=10;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -100, 0);
				glVertex3f(i*10,  100, 0);

				//lines parallel to X-axis
				glVertex3f(-100, i*10, 0);
				glVertex3f( 100, i*10, 0);
			}
		}glEnd();

	// draw the two AXES
	glColor3f(1, 1, 1);	//100% white
	glBegin(GL_LINES);{
		//Y axis
		glVertex3f(0, -150, 0);	// intentionally extended to -150 to 150, no big deal
		glVertex3f(0,  150, 0);

		//X axis
		glVertex3f(-150, 0, 0);
		glVertex3f( 150, 0, 0);
	}glEnd();

}
void LoadCamera(){
	glGetFloatv(GL_MODELVIEW_MATRIX, mat);
	for(int i=0;i<15;i++)
    {
        cout<<"Mat"<<i<<":"<<mat[i]<<" ";
    }
	double l=Len(mat[2], mat[6], mat[10]);
	Front = Point(mat[2]/l, mat[6]/l, mat[10]/l);
	l=Len(mat[0], mat[4], mat[8]);
	Side = Point(mat[0]/l, mat[4]/l, mat[8]/l);
	l=Len(mat[1], mat[5], mat[9]);
	Up = Point(mat[1]/l, mat[5]/l, mat[9]/l);
/*	double l=Len(-1, 0, 0);
	Front = Point(-1/l, 0/l, 0/l);
	l=Len(0, -1, 0);
	Side = Point(0/l, -1/l, 0/l);
	l=Len(0, 0, 1);
	Up = Point(0/l, 0/l, 1/l);*/
}
void initCamera(){
	CameraPos.setPoint(-100.0,0,20.0);
	CameraLook.setPoint(0.0,0.0,20.0);
	Up.setPoint(0,0,1.0);
}
void Walking(double f)
{
		LoadCamera();
		CameraPos.setPoint(CameraPos.x+Front.x*f,CameraPos.y+Front.y*f,CameraPos.z+Front.z*f);
		CameraLook.setPoint(CameraLook.x+Front.x*f,CameraLook.y+Front.y*f,CameraLook.z+Front.z*f);
}
void Strafing(double f){
	  LoadCamera();
	  CameraPos.setPoint(CameraPos.x+Side.x*f,CameraPos.y+Side.y*f,CameraPos.z+Side.z*f);
	  CameraLook.setPoint(CameraLook.x+Side.x*f,CameraLook.y+Side.y*f,CameraLook.z+Side.z*f);
}
void Flying(double f){
	LoadCamera();
	CameraPos.setPoint(CameraPos.x+Up.x*f,CameraPos.y+Up.y*f,CameraPos.z+Up.z*f);
	CameraLook.setPoint(CameraLook.x+Up.x*f,CameraLook.y+Up.y*f,CameraLook.z+Up.z*f);
}
void Pitching(double th){
	glPushMatrix();
		{
			glRotatef(th,Side.x,Side.y,Side.z);
			LoadCamera();
		}
	glPopMatrix();
	CameraLook.setPoint(CameraPos.x-Front.x,CameraPos.y-Front.y,CameraPos.z-Front.z);
}
void Yawing(double th){
	glPushMatrix();
		{
			glRotatef(th,Up.x,Up.y,Up.z);
			LoadCamera();
		}
	glPopMatrix();
	CameraLook.setPoint(CameraPos.x-Front.x,CameraPos.y-Front.y,CameraPos.z-Front.z);
}
void Rolling(double th){
	glPushMatrix();
		{
			glRotatef(th,Front.x,Front.y,Front.z);
			LoadCamera();
		}
	glPopMatrix();
	CameraLook.setPoint(CameraPos.x-Front.x,CameraPos.y-Front.y,CameraPos.z-Front.z);
}





class Pyramid:public Shape{
public:
    Point lowest_corner;
    float width,height;

    void setParametres(Point c,float w,float h)
    {
        lowest_corner.x=c.x;
        lowest_corner.y=c.y;
        lowest_corner.z=c.z;
        width=w;
        height=h;
    }
    void drawPyramid()
    {
        Triangle t1;
        Triangle t2;
        Triangle t3;
        Triangle t4;
        float a = lowest_corner.x+30;
        float b = lowest_corner.y-40;
        float c = lowest_corner.z;
/*        t1.seta(Point(a,b,c));
			t1.setb(Point(a,b+width,c));
			t1.setc(Point(a+width/2,b+width/2,c+height));
			t1.setColor(color[0],color[1],color[2]);
			t1.setambCoeff(this->ambCoeff);
			t1.setdifCoeff(this->difCoeff);
			t1.setspecCoeff(this->specCoeff);
			t1.setrefCoeff(this->refCoeff);
            //rt.TL.insertTriangle(t1);
           // sT.push(t1);

            t2.seta(Point(a,b-width,c));
			t2.setb(Point(a+width,b-width,c));
			t2.setc(Point(a+width/2,b+width/2,c+height));
			t2.setColor(color[0],color[1],color[2]);
			t2.setambCoeff(this->ambCoeff);
			t2.setdifCoeff(this->difCoeff);
			t2.setspecCoeff(this->specCoeff);
			t2.setrefCoeff(this->refCoeff);
            //rt.TL.insertTriangle(t2);
            //sT.push(t2);

            t3.seta(Point(a+width,b-width,c));
			t3.setb(Point(a+width,0,c));
			t3.setc(Point(a+width/2,b+width/2,c+height));
			t3.setColor(color[0],color[1],color[2]);
			t3.setambCoeff(this->ambCoeff);
			t3.setdifCoeff(this->difCoeff);
			t3.setspecCoeff(this->specCoeff);
			t3.setrefCoeff(this->refCoeff);
            //rt.TL.insertTriangle(t3);
            //sT.push(t3);

            t4.seta(Point(a+width,0,c));
			t4.setb(Point(a,b,c));
			t4.setc(Point(a+width/2,b+width/2,c+height));
			t4.setColor(color[0],color[1],color[2]);
			t4.setambCoeff(this->ambCoeff);
			t4.setdifCoeff(this->difCoeff);
			t4.setspecCoeff(this->specCoeff);
			t4.setrefCoeff(this->refCoeff);*/
        glColor3f(color[0],color[1],color[2]);

        glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(a,b,c);
				glVertex3f(a,b+width,c);
				glVertex3f(a+width,b+width,c);
				glVertex3f(a+width,b,c);
			}glEnd();

			//float m = (sqrt(pow(30,2)+pow(30,2)))/2;


            //rt.TL.insertTriangle(t4);
            //sT.push(t4);

            glBegin(GL_TRIANGLES);{

                glVertex3f(a,b,c);
				glVertex3f(a,b+width,c);
                glVertex3f(a+width/2,b+width/2,c+height);

				glVertex3f(a,b+width,c);
				glVertex3f(a+width,b+width,c);
                glVertex3f(a+width/2,b+width/2,c+height);

                glVertex3f(a+width,b+width,c);
				glVertex3f(a+width,b,c);
                glVertex3f(a+width/2,b+width/2,c+height);

                glVertex3f(a+width,b,c);
                glVertex3f(a,b,c);
                glVertex3f(a+width/2,b+width/2,c+height);
            }glEnd();
    }

};



class PyramidList{
public:
    Pyramid P[N];
    int sz;
    PyramidList(){sz=0;}
    void insertPyramid(Pyramid s){
        P[sz].setCoeff(s.ambCoeff,s.difCoeff,s.refCoeff,s.specCoeff,s.specExp);
        P[sz].setColor(s.color[0],s.color[1],s.color[2]);
        P[sz].setParametres(s.lowest_corner,s.width,s.height);
        sz++;
    }
};





class RayTracer{
public:
	int filecnt;
	int width,height,recdepth,CO[3];
	double disEyePixel;
	checkerBoard CB;
	SphereList SL;

	TriangleList TL;
	LightSource LS;
	PyramidList PL;

	RayTracer(){
	    filecnt=1;
	    int k = sT.size();
	    for(int i=0;i<k;i++)
        {
            this->TL.insertTriangle(sT.front());
            sT.pop();
        }

	}
	static bool cmp( const pair<int,double> &a, const pair<int, double> &b){
		if(a.second==b.second) return a.first < b.first;
		else return a.second < b.second;
	}
	void setDimension(int w,int h){
		this->width=w;
		this->height=h;
		disEyePixel=this->height/(2*tan(VIEW*PI/360));
	}
	void setrecdepth(int r){this->recdepth=r;}
	void convertColor(float r,float g,float b){
		float f2;
		f2 = max(0.0f, min(1.0f, r));
		this->CO[0] =(int) floor(f2 == 1.0 ? 255 : f2 * 256.0);
		f2 =  max(0.0f, min(1.0f, g));
		this->CO[1] =(int) floor(f2 == 1.0 ? 255 : f2 * 256.0);
		f2 =  max(0.0f, min(1.0f, b));
		this->CO[2] =(int) floor(f2 == 1.0 ? 255 : f2 * 256.0);

	}
	void GenerateTriangle(Point lowest_corner, Pyramid p)
	{
	    Triangle t1,t2,t3,t4;
	    float a = lowest_corner.x+30;
        float b = lowest_corner.y-40;
        float c = lowest_corner.z;
        float ambCoeff = p.ambCoeff;
        float difCoeff = p.difCoeff;
        float specCoeff = p.specCoeff;
        float refCoeff = p.refCoeff;
        float width = p.width;
        float height = p.height;



        //float m = (sqrt(pow(30,2)+pow(30,2)))/2;
			t1.seta(Point(a,b,c));
			t1.setb(Point(a,b+width,c));
			t1.setc(Point(a+width/2,b+width/2,c+height));
			t1.setColor(p.color[0],p.color[1],p.color[2]);
			t1.setambCoeff(p.ambCoeff);
			t1.setdifCoeff(p.difCoeff);
			t1.setspecCoeff(p.specCoeff);
			t1.setrefCoeff(p.refCoeff);
            this->TL.insertTriangle(t1);
            //sT.push(t1);

            t2.seta(Point(a,b+width,c));
			t2.setb(Point(a+width,b+width,c));
			t2.setc(Point(a+width/2,b+width/2,c+height));
			t2.setColor(p.color[0],p.color[1],p.color[2]);
			t2.setambCoeff(p.ambCoeff);
			t2.setdifCoeff(p.difCoeff);
			t2.setspecCoeff(p.specCoeff);
			t2.setrefCoeff(p.refCoeff);
            this->TL.insertTriangle(t2);
            //sT.push(t2);

            t3.seta(Point(a+width,b+width,c));
			t3.setb(Point(a+width,b,c));
			t3.setc(Point(a+width/2,b+width/2,c+height));
			t3.setColor(p.color[0],p.color[1],p.color[2]);
			t3.setambCoeff(p.ambCoeff);
			t3.setdifCoeff(p.difCoeff);
			t3.setspecCoeff(p.specCoeff);
			t3.setrefCoeff(p.refCoeff);
            this->TL.insertTriangle(t3);
            //sT.push(t3);

            t4.seta(Point(a+width,b,c));
			t4.setb(Point(a,b,c));
			t4.setc(Point(a+width/2,b+width/2,c+height));
			t4.setColor(p.color[0],p.color[1],p.color[2]);
			t4.setambCoeff(p.ambCoeff);
			t4.setdifCoeff(p.difCoeff);
			t4.setspecCoeff(p.specCoeff);
			t4.setrefCoeff(p.refCoeff);
            this->TL.insertTriangle(t4);



	}
	Point getReflectedColor(int depth,Point hitPoint,Ray R,Point normal,Point objColor,float refC){
		if(depth==this->recdepth || (objColor.x==0 && objColor.y==0 && objColor.z==0)){
			return objColor;
		}
		else{
			Ray ref=R.getReflectedRay(hitPoint,normal);
			int i,c=-1,minidx,minTyp;
			double mindis=INF;
			Point minHitPoint;
			for(i=0;i<SL.sz;i++){
				Point ip=SL.SP[i].findintersection(ref); //interseting point
				//cout<<"find intersection point of sphere "<<ip.x<<" "<<ip.y<<" "<<ip.z<<endl;
				if(SL.SP[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					//cout<<"sphere hit dis:"<<tmp<<endl;
					if(tmp<mindis){
						mindis=tmp;
						minTyp=1;
						minidx=i;
						minHitPoint.setPoint(ip);
					}
				}
			}
			for(i=0;i<TL.sz;i++){
				Point ip=TL.TR[i].findintersection(ref);
				if(TL.TR[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					if(tmp<mindis){
						mindis=tmp;
						minTyp=2;
						minidx=i;
						minHitPoint.setPoint(ip);
					}
				}
			}

		//checkerboard intersection
			Point ip=CB.findintersection(ref);
			if(CB.doesIntersect(ip)){
				double tmp= (ref.origin-ip).magnitude2();
				if(tmp<mindis){
						mindis=tmp;
						minTyp=4;
						minidx=0;
						minHitPoint.setPoint(ip);
						c=CB.getColor(ip);
				}
			}
			if(mindis<INF) {
				Point n;
				Point colr;
				float tmprefC;
				if(minTyp==1){ //sphere
					n.setPoint(SL.SP[minidx].getNormal(minHitPoint));
					colr.setPoint(SL.SP[minidx].color[0],SL.SP[minidx].color[1],SL.SP[minidx].color[2]);
					tmprefC=SL.SP[minidx].refCoeff;
				}
				else if(minTyp==2){ //triangle
					n.setPoint(TL.TR[minidx].getNormal());
					colr.setPoint(TL.TR[minidx].color[0],TL.TR[minidx].color[1],TL.TR[minidx].color[2]);
					tmprefC=TL.TR[minidx].refCoeff;
				}

				else { //checkerboard
					n.setPoint(CB.getNormal());
					if(c==1) colr.setPoint(CB.color[0],CB.color[1],CB.color[2]);
					else if(c==2) colr.setPoint(CB.color2[0],CB.color2[1],CB.color2[2]);
					tmprefC=CB.refCoeff;
				}
				Point recF=objColor + refC * getReflectedColor(depth+1,minHitPoint,ref,n,colr,tmprefC);
				return recF;
			}
			else{
				Point INFP(INF,INF,INF);
				Point Black(0.0f,0.0f,0.0f);
				Point recF=objColor + refC * getReflectedColor(depth+1,INFP,ref,INFP,Black,0);
				return recF;
			}


		}
	}
	Point getRefractedColor(int depth,Point hitPoint,Ray R,Point normal,Point objColor,float refIdx,float specf){
		if(depth==this->recdepth || (objColor.x==0 && objColor.y==0 && objColor.z==0)){
			return objColor;
		}
		else{
			Ray ref=R.getRefractedRay(hitPoint,normal,refIdx);
			int i,c=-1,minidx,minTyp;
			double mindis=INF;
			Point minHitPoint;
			for(i=0;i<SL.sz;i++){
				Point ip=SL.SP[i].findintersection(ref); //interseting point
				//cout<<"find intersection point of sphere "<<ip.x<<" "<<ip.y<<" "<<ip.z<<endl;
				if(SL.SP[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					//cout<<"sphere hit dis:"<<tmp<<endl;
					if(tmp<mindis){
						mindis=tmp;
						minTyp=1;
						minidx=i;
						minHitPoint.setPoint(ip);
					}
				}
			}
			for(i=0;i<TL.sz;i++){
				Point ip=TL.TR[i].findintersection(ref);
				if(TL.TR[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					if(tmp<mindis){
						mindis=tmp;
						minTyp=2;
						minidx=i;
						minHitPoint.setPoint(ip);
					}
				}
			}

		//checkerboard intersection
			Point ip=CB.findintersection(ref);
			if(CB.doesIntersect(ip)){
				double tmp= (ref.origin-ip).magnitude2();
				if(tmp<mindis){
						mindis=tmp;
						minTyp=3;
						minidx=0;
						minHitPoint.setPoint(ip);
						c=CB.getColor(ip);
				}
			}
			if(mindis<INF) {
				Point n;
				Point colr;
				float tmprefC;
				if(minTyp==1){ //sphere
					n.setPoint(SL.SP[minidx].getNormal(minHitPoint));
					colr.setPoint(SL.SP[minidx].color[0],SL.SP[minidx].color[1],SL.SP[minidx].color[2]);
					tmprefC=SL.SP[minidx].specCoeff;
				}
				else if(minTyp==2){ //triangle
					n.setPoint(TL.TR[minidx].getNormal());
					colr.setPoint(TL.TR[minidx].color[0],TL.TR[minidx].color[1],TL.TR[minidx].color[2]);
					tmprefC=TL.TR[minidx].specCoeff;
				}

				else { //checkerboard
					n.setPoint(CB.getNormal());
					if(c==1) colr.setPoint(CB.color[0],CB.color[1],CB.color[2]);
					else if(c==2) colr.setPoint(CB.color2[0],CB.color2[1],CB.color2[2]);
					tmprefC=CB.specCoeff;
				}
				Point recF=objColor +specf * getRefractedColor(depth+1,minHitPoint,ref,n,colr,refIdx,tmprefC);
				return recF;
			}
			else{
				Point INFP(INF,INF,INF);
				Point Black(0.0f,0.0f,0.0f);
				Point recF=objColor + specf * getRefractedColor(depth+1,INFP,ref,INFP,Black,0,0);
				return recF;
			}


		}
	}

	Point getFinalColor(Point hitPoint,Ray R,Point normal,Point objColor,float refC,float refIdx,float specC){
		Point Cr= this->getReflectedColor(1,hitPoint,R,normal,objColor,refC);
		if(refIdx>0){
			double Exp,Cosi,Rs=0.3,R0=(1-refIdx)/(1+refIdx);
			R0=R0*R0;
			Cosi= R.direction.VectDot(normal);
			Cosi= 1-Cosi;
			Exp= pow(Cosi,5);
			Rs= R0+(1-R0) * Exp;
			Point Cf=this->getRefractedColor(1,hitPoint,R,normal,objColor,refIdx,specC);
			Point Cfinal= Rs*Cr + (1-Rs)*Cf;
			return Cfinal;
		}
		return Cr;
	}
	void findPixelColor(Ray R, int depth, Point eye ){
		int i,cnt=0,c=-1;
		vector< pair<int, double> > Obj;
		for(i=0;i<SL.sz;i++){
			Point ip=SL.SP[i].findintersection(R); //interseting point
			//cout<<"find intersection point of sphere "<<ip.x<<" "<<ip.y<<" "<<ip.z<<endl;
			if(SL.SP[i].doesIntersect(ip)){
				SL.SP[i].setIP(ip);
				double tmp= (R.origin-ip).magnitude2();
				//cout<<"sphere hit dis:"<<tmp<<endl;
				Obj.push_back(make_pair(i,tmp));
			}
		}
		//cout<<"sphere hit:"<<hitcnt<<endl;
		cnt+=SL.sz;
		for(i=0;i<TL.sz;i++){
			Point ip=TL.TR[i].findintersection(R);
			if(TL.TR[i].doesIntersect(ip)){
				TL.TR[i].setIP(ip);
				double tmp= (R.origin-ip).magnitude2();
				Obj.push_back(make_pair(i+cnt,tmp));
			}
		}
		cnt+=TL.sz;

		//checkerboard intersection
		Point ip=CB.findintersection(R);

		if(CB.doesIntersect(ip)){
			CB.setIP(ip);
			double tmp= (R.origin-ip).magnitude2();
			Obj.push_back(make_pair(cnt,tmp));
			c=CB.getColor(ip);
		}
		Point I;
		if(Obj.size()!=0) {
			sort(Obj.begin(),Obj.end(),cmp);
			if(Obj[0].first<SL.sz){
				int idx=Obj[0].first;
				//cout<<"pixel sphere hit dis:"<<Obj[0].first<<", "<<Obj[0].second<<endl;
				Point ip=SL.SP[idx].getIP();
				Point normal=SL.SP[idx].getNormal(ip);
				Point objColor(SL.SP[idx].color[0],SL.SP[idx].color[1],SL.SP[idx].color[2]);
				Point tmp=LS.findIntensity(ip,normal,eye,SL.SP[idx].ambCoeff,
					SL.SP[idx].difCoeff,SL.SP[idx].specCoeff,SL.SP[idx].specExp/*,objColor*/);
				I.setPoint(tmp);
				Point Cf=getFinalColor(ip,R,normal,objColor,SL.SP[idx].refCoeff,0,SL.SP[idx].specCoeff);
				convertColor(I.x*Cf.x,I.y*Cf.y,I.z*Cf.z);
				//cout<<"CO:("<<SL.SP[idx].color[0]<<","<<SL.SP[idx].color[1]<<","<<SL.SP[idx].color[2]<<")"<<endl;
			}
			else if(Obj[0].first>=SL.sz && Obj[0].first<(SL.sz+TL.sz)){
				int idx=Obj[0].first-SL.sz;
				Point ip=TL.TR[idx].getIP();
				Point normal=TL.TR[idx].getNormal();
				Point objColor(TL.TR[idx].color[0],TL.TR[idx].color[1],TL.TR[idx].color[2]);
				Point tmp=LS.findIntensity(ip,normal,eye,TL.TR[idx].ambCoeff,
					TL.TR[idx].difCoeff,TL.TR[idx].specCoeff,TL.TR[idx].specExp/*,objColor*/);
				I.setPoint(tmp);
				Point Cf=getFinalColor(ip,R,normal,objColor,TL.TR[idx].refCoeff,TL.TR[idx].refIdx,TL.TR[idx].specCoeff);
				convertColor(I.x*Cf.x,I.y*Cf.y,I.z*Cf.z);
				//cout<<"pixel Triangle hit dis:"<<Obj[0].first<<", "<<Obj[0].second<<endl;
				//cout<<"CO:("<<SL.SP[idx].color[0]<<","<<SL.SP[idx].color[1]<<","<<SL.SP[idx].color[2]<<")"<<endl;
			}

			else if(Obj[0].first>=(SL.sz+TL.sz)){
				int idx=Obj[0].first-(SL.sz+TL.sz);
				Point ip=CB.getIP();
				Point normal=CB.getNormal();
				Point objColor;
				if(c==1) objColor.setPoint(CB.color[0],CB.color[1],CB.color[2]);
				else if(c==2) objColor.setPoint(CB.color2[0],CB.color2[1],CB.color2[2]);
				Point tmp=LS.findIntensity(ip,normal,eye,CB.ambCoeff,
					CB.difCoeff,CB.specCoeff,CB.specExp/*,objColor*/);
				I.setPoint(tmp);
				Point Cf=getFinalColor(ip,R,normal,objColor,CB.refCoeff,0,CB.specCoeff);
				convertColor(I.x*Cf.x,I.y*Cf.y,I.z*Cf.z);
				//cout<<"CO:("<<CO[0]<<","<<CO[1]<<","<<CO[2]<<")"<<endl;
			}
		}
		else{
			convertColor(0,0,0); //background color BLACK
		}
		//return Point(255,255,255);
	}
	void drawBMPImage(){
		bitmap_image image(this->width,this->height);
		cout<<"starting ouput.bmp"<<endl;
		image.clear();
		int i,j;
		Ray R;
		LoadCamera();
		//Point leftpixel=(Front * this->disEyePixel) + ((-Side) * ((this->width-1) / 2)) +  (Up * ((this->height-1)/2));
		cout<<"camaraPosition:("<<CameraPos.x<<", "<<CameraPos.y<<", "<<CameraPos.z<<")"<<endl;
		//cout<<"camaraFront:("<<Front.x<<", "<<Front.y<<", "<<Front.z<<")"<<endl;
		Point dir(-Front.x,-Front.y,-Front.z);
		R.setRay(CameraPos,dir);
		cout<<"direction:("<<R.direction.x<<", "<<R.direction.y<<", "<<R.direction.z<<")"<<endl;
		Point midPoint= R.origin + VIEW*R.direction;
		cout<<"midpoint:("<<midPoint.x<<", "<<midPoint.y<<", "<<midPoint.z<<")"<<endl;
		//cout<<"Side:("<<Side.x<<", "<<Side.y<<", "<<Side.z<<")"<<endl;
		//cout<<"Up:("<<Up.x<<", "<<Up.y<<", "<<Up.z<<")"<<endl;
		Point leftpixel=midPoint + (-Side) * (this->width/2) + (Up * (this->height/2));
		cout<<"leftpixel:("<<leftpixel.x<<", "<<leftpixel.y<<", "<<leftpixel.z<<")"<<endl;
		for(i=0;i<this->height;i++){
			for(j=0;j<this->width;j++){
				Point P = leftpixel - (Up * i) + (Side * j);
				P= P + 0.5*(-Up) + 0.5 * Side;
				//cout<<"pixel dimension:("<<P.x<<","<<P.y<<","<<P.z<<")"<<endl;
				P.setPoint(P-CameraPos);
				R.setRay(CameraPos,P/*.UnitVect()*/);
				this->findPixelColor(R,this->recdepth,CameraPos);
				image.set_pixel(j,i,this->CO[0],this->CO[1],this->CO[2]);
				//cout<<"Ray P:"<<P.x<<" "<<P.y<<" "<<P.z<<endl;

			}
			//cout<<"Ray origin:"<<R.origin.x<<" "<<R.origin.y<<" "<<R.origin.z<<endl;
			//cout<<"Ray direction:"<<R.direction.x<<" "<<R.direction.y<<" "<<R.direction.z<<endl;

		}
		ostringstream ss;
		ss << this->filecnt;
		string fileName="output" + ss.str() + ".bmp";
		image.save_image(fileName);
		this->filecnt++;
		image.clear();
		cout<<fileName<<" file created."<<endl;
	}

};

RayTracer rt;



void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(BLACK, 0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//gluLookAt(r*cos(cameraAngle), r*sin(cameraAngle), cameraZ,		0,0,5,		0,0,1);//150 is radius
	//cout<<"Camera Pos:"<<CameraPos.x<<" "<<CameraPos.y<<" "<<CameraPos.z<<endl;
	gluLookAt(CameraPos.x, CameraPos.y, CameraPos.z,  CameraLook.x,CameraLook.y,CameraLook.z, Up.x,Up.y,Up.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);
	//drawGrid();

	/*checkerBoard board;
	board.setCoeff(0.4,0.1,0.4,0.1,0.0);
	board.setColor(1,1,1);
	board.setColor2(0,0,0);*/
	rt.CB.drawCheckerBoard();
	/*Cylinder C;
	C.setvalues(0,0,5,30,50);
	C.setColor(1.0,1.0,0.0);
	C.setCoeff(0.4,0.2,0.2,0.2,5.0);*/
	int i;
	/*for(i=0;i<rt.CL.sz;i++){
		rt.CL.CY[i].drawCylinder();
	}*/
	int size = spheres.size();
	for(i=0;i<size;i++){
            Sphere s = spheres.front();
        spheres.pop();
		s.drawSphere();
		spheres.push(s);
	}
	for(i=0;i<rt.PL.sz;i++){
		rt.PL.P[i].drawPyramid();
	}
	/*for(i=0;i<rt.TL.sz;i++){
		rt.TL.TR[i].drawTriangle();
	}*/
	glutSwapBuffers();
}

void animate(){
	//codes for any changes in Camera
	//codes for any changes in Models


	//MISSING SOMETHING? -- YES: add the following
	glutPostRedisplay();	//this will call the display AGAIN
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){
    case '0':
        rt.drawBMPImage();

		case '1':
			Yawing(-2); //left
			break;
		case '2':	//right
			Yawing(2);
			break;
		case '3':	//up
			Pitching(2);
			break;
		case '4':	//down
			Pitching(-2);
			break;
	    case '5':	//twist left
			Rolling(2);
			break;
		case '6':	//twist right
			Rolling(-2);
			break;
	}
}

void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//back
			Walking(7);
			break;
		case GLUT_KEY_UP:
			Walking(-7);
			break; //forward
		case GLUT_KEY_RIGHT:
			Strafing(7);
			break; //right
		case GLUT_KEY_LEFT:
			Strafing(-7);
			break; //left
		case GLUT_KEY_PAGE_UP:
			Flying(7);
			break; //up
		case GLUT_KEY_PAGE_DOWN:
			Flying(-7);
			break; //down
		default:
			break;
	}
}




void input(){
	ifstream in("description1.txt");
	in>>recursion_level;
	rt.setrecdepth(recursion_level);
	in>>window_size;
	rt.setDimension(window_size,window_size);
	in>>no_of_obj;
	string objType;
	float r,g,b,radius,ambient,diffuse,specular,reflection,shineness,w,h;
	Point P;
	for (int i = 0; i < no_of_obj; i++)
	{
        in>>objType;
        cout<<objType;
        if(objType=="sphere")
        {
            cout<<"ddddddd"<<endl;
            Sphere s;
            in>>P.x>>P.y>>P.z;
            s.setcenter(P);
            in>>radius;
            s.setradius(radius);
            in>>r>>g>>b;
            s.setColor(r,g,b);
            in>>ambient>>diffuse>>specular>>reflection;
            s.setambCoeff(ambient);
            s.setdifCoeff(diffuse);
            s.setspecCoeff(specular);
            s.setrefCoeff(reflection);
            in>>shineness;
            s.setspecExp(shineness);
            rt.SL.insertSphere(s);
            spheres.push(s);

        }
        else if(objType=="pyramid")
        {
            cout<<"ppppp"<<endl;
            Pyramid p;
            in>>P.x>>P.y>>P.z;
            in>>w>>h;
            p.setParametres(P,w,h);
            in>>r>>g>>b;
            p.setColor(r,g,b);
            in>>ambient>>diffuse>>specular>>reflection;
            p.setambCoeff(ambient);
            p.setdifCoeff(diffuse);
            p.setspecCoeff(specular);
            p.setrefCoeff(reflection);
            in>>shineness;
            p.setspecExp(shineness);
            rt.PL.insertPyramid(p);
            rt.GenerateTriangle(P,p);


        }
        rt.CB.setColor(1.0,1.0,1.0);
        rt.CB.setColor2(0,0,0);
        rt.CB.setambCoeff(0.4);
        rt.CB.setdifCoeff(0.1);
        rt.CB.setspecCoeff(0.4);
        rt.CB.setrefCoeff(0.1);

		/*in>>obj[i].type;
		in>>obj[i].pos_x>>obj[i].pos_y>>obj[i].pos_z;
		in>>obj[i].rad;
		in>>obj[i].col_r>>obj[i].col_g>>obj[i].col_b;
		in>>obj[i].coef_a>>obj[i].coef_d>>obj[i].coef_s>>obj[i].coef_r;
		in>>obj[i].shine;*/
	}
	in>>no_of_lights;
	for (int i = 0; i < no_of_lights; ++i)
	{
		in>>P.x>>P.y>>P.z;
		rt.LS.insertLight(P.x,P.y,P.z);
	}
}



void initLighting(){
	//diffusecolor=1.0;
	//diffuseangle=0;
	glShadeModel(GL_SMOOTH);
	GLfloat amb0[] = {1, 1, 1, 1}; // define some colors
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb0);
	GLfloat diffusePoint[] = {1.0, 1.0, 1.0, 1.0}; //Color (0.5, 0.5, 0.5)
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffusePoint); //

    glEnable(GL_NORMALIZE); //Automatically normalize normals needed by the illumination model
    glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	//setSpotLight();
}

void init(){
	//codes for initialization
	initCamera();
	//clear the screen
	glClearColor(BLACK, 0);

	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();
	gluPerspective(70,	1,	0.1,	10000.0);
	//initLighting();
	input();
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	//glutInitWindowSize(window_size, window_size);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();


	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	//ADD keyboard listeners:
	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
