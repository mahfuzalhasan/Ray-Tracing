
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

#include<GL/glut.h>
#include "bitmap_image.hpp"
#include "Point.hpp"
#include "Shape.hpp"
#include "checkerBoard.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Triangle.hpp"
#include "Ray.hpp"
#include "LightSource.hpp"
#define BLACK 0, 0, 0
#define N 10
#define VIEW 548.4088346
#define PI 3.141592653589793238462643383
#define Cos(th) cos(PI/180.0*th)
#define Sin(th) sin(PI/180.0*th)
#define INF 5000
using namespace std;

//make a global variable -- for tracking the anglular position of camera //in radian
//GLfloat unset[]={0,0,0,1};
int recDepth,pixels;
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
class CylinderList{
public:
	Cylinder CY[N];
	int sz;
	CylinderList(){sz=0;}
	void insertCylinder(Cylinder c){
		CY[sz].setCoeff(c.ambCoeff,c.difCoeff,c.refCoeff,c.specCoeff,c.specExp);
		CY[sz].setColor(c.color[0],c.color[1],c.color[2]);
		CY[sz].setvalues(c.xCenter, c.yCenter, c.radius, c.zMin,c.zMax);
		sz++;
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
	}
};
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
	double l=Len(mat[2], mat[6], mat[10]);
	Front = Point(mat[2]/l, mat[6]/l, mat[10]/l);
	l=Len(mat[0], mat[4], mat[8]);
	Side = Point(mat[0]/l, mat[4]/l, mat[8]/l);
	l=Len(mat[1], mat[5], mat[9]);
	Up = Point(mat[1]/l, mat[5]/l, mat[9]/l);
}
void initCamera(){
	CameraPos.setPoint(100.0,0,20.0);
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


class RayTracer{
public:
	int filecnt;
	int width,height,recdepth,CO[3];
	double disEyePixel;
	checkerBoard CB;
	SphereList SL;
	CylinderList CL;
	TriangleList TL;
	LightSource LS;

	RayTracer(){filecnt=1;}
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
			for(i=0;i<CL.sz;i++){
				Point ip=CL.CY[i].findintersection(ref);
				if(CL.CY[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					if(tmp<mindis){
						mindis=tmp;
						minTyp=3;
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
				else if(minTyp==3){ //cylinder
					n.setPoint(CL.CY[minidx].getNormal(minHitPoint));
					colr.setPoint(CL.CY[minidx].color[0],CL.CY[minidx].color[1],CL.CY[minidx].color[2]);
					tmprefC=CL.CY[minidx].refCoeff;
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
			for(i=0;i<CL.sz;i++){
				Point ip=CL.CY[i].findintersection(ref);
				if(CL.CY[i].doesIntersect(ip)){
					double tmp= (ref.origin-ip).magnitude2();
					if(tmp<mindis){
						mindis=tmp;
						minTyp=3;
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
					tmprefC=SL.SP[minidx].specCoeff;
				}
				else if(minTyp==2){ //triangle
					n.setPoint(TL.TR[minidx].getNormal());
					colr.setPoint(TL.TR[minidx].color[0],TL.TR[minidx].color[1],TL.TR[minidx].color[2]);
					tmprefC=TL.TR[minidx].specCoeff;
				}
				else if(minTyp==3){ //cylinder
					n.setPoint(CL.CY[minidx].getNormal(minHitPoint));
					colr.setPoint(CL.CY[minidx].color[0],CL.CY[minidx].color[1],CL.CY[minidx].color[2]);
					tmprefC=CL.CY[minidx].specCoeff;
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
		for(i=0;i<CL.sz;i++){
			Point ip=CL.CY[i].findintersection(R);
			if(CL.CY[i].doesIntersect(ip)){
				CL.CY[i].setIP(ip);
				double tmp= (R.origin-ip).magnitude2();
				Obj.push_back(make_pair(i+cnt,tmp));
			}
		}
		cnt+=CL.sz;
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
			else if(Obj[0].first>=(SL.sz+TL.sz) && Obj[0].first<(SL.sz+TL.sz+CL.sz)){
				int idx=Obj[0].first-(SL.sz+TL.sz);
				Point ip=CL.CY[idx].getIP();
				Point normal=CL.CY[idx].getNormal(ip);
				Point objColor(CL.CY[idx].color[0],CL.CY[idx].color[1],CL.CY[idx].color[2]);
				Point tmp=LS.findIntensity(ip,normal,eye,CL.CY[idx].ambCoeff,
					CL.CY[idx].difCoeff,CL.CY[idx].specCoeff,CL.CY[idx].specExp/*,objColor*/);
				I.setPoint(tmp);
				Point Cf=getFinalColor(ip,R,normal,objColor,CL.CY[idx].refCoeff,0,CL.CY[idx].specCoeff);
				convertColor(I.x*Cf.x,I.y*Cf.y,I.z*Cf.z);
				//convertColor(I.x*CL.CY[idx].color[0],I.y*CL.CY[idx].color[1],I.z*CL.CY[idx].color[2]);
			}
			else if(Obj[0].first>=(SL.sz+TL.sz+CL.sz)){
				int idx=Obj[0].first-(SL.sz+TL.sz+CL.sz);
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
	for(i=0;i<rt.CL.sz;i++){
		rt.CL.CY[i].drawCylinder();	
	}
	for(i=0;i<rt.SL.sz;i++){
		rt.SL.SP[i].drawSphere();	
	}
	for(i=0;i<rt.TL.sz;i++){
		rt.TL.TR[i].drawTriangle();
	}
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
		case '7':
			rt.drawBMPImage();
		default:
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
void parseSDL(){
	ifstream ReadFile;
	ReadFile.open("Spec.txt");
	if (ReadFile.is_open()) {
		string output;
		float r,g,b,tmp;
		Point P;
		while (!ReadFile.eof()) {
			ReadFile >> output;
			if(ReadFile.eof()) break;
			//cout<<output<<endl;
			if(output=="objStart" || output=="objEnd") continue;
			else if(output== "recDepth"){
					ReadFile>>recDepth;
					rt.setrecdepth(recDepth);
					//cout<<recDepth<<endl;
					
				}
			else if(output== "pixels"){
					ReadFile>>pixels;
					rt.setDimension(pixels,pixels);
					//cout<<pixels<<endl;
			}
			else if(output== "CHECKERBOARD"){
					string output2;
					while(!ReadFile.eof()){
						ReadFile>>output2;
						if(output2=="objEnd") break;
						else if(output2== "colorOne"){
								ReadFile>>r>>g>>b;
								//cout<<r<<" "<<g<<" "<<b<<endl;
								rt.CB.setColor(r,g,b);
							}
						else if(output2== "colorTwo"){
								ReadFile>>r>>g>>b;
								//cout<<r<<" "<<g<<" "<<b<<endl;
								rt.CB.setColor2(r,g,b);
							}
						else if(output2== "ambCoeff"){
								ReadFile>>tmp;
								rt.CB.setambCoeff(tmp);
								//cout<<CB.ambCoeff<<endl;
							}
						else if(output2== "difCoeff"){
								ReadFile>>tmp;
								rt.CB.setdifCoeff(tmp);
								//cout<<CB.difCoeff<<endl;
								
							}
						else if(output2== "refCoeff"){
								ReadFile>>tmp;
								rt.CB.setrefCoeff(tmp);
								//cout<<CB.refCoeff<<endl;
								
							}
						else if(output2== "specCoeff"){
								ReadFile>>tmp;
								rt.CB.setspecCoeff(tmp);
								//cout<<CB.specCoeff<<endl;
								
							}
						else if(output2== "specExp"){
								ReadFile>>tmp;
								rt.CB.setspecExp(tmp);
								//cout<<CB.specExp<<endl;
								
							}
					}
					
				}
			else if(output== "SPHERE"){
					string output2;
					Sphere s;
					while(!ReadFile.eof()){
						ReadFile>>output2;
						if(output2=="objEnd") break;
						else if(output2== "color"){
								ReadFile>>r>>g>>b;
								//cout<<r<<" "<<g<<" "<<b<<endl;
								s.setColor(r,g,b);
								
							}
						else if(output2== "center"){
								ReadFile>>P.x>>P.z>>P.y;
								//cout<<P.x<<" "<<P.y<<" "<<P.z<<endl;
								s.setcenter(P);
								
							}
						else if(output2== "radius"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								s.setradius(tmp);
								
							}
						else if(output2== "ambCoeff"){
								ReadFile>>tmp;
								s.setambCoeff(tmp);
								//cout<<s.ambCoeff<<endl;
								
							}
						else if(output2== "difCoeff"){
								ReadFile>>tmp;
								s.setdifCoeff(tmp);
								//cout<<s.difCoeff<<endl;
								
							}
						else if(output2== "refCoeff"){
								ReadFile>>tmp;
								s.setrefCoeff(tmp);
								//cout<<s.refCoeff<<endl;
								
							}
						else if(output2== "specCoeff"){
								ReadFile>>tmp;
								s.setspecCoeff(tmp);
								//cout<<s.specCoeff<<endl;
								
							}
						else if(output2== "specExp"){
								ReadFile>>tmp;
								s.setspecExp(tmp);
								//cout<<s.specExp<<endl;
								
						}
							
					}
					rt.SL.insertSphere(s);
				}
			else if(output== "CYLINDER"){
					string output2;
					Cylinder c;
					while(!ReadFile.eof()){
						ReadFile>>output2;
						if(output2=="objEnd") break;
						else if(output2== "color"){
								ReadFile>>r>>g>>b;
								//cout<<r<<" "<<g<<" "<<b<<endl;
								c.setColor(r,g,b);
								
							}
						else if(output2== "xCenter"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								c.setxCenter(tmp);
								
							}
						else if(output2== "zCenter"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								c.setyCenter(tmp);
								
							}
						else if(output2== "yMin"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								c.setzMin(tmp);
								
							}
						else if(output2== "yMax"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								c.setzMax(tmp);
								
							}
						else if(output2== "radius"){
								ReadFile>>tmp;
								//cout<<tmp<<endl;
								c.setradius(tmp);
								
							}
						else if(output2== "ambCoeff"){
								ReadFile>>tmp;
								c.setambCoeff(tmp);
								//cout<<c.ambCoeff<<endl;
								
							}
						else if(output2== "difCoeff"){
								ReadFile>>tmp;
								c.setdifCoeff(tmp);
								//cout<<c.difCoeff<<endl;
								
							}
						else if(output2== "refCoeff"){
								ReadFile>>tmp;
								c.setrefCoeff(tmp);
								//cout<<c.refCoeff<<endl;
								
							}
						else if(output2== "specCoeff"){
								ReadFile>>tmp;
								c.setspecCoeff(tmp);
								//cout<<c.specCoeff<<endl;
								
							}
						else if(output2== "specExp"){
								ReadFile>>tmp;
								c.setspecExp(tmp);
								//cout<<c.specExp<<endl;
							
							}
					}
					rt.CL.insertCylinder(c);
				}
			else if(output== "TRIANGLE"){
					string output2;
					Triangle t;
					while(!ReadFile.eof()){
						ReadFile>>output2;
						if(output2=="objEnd") break;
						else if(output2== "a"){
								ReadFile>>P.x>>P.z>>P.y;
								t.seta(P);
								//cout<<P.x<<" "<<P.y<<" "<<P.z<<endl;
								
							}
						  else if(output2== "b"){
								ReadFile>>P.x>>P.z>>P.y;
								t.setb(P);
								//cout<<P.x<<" "<<P.y<<" "<<P.z<<endl;
								
							}
						   else if(output2== "c"){
								ReadFile>>P.x>>P.z>>P.y;
								t.setc(P);
								//cout<<P.x<<" "<<P.y<<" "<<P.z<<endl;
								
							}
						   else if(output2== "color"){
								ReadFile>>r>>g>>b;
								//cout<<r<<" "<<g<<" "<<b<<endl;
								t.setColor(r,g,b);
								
							}
						   else if(output2=="ambCoeff"){
								ReadFile>>tmp;
								t.setambCoeff(tmp);
								//cout<<t.ambCoeff<<endl;
								
							}
							else if(output2== "difCoeff"){
								ReadFile>>tmp;
								t.setdifCoeff(tmp);
								//cout<<t.difCoeff<<endl;
								
							}
							else if(output2=="refCoeff"){
								ReadFile>>tmp;
								t.setrefCoeff(tmp);
								//cout<<t.refCoeff<<endl;
								
							}
							else if(output2=="specCoeff"){
								ReadFile>>tmp;
								t.setspecCoeff(tmp);
								//cout<<t.specCoeff<<endl;
								
							}
							else if(output2== "specExp"){
								ReadFile>>tmp;
								t.setspecExp(tmp);
								//cout<<t.specExp<<endl;
								
							}
							else if(output2== "refractiveIndex"){
								ReadFile>>tmp;
								t.setrefIdx(tmp);
								//cout<<t.refIdx<<endl;
								
							}
					}
					rt.TL.insertTriangle(t);
				}
				else if(output== "light"){
					Point P;
					ReadFile>>P.x>>P.z>>P.y;
					rt.LS.insertLight(P.x,P.y,P.z);
					//cout<<P.x<<" "<<P.y<<" "<<P.z<<endl;
					
				}
			
		}
	}
	ReadFile.close();
	cout<<"file parsed successfully"<<endl;
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
	parseSDL();
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
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
