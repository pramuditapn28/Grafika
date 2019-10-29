#include <math.h>
#include<windows.h>
#include <GL/glut.h>
//#include <fstream.h>
#include <stdlib.h>

typedef struct {
	float m[4][4];
} matrix3D_t;

typedef struct {
	float v[4];
} vector3D_t;

typedef struct {
	float x;
	float y;
	float z;
} point3D_t;

typedef struct {
	float x;
	float y;
} point2D_t;

typedef struct {
	float r;
	float g;
	float b;
} color_t;

////////////////// matrices and vectors 3D ver 2 /////////////////
matrix3D_t createIdentity(void){
	matrix3D_t u;
	int i,j;
	for (i=0;i<4;i++) {
		for(j=0;j<4;j++) u.m[i][j]=0.;
		u.m[i][i]=1.;
	}
	return u;
}

matrix3D_t operator * (matrix3D_t a,matrix3D_t b){
	matrix3D_t c;//c=a*b
	int i,j,k;
	for (i=0;i<4;i++) for (j=0;j<4;j++) {
		c.m[i][j]=0;
		for (k=0;k<4;k++) c.m[i][j]+=a.m[i][k]*b.m[k][j];
	}
	return c;
}

vector3D_t operator * (matrix3D_t a, vector3D_t b){
	vector3D_t c;//c=a*b
	int i,j;
	for (i=0;i<4;i++) {
		c.v[i]=0;
		for (j=0;j<4;j++) c.v[i]+=a.m[i][j]*b.v[j];
	}
	return c;
}

matrix3D_t translationMTX(float dx,float dy,float dz){
	matrix3D_t trans=createIdentity();
	trans.m[0][3]=dx;
	trans.m[1][3]=dy;
	trans.m[2][3]=dz;
	return trans;
}

matrix3D_t rotationXMTX(float theta){
	matrix3D_t rotate=createIdentity();
	float cs=cos(theta);
	float sn=sin(theta);
	rotate.m[1][1]=cs; rotate.m[1][2]=-sn;
	rotate.m[2][1]=sn; rotate.m[2][2]=cs;
	return rotate;
}

matrix3D_t rotationYMTX(float theta){
	matrix3D_t rotate=createIdentity();
	float cs=cos(theta);
	float sn=sin(theta);
	rotate.m[0][0]=cs; rotate.m[0][2]=sn;
	rotate.m[2][0]=-sn; rotate.m[2][2]=cs;
	return rotate;
}

matrix3D_t rotationZMTX(float theta){
	matrix3D_t rotate=createIdentity();
	float cs=cos(theta);
	float sn=sin(theta);
	rotate.m[0][0]=cs; rotate.m[0][1]=-sn;
	rotate.m[1][0]=sn; rotate.m[1][1]=cs;
	return rotate;
}

matrix3D_t scalingMTX(float factorx,float factory,float factorz){
	matrix3D_t scale=createIdentity();
	scale.m[0][0]=factorx;
	scale.m[1][1]=factory;
	scale.m[2][2]=factorz;
	return scale;
}

matrix3D_t perspectiveMTX(float eyelength){
	matrix3D_t perspective=createIdentity();
	perspective.m[3][2]=-1./eyelength;
	return perspective;
}

point2D_t Vector2Point2D(vector3D_t vec){
	point2D_t pnt;
	pnt.x=vec.v[0];
	pnt.y=vec.v[1];
	return pnt;
}

point3D_t Vector2Point3D(vector3D_t vec){
	point3D_t pnt;
	pnt.x=vec.v[0];
	pnt.y=vec.v[1];
	pnt.z=vec.v[2];
	return pnt;
}

vector3D_t Point2Vector(point3D_t pnt){
	vector3D_t vec;
	vec.v[0]=pnt.x;
	vec.v[1]=pnt.y;
	vec.v[2]=pnt.z;
	vec.v[3]=1.;
	return vec;
}

vector3D_t homogenizeVector(vector3D_t vec){
	int i;
	for (i=0;i<3;i++) {
		vec.v[i]/=vec.v[3];
	}
	vec.v[3]=1.;
	return vec;
}

vector3D_t unitVector(vector3D_t vec){
	int i;
	float vec2=0.;
	float vec1,invvec1;
	for (i=0;i<3;i++) {
		vec2+=vec.v[i]*vec.v[i];
	}
	vec1=sqrt(vec2);
	if (vec1!=0.) {
		invvec1=1./vec1;
		for (i=0;i<3;i++) {
			vec.v[i]*=invvec1;
		}
	}
	vec.v[3]=1.;
	return vec;
}

// inner product (dot product) of homogeneous vector
float operator * (vector3D_t a, vector3D_t b){
	float c;//c=a*b
	int i;
	c=0;
	for (i=0;i<3;i++) {
		c+=a.v[i]*b.v[i];
	}
	return c;
}

// outer product (cross product ) of homogeneous vector
//       i         j         k
//       a0       a1        a2
//       b0       b1        b2
vector3D_t operator ^ (vector3D_t a, vector3D_t b){
	vector3D_t c;//c=a*b
	c.v[0]=a.v[1]*b.v[2]-a.v[2]*b.v[1];
	c.v[1]=a.v[2]*b.v[0]-a.v[0]*b.v[2];
	c.v[2]=a.v[0]*b.v[1]-a.v[1]*b.v[0];
	c.v[3]=1.;
	return c;
}

vector3D_t operator - (vector3D_t v1,vector3D_t v0){
	vector3D_t c;//c=v1-v0
	c.v[0]=v1.v[0]-v0.v[0];
	c.v[1]=v1.v[1]-v0.v[1];
	c.v[2]=v1.v[2]-v0.v[2];
	c.v[3]=1.;
	return c;
}

vector3D_t operator - (vector3D_t v){
	vector3D_t c;//c=-v
	c.v[0]=-v.v[0];
	c.v[1]=-v.v[1];
	c.v[2]=-v.v[2];
	c.v[3]=1.;
	return c;
}

vector3D_t operator * (float r, vector3D_t b){
	vector3D_t c;//c=r*b
	int i;
	for (i=0;i<3;i++) {
		c.v[i]=r*b.v[i];
	}
	c.v[3]=1.;
	return c;
}

vector3D_t operator * (vector3D_t b, float r){
	vector3D_t c;//c=r*b
	int i;
	for (i=0;i<3;i++) {
		c.v[i]=r*b.v[i];
	}
	c.v[3]=1.;
	return c;
}

float funcPositive(float x){
	if (0.<x) return x;
	else return 0.;
}

// x to yth power
float power(float x,float y){
	//ln z = y ln x        z = exp (y ln x)
	if (x==0.) return 0;
	return exp(y*log(x));
}

color_t operator + (color_t c1, color_t c2){
	color_t col;
	col.r=c1.r+c2.r;
	col.g=c1.g+c2.g;
	col.b=c1.b+c2.b;
	return col;
}

color_t operator * (float r, color_t c){
	color_t col;
	col.r=r*c.r;
	col.g=r*c.g;
	col.b=r*c.b;
	return col;
}

color_t operator * (color_t c, float r){
	color_t col;
	col.r=r*c.r;
	col.g=r*c.g;
	col.b=r*c.b;
	return col;
}

//PhongModel color calculation
// LightVector, NormalVector, ViewVector, ColorofObject
color_t PhongModel(vector3D_t Light,vector3D_t Normal,vector3D_t View,color_t col){
	float kspe=0.8; // specular reflection coefficient
	float kdif=0.5; // diffuse reflection coefficient
	float kamb=0.3; // ambient light coefficient
	float tmp,NL,RV;
	color_t ColWhite={1,0,0};
	vector3D_t ReflectionVector=(2.*(Light*Normal)*Normal)-Light;
	tmp=Normal*Light;
	NL=funcPositive(tmp);
	tmp=ReflectionVector*View;
	RV=funcPositive(tmp);
	return kdif*NL*col+kspe*power(RV,4)*ColWhite+kamb*col;
	//return kdif*NL*col+kamb*col;
}

////////////// End of matrices and vectors 3D ver 2 //////////////

////////////// OpenGL drawShape Functions ver 1 /////////////////
void setColor(float red,float green,float blue){
	glColor3f(red, green, blue);
}

void setColor(color_t col){
	glColor3f(col.r, col.g, col.b);
}

void drawDot(float x,float y){
	glBegin(GL_POINTS);
		glVertex2f(x, y);
	glEnd();
}

void drawLine(float x1, float y1, float x2, float y2){
	glBegin(GL_LINES);
		glVertex2f(x1, y1);
		glVertex2f(x2, y2);
	glEnd();
}

void drawLine(point2D_t p1,point2D_t p2){
	drawLine(p1.x,p1.y,p2.x,p2.y);
}

//n: number of points
void drawPolyline(point2D_t pnt[],int n){
	int i;
	glBegin(GL_LINE_STRIP);
		for (i=0;i<n;i++) {
			glVertex2f(pnt[i].x, pnt[i].y);
		}
	glEnd();
}

//n: number of vertices
void drawPolygon(point2D_t pnt[],int n){
	int i;
	glBegin(GL_LINE_LOOP);
		for (i=0;i<n;i++) {
			glVertex2f(pnt[i].x, pnt[i].y);
		}
	glEnd();
}

// The function fillPolygon can fills only convex polygons
//n: number of vertices
void fillPolygon(point2D_t pnt[],int n,color_t color){
	int i;
	setColor(color);
	glBegin(GL_POLYGON);
		for (i=0;i<n;i++) {
			glVertex2f(pnt[i].x, pnt[i].y);
		}
	glEnd();
}

// The function gradatePolygon can fills only convex polygons
// The vertices will be painted with corresponding given colors.
// The points inside the polygon will be painted with the mixed color.
//n: number of vertices
void gradatePolygon(point2D_t pnt[],int num,color_t col[]){
	int i;
	glBegin(GL_POLYGON);
		for (i=0;i<num;i++) {
			setColor(col[i]);
			glVertex2f(pnt[i].x, pnt[i].y);
		}
	glEnd();
}

//////////// End of OpenGL drawShape Functions ver 1 ////////////

void userdraw(void);

void display(void){
	glClear( GL_COLOR_BUFFER_BIT);
	userdraw();
	glutSwapBuffers();
}

//////////////////////////////////////////////////////////////////
void drawcharX(float x,float y){
	drawLine(x,y,x+10,y+12);drawLine(x,y+12,x+10,y);
}

void drawcharY(float x,float y){
	drawLine(x+5,y,x+5,y+7);drawLine(x,y+12,x+5,y+7);drawLine(x+10,y+12,x+5,y+7);
}

void drawcharZ(float x,float y){
	drawLine(x,y+12,x+10,y+12);drawLine(x+10,y+12,x,y);drawLine(x,y,x+10,y);
}

void drawAxes(matrix3D_t view){
#define HALFAXIS  220
#define HALFAXIS1 (HALFAXIS-10)
	point3D_t axes[14]={
		{-HALFAXIS,0,0},{HALFAXIS,0,0},{HALFAXIS1,5,0},{HALFAXIS1,0,0},{0,0,0},
		{0,-HALFAXIS,0},{0,HALFAXIS,0},{0,HALFAXIS1,5},{0,HALFAXIS1,0},{0,0,0},
		{0,0,-HALFAXIS},{0,0,HALFAXIS},{5,0,HALFAXIS1},{0,0,HALFAXIS1}
	};
	vector3D_t vec[14];
	point2D_t buff[14];
	int i;
	for (i=0;i<14;i++) {
		vec[i]=Point2Vector(axes[i]);
		vec[i]=view*vec[i];
		buff[i]=Vector2Point2D(vec[i]);
	}
	drawPolyline(buff,14);
	drawcharX(buff[1].x,buff[1].y);
	drawcharY(buff[6].x,buff[6].y);
	drawcharZ(buff[11].x-14,buff[11].y);
}

//////////////////////////////////////////////////////////////////
typedef struct {
	int NumberofVertices; //in the face
	short int pnt[50];
	color_t col;
} face_t;
typedef struct {
	int NumberofVertices; //of the object
	point3D_t pnt[1600];
	color_t col[1600];
	int NumberofFaces; //of the object
	face_t fc[1000];
} object3D_t;

float zRata(object3D_t obyek,matrix3D_t mat){
	int i;
	float z=0;
	vector3D_t vec;

	for(i=0;i<obyek.NumberofVertices;i++){
		vec=Point2Vector(obyek.pnt[i]);
		vec=mat*vec;
		z=z+vec.v[2];
	}
	z=z/obyek.NumberofVertices;
	return z;
}


void draw3D(object3D_t obyek,matrix3D_t mat){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	point2D_t p[50];
	int i,j;
	for(i=0;i<obyek.NumberofVertices;i++){
		vec[i]=Point2Vector(obyek.pnt[i]);
		vec[i]=mat*vec[i];
	}
	setColor(1,0,0);
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]<0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			drawPolygon(p,obyek.fc[i].NumberofVertices);
		}
	}
	setColor(1,1,1);
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]>=0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			drawPolygon(p,obyek.fc[i].NumberofVertices);
		}
	}
}

void draw3Da(object3D_t obyek,matrix3D_t mat){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	point2D_t p[50];
	int i,j;
	for(i=0;i<obyek.NumberofVertices;i++){
		vec[i]=Point2Vector(obyek.pnt[i]);
		vec[i]=mat*vec[i];
	}
	setColor(1,1,1);
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			p[j]=Vector2Point2D(vecbuff[j]);
		drawPolygon(p,obyek.fc[i].NumberofVertices);
	}
}

void draw3Dw(object3D_t obyek,matrix3D_t mat,color_t col){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	vector3D_t lightVector={0,0,1,1},viewVector={0,0,1,1};
	color_t colbuff;
	point2D_t p[50];
	int i,j;
	for(i=0;i<obyek.NumberofVertices;i++){
		vec[i]=Point2Vector(obyek.pnt[i]);
		vec[i]=mat*vec[i];
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]<0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			vecNormal=unitVector(vecNormal);
			colbuff=PhongModel(lightVector,vecNormal,viewVector,col);
			fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]>=0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			vecNormal=unitVector(vecNormal);
			colbuff=PhongModel(lightVector,vecNormal,viewVector,col);
			fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
}

void draw3Dc(object3D_t obyek,matrix3D_t mat){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	vector3D_t lightVector={0,0,1,1},viewVector={0,0,1,1};
	color_t colbuff;
	point2D_t p[50];
	int i,j;
	for(i=0;i<obyek.NumberofVertices;i++){
		vec[i]=Point2Vector(obyek.pnt[i]);
		vec[i]=mat*vec[i];
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]<0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			vecNormal=unitVector(vecNormal);
			colbuff=PhongModel(lightVector,vecNormal,viewVector,obyek.fc[i].col);
			fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]>=0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
			}
			vecNormal=unitVector(vecNormal);
			colbuff=PhongModel(lightVector,vecNormal,viewVector,obyek.fc[i].col);
			fillPolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
}

void draw3Dcg(object3D_t obyek,matrix3D_t mat){
	vector3D_t vec[1600], vecbuff[50];
	vector3D_t vecNormal;
	vector3D_t lightVector={0,0,1,1},viewVector={0,0,1,1};
	color_t colbuff[50];
	point2D_t p[50];
	int i,j;
	for(i=0;i<obyek.NumberofVertices;i++){
		vec[i]=Point2Vector(obyek.pnt[i]);
		vec[i]=mat*vec[i];
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]<0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
				vecNormal=unitVector(vecbuff[j]);
				colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
			}
			gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
	for(i=0;i<obyek.NumberofFaces;i++){
		for(j=0;j<obyek.fc[i].NumberofVertices;j++)
			vecbuff[j]=vec[obyek.fc[i].pnt[j]];
		vecNormal=(vecbuff[1]-vecbuff[0])^(vecbuff[2]-vecbuff[0]);
		if(vecNormal.v[2]>=0){
			for(j=0;j<obyek.fc[i].NumberofVertices;j++){
				p[j]=Vector2Point2D(vecbuff[j]);
				vecNormal=unitVector(vecbuff[j]);
				colbuff[j]=PhongModel(lightVector,vecNormal,viewVector,obyek.col[obyek.fc[i].pnt[j]]);
			}
			gradatePolygon(p,obyek.fc[i].NumberofVertices,colbuff);
		}
	}
}

void makeCube(object3D_t &kubus, float d){
	object3D_t obyek={8,
	{{0,0,0},{d,0,0},{d,0,d},{0,0,d},{0,d,0},{d,d,0},{d,d,d},{0,d,d}},
	{{1,1,0},{1,1,0},{1,1,0},{1,1,0},{1,1,0},{1,1,0},{1,1,0}},
	6,
	{{4,{0,1,2,3},{1,1,0}},{4,{4,7,6,5},{1,1,0}},{4,{1,5,6,2},{1,1,0}},{4,{0,3,7,4},{1,1,0}},{4,{2,6,7,3},{1,1,0}},{4,{0,4,5,1},{1,1,0}}}};
	matrix3D_t mat=translationMTX(-d/2,-d/2,-d/2);
	vector3D_t vec;
	kubus=obyek;
	for(int i=0;i<kubus.NumberofVertices;i++){
		vec=Point2Vector(kubus.pnt[i]);
		vec=mat*vec;
		kubus.pnt[i]=Vector2Point3D(vec);
	}
}

void makeCylinder(object3D_t &silinder, int n, float r, float h){
	float a=6.28/n;
	int i;
	for(i=0;i<n;i++){
		silinder.pnt[i].x=r*cos(i*a);
		silinder.pnt[i].y=0;
		silinder.pnt[i].z=r*sin(i*a);
		silinder.pnt[n+i].x=r*cos(i*a);
		silinder.pnt[n+i].y=h;
		silinder.pnt[n+i].z=r*sin(i*a);
	}
	silinder.NumberofVertices=2*n;
	for(i=0;i<n;i++){
		silinder.fc[i].NumberofVertices=4;
		silinder.fc[i].pnt[0]=i;
		silinder.fc[i].pnt[1]=n+i;
		silinder.fc[i].pnt[2]=n+i+1;
		silinder.fc[i].pnt[3]=i+1;
		if(i==(n-1)){
			silinder.fc[i].pnt[2]=n;
			silinder.fc[i].pnt[3]=0;
		}
	}
	silinder.fc[n].NumberofVertices=n;
	for(i=0;i<n;i++) silinder.fc[n].pnt[i]=i;
	silinder.fc[n+1].NumberofVertices=n;
	for(i=0;i<n;i++) silinder.fc[n+1].pnt[i]=2*n-1-i;
	silinder.NumberofFaces=n+2;
	color_t c={1,1,0};
	for(i=0;i<silinder.NumberofFaces;i++) silinder.fc[i].col=c;
	for(i=0;i<silinder.NumberofVertices;i++)
		silinder.col[i]=c;
}

void makeCylinderAwur(object3D_t &silinder, int n, float r, float h){
	float a=6.28/n;
	int i;
	float f;
	for(i=0;i<n;i++){
		f=3*i/n;
		silinder.pnt[i].x=r*f*cos(i*a);
		silinder.pnt[i].y=0;
		silinder.pnt[i].z=r*f*sin(i*a);
		silinder.pnt[n+i].x=r*f*cos(i*a);
		silinder.pnt[n+i].y=h;
		silinder.pnt[n+i].z=r*f*sin(i*a);
	}
	silinder.NumberofVertices=2*n;
	for(i=0;i<n;i++){
		silinder.fc[i].NumberofVertices=4;
		silinder.fc[i].pnt[0]=i;
		silinder.fc[i].pnt[1]=n+i;
		silinder.fc[i].pnt[2]=n+i+1;
		silinder.fc[i].pnt[3]=i+1;
		if(i==(n-1)){
			silinder.fc[i].pnt[2]=n;
			silinder.fc[i].pnt[3]=0;
		}
	}
	silinder.fc[n].NumberofVertices=n;
	for(i=0;i<n;i++) silinder.fc[n].pnt[i]=i;
	silinder.fc[n+1].NumberofVertices=n;
	for(i=0;i<n;i++) silinder.fc[n+1].pnt[i]=2*n-1-i;
	silinder.NumberofFaces=n+2;
	color_t c={1,1,0};
	for(i=0;i<silinder.NumberofFaces;i++) silinder.fc[i].col=c;
	for(i=0;i<silinder.NumberofVertices;i++)
		silinder.col[i]=c;
}

void makeCone(object3D_t &kerucut, int n, float r,float h){
	float a=6.28/n;
	int i;
	kerucut.pnt[0].x=0;
	kerucut.pnt[0].y=h;
	kerucut.pnt[0].z=0;
	for(i=1;i<=n;i++){
		kerucut.pnt[i].x=r*cos(i*a);
		kerucut.pnt[i].y=0;
		kerucut.pnt[i].z=r*sin(i*a);
	}
	for(i=0;i<n;i++){
		kerucut.fc[i].NumberofVertices=3;
		kerucut.fc[i].pnt[0]=0;
		kerucut.fc[i].pnt[1]=i+2;
		kerucut.fc[i].pnt[2]=i+1;
		if(i==(n-1)) kerucut.fc[i].pnt[1]=1;
	}
	kerucut.fc[n].NumberofVertices=n;
	for(i=0;i<n;i++) kerucut.fc[n].pnt[i]=i+1;
	kerucut.NumberofVertices=n+1;
	kerucut.NumberofFaces=n+1;
	color_t c={1,1,0};
	for(i=0;i<kerucut.NumberofFaces;i++) kerucut.fc[i].col=c;
	for(i=0;i<kerucut.NumberofVertices;i++)
		kerucut.col[i]=c;
}

void makeCylinderN(object3D_t &silinder,int m,int n,float r[],float h[],int sw){
	float a=6.26/n;
	float b=0;
	int i,j;
	silinder.NumberofVertices=(m+1)*n;
	for(i=0;i<=m;i++){
		if(i>0) b=b+h[i-1];
		for(j=0;j<n;j++){
			silinder.pnt[i*n+j].x=r[i]*cos(j*a);
			silinder.pnt[i*n+j].y=b;
			silinder.pnt[i*n+j].z=r[i]*sin(j*a);
		}
	}
	silinder.NumberofFaces=m*n+2;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			silinder.fc[i*n+j].NumberofVertices=4;
			silinder.fc[i*n+j].pnt[0]=i*n+j;
			silinder.fc[i*n+j].pnt[1]=(i+1)*n+j;
			silinder.fc[i*n+j].pnt[2]=(i+1)*n+j+1;
			silinder.fc[i*n+j].pnt[3]=i*n+j+1;
			if(j==(n-1)){
				silinder.fc[i*n+j].pnt[2]=i*n+j+1;
				silinder.fc[i*n+j].pnt[3]=(i-1)*n+j+1;
			}
		}
	}
	if(sw==0 || sw==1){
		silinder.fc[m*n].NumberofVertices=n;
		for(i=0;i<n;i++) silinder.fc[m*n].pnt[i]=i;
	}
	if(sw==0 || sw==2){
		silinder.fc[m*n+1].NumberofVertices=n;
		for(i=0;i<n;i++) silinder.fc[m*n+1].pnt[i]=(m+1)*n-1-i;
	}
	color_t c={1,1,0};
	for(i=0;i<silinder.NumberofFaces;i++) silinder.fc[i].col=c;
	for(i=0;i<silinder.NumberofVertices;i++)
		silinder.col[i]=c;
}

void makeSphere(object3D_t &sphere,int n,float r){
	float a=6.28/n;
	float b=6.28/n;
	int i,j;
	sphere.NumberofVertices=(n+1)*n;
	for(i=0;i<=n;i++){
		for(j=0;j<n;j++){
			sphere.pnt[i*n+j].x=r*cos(j*a)*sin(i*b);
			sphere.pnt[i*n+j].y=r*cos(i*b);
			sphere.pnt[i*n+j].z=r*sin(j*a)*sin(i*b);
		}
	}
	sphere.NumberofFaces=n*n+2;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			sphere.fc[i*n+j].NumberofVertices=4;
			sphere.fc[i*n+j].pnt[0]=i*n+j;
			sphere.fc[i*n+j].pnt[1]=(i+1)*n+j;
			sphere.fc[i*n+j].pnt[2]=(i+1)*n+j+1;
			sphere.fc[i*n+j].pnt[3]=i*n+j+1;
			if(j==(n-1)){
				sphere.fc[i*n+j].pnt[2]=i*n+j+1;
				sphere.fc[i*n+j].pnt[3]=(i-1)*n+j+1;
			}
		}
	}
	sphere.fc[n*n].NumberofVertices=n;
	for(i=0;i<n;i++) sphere.fc[n*n].pnt[i]=i;
	sphere.fc[n*n+1].NumberofVertices=n;
	for(i=0;i<n;i++) sphere.fc[n*n+1].pnt[i]=(n+1)*n-1-i;
	color_t c={1,1,0};
	for(i=0;i<sphere.NumberofFaces;i++) sphere.fc[i].col=c;
	for(i=0;i<sphere.NumberofVertices;i++)
		sphere.col[i]=c;
}

void makeApple(object3D_t &sphere,int n,float r){
	float a=6.28/n;
	float b=6.28/n;
	int i,j,k;
	sphere.NumberofVertices=(n+1)*n;
	for(i=0;i<=n;i++){
		for(j=0;j<n;j++){
			sphere.pnt[i*n+j].x=r*cos(j*a)*sin(i*b);
			sphere.pnt[i*n+j].y=r*cos(i*b);
			if(i>n-1){
				k=i-(n-1);
				sphere.pnt[i*n+j].y-=20*k;
			}
			sphere.pnt[i*n+j].z=r*sin(j*a)*sin(i*b);
		}
	}
	sphere.NumberofFaces=n*n+2;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			sphere.fc[i*n+j].NumberofVertices=4;
			sphere.fc[i*n+j].pnt[0]=i*n+j;
			sphere.fc[i*n+j].pnt[1]=(i+1)*n+j;
			sphere.fc[i*n+j].pnt[2]=(i+1)*n+j+1;
			sphere.fc[i*n+j].pnt[3]=i*n+j+1;
			if(j==(n-1)){
				sphere.fc[i*n+j].pnt[2]=i*n+j+1;
				sphere.fc[i*n+j].pnt[3]=(i-1)*n+j+1;
			}
		}
	}
	sphere.fc[n*n].NumberofVertices=n;
	for(i=0;i<n;i++) sphere.fc[n*n].pnt[i]=i;
	sphere.fc[n*n+1].NumberofVertices=n;
	for(i=0;i<n;i++) sphere.fc[n*n+1].pnt[i]=(n+1)*n-1-i;
	color_t c={1,1,0};
	for(i=0;i<sphere.NumberofFaces;i++) sphere.fc[i].col=c;
	for(i=0;i<sphere.NumberofVertices;i++)
		sphere.col[i]=c;
}

void makeCheese(object3D_t &o,int n,float d){
	o.NumberofVertices=(n+1)*(n+1);
	int i,j,k;
	for(i=0;i<n+1;i++)
		for(j=0;j<n+1;j++){
			k=(n+1)*i+j;
			o.pnt[k].x=j*d;
			o.pnt[k].y=0;
			o.pnt[k].z=i*d;
		}
	o.NumberofFaces = n*n;
	color_t w1={1,1,1},w2={0,0,0};
	for(i=0;i<n;i++)
		for(j=0;j<n;j++){
			k=n*i+j;
			o.fc[k].NumberofVertices=4;
			o.fc[k].pnt[0]=(n+1)*i+j;
			o.fc[k].pnt[1]=(n+1)*i+j+n+1;
			o.fc[k].pnt[2]=(n+1)*i+j+n+2;
			o.fc[k].pnt[3]=(n+1)*i+j+1;
			if((i+j)%2==0) o.fc[k].col=w1;
			else o.fc[k].col=w2;
		}
	for(i=0;i<o.NumberofVertices;i++)
		o.col[i]=w1;
}

point3D_t interpolate(point3D_t p1,point3D_t p2,float a){
	point3D_t p;
	p.x=(1-a)*p1.x+a*p2.x;
	p.y=(1-a)*p1.y+a*p2.y;
	p.z=(1-a)*p1.z+a*p2.z;
	return p;
}

void userdraw(void){
    static float a=0;
    setColor(1,1,1);
    matrix3D_t tilting = rotationYMTX(a)*rotationXMTX(0.5);
    setColor(1,1,0);
    drawAxes(tilting);

    /*
    //draw points
    point3D_t p={100,100,50};
    vector3D_t vec;
    point2D_t p2;
    vec= Point2Vector(p);
    vec=tilting*vec;
    p2=Vector2Point2D(vec);
    glPointSize(4);
    glBegin(GL_POINTS);
    glVertex2f(p2.x,p2.y);
    glEnd();

    //draw line
    point3D_t q[2]= {{100,0,0},{0,100,0}};
    point2D_t q2[2];
    for(int i=0;i<2;i++){
        vec= Point2Vector(q[i]);
        vec=tilting*vec;
        q2[i]=Vector2Point2D(vec);
        drawLine(q2[0],q2[1]);
    }

    //draw triangle
    point3D_t qs[3]= {{0,0,0},{-100,0,0},{0,-100,0}};
    point2D_t qs2[3];
    for(int i=0;i<3;i++){
        vec= Point2Vector(qs[i]);
        vec=tilting*vec;
        qs2[i]=Vector2Point2D(vec);
    }
    setColor(1,0,0);
    glBegin(GL_TRIANGLES);
    glVertex2f(qs2[0].x,qs2[0].y);
    glVertex2f(qs2[1].x,qs2[1].y);
    glVertex2f(qs2[2].x,qs2[2].y);
    glEnd();

    //draw cube wirh function
    object3D_t obj;
    makeCube(obj,50);
    draw3D(obj,tilting);
    */
    object3D_t kubus={8,
        {{0,0,0},{100,0,0},{100,0,100},{0,0,100},
            {0,100,0},{100,100,0},{100,100,100},
            {0,100,100}},
        {{1,1,1},{1,1,1},{1,1,1},{1,1,1},{1,1,1},
            {1,1,1},{1,1,1},{1,1,1}},
        6,
        {{4,{0,1,2,3},{1,1,1}},{4,{7,6,5,4},{1,1,1}},{4,{6,7,3,2},{1,1,1}},
            {4,{4,5,1,0},{1,1,1}},{4,{0,3,7,4},{1,1,1}},{4,{2,1,5,6},{1,1,1}}}
    };
    setColor(1,1,1);
    draw3Dcg(kubus,tilting);
    a+=0.001;

}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB );
	glutInitWindowPosition(100,100);
	glutInitWindowSize(640,480);
	glutCreateWindow ("cube");
	glClearColor(0.0, 0.0, 0.0, 0.0);
	gluOrtho2D(-320., 320., -240.0, 240.0);
	  // Define the dimensions of the Orthographic Viewing Volume
	glutIdleFunc(display); // idle event call back
	glutDisplayFunc(display);
	glutMainLoop();
	return 0;
}
