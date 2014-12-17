 

 
 

#include <windows.h>
// OpenGL Graphics includes
#include <GL/glew.h>
#include <GL/freeglut.h>


// includes, system
#include<stdlib.h>
#include <stdio.h> 
#include <math.h>

// includes, project
#include <cutil_inline.h>
//#include <shrQATest.h>
#include<iostream>
using namespace std;
 
#include "KDTree.h"
#include "Mesh.h"

#define WIDTH 200
#define HEIGH 200
#define INF 100000			
#define EPSILON 0.001
#define PI 3.14159265

#define PHON 30
#define AMBIENT 0.8
#define SPECULAR 0.0
#define DIFFUSE 0.8
#define RGB_R 1.0
#define RGB_G 1.0
#define RGB_B 1.0
#define GLUT_WHEEL_UP 3
#define GLUT_WHEEL_DOWN 4



double anger_x = 0;
double anger_y = 0;
double anger_move_x = 0;
double anger_move_y = 0;

double shift_x = 0;
double shift_y = 0;
double shift_move_x = 0;
double shift_move_y = 0;

double scale_x = 0;
double sacle_y = 0;

double scale = 50;
bool left_tag = false;
bool right_tag = false;
bool middle_tag = false;
int mouse_x = 0;
int mouse_y = 0;
void processMouseActiveMotion(int x, int y);
void processMouse(int button, int state, int x, int y);

void keyboard(unsigned char key,int x, int y);
void draw();
void drawGraph();
void drawBox(Point& p1,Point& p2);


KDTree tree;
Mesh mesh;
GLubyte* pixelData;
int index_node=0;
KDNode* kdnode;

struct Stack{
	int index;
	bool leafNode;
};

struct Ray
{
	Point pos;
	Point dir;

};
struct Light
{
	Point pos;
	Point col;

};

void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(WIDTH,HEIGH,GL_BGR_EXT,GL_UNSIGNED_BYTE,pixelData);
	glutSwapBuffers();
}

__device__ float Dot(Point& a, Point& b)
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}

__device__ Point CrossProduct(Point& a, Point& b)
{
	Point ret;
	ret.x = a.y * b.z - a.z * b.y;
	ret.y = a.z * b.x - a.x * b.z;
	ret.z = a.x * b.y - a.y * b.x;
	return ret;
}

__device__ void Normalize(Point& vector)
{
	float v=sqrt(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
	vector.x /= v; vector.y /= v; vector.z /= v;
}
__device__ Point minus(Point& p1, Point& p2)
{
	Point ret;
	ret.x = p1.x - p2.x;
	ret.y = p1.y - p2.y;
	ret.z = p1.z - p2.z;

	return ret;

}



__device__ Point multi(Point& p1, Point& p2)  
{
	Point ret;
	ret.x = p1.x*p2.x;
	ret.y = p1.y*p2.y;
	ret.z = p1.z*p2.z;
	return ret;
}

__device__ Point make(float a, float b, float c)
{
	Point ret;
	ret.x = a;
	ret.y = b;
	ret.z= c;
	return ret;
}

__device__ Point add(Point& p1, Point& p2)
{
	Point ret;
	ret.x = p1.x+p2.x;
	ret.y = p1.y+p2.y;
	ret.z = p1.z+p2.z;
	return ret;

}

__device__ Point getNormal(Point& rkPoint0, Point& rkPoint1, Point& rkPoint2)
{
	Point kEdge1 = minus(rkPoint1 , rkPoint0);
	Point kEdge2 = minus(rkPoint2 , rkPoint0);
	Point normal = CrossProduct(kEdge1 , kEdge2);
	Normalize(normal);
	return normal;
}

__device__ bool pointInTriangle(Point& a, Point& b, Point& c, Point& p)
{
	Point AB = minus(b , a);
	Point AC = minus(c , a);
	Point AP = minus(p , a);
	Point BC = minus(c ,  b);
	Point BP = minus(p , b);
	Point temp;
	temp = CrossProduct(AB, AC);
	float left = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
	temp = CrossProduct(AB, AP);

	float right = 0;//distance(CrossProduct(AB, AP), temp)+ distance(CrossProduct(AP, AC),temp) + distance(CrossProduct(BC, BP),temp);
	right += sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
	temp = CrossProduct(AP, AC);
	right += sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
	temp = CrossProduct(BC, BP);
	right += sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);


	return (left-EPSILON<right)&&(right<left+EPSILON);
}

__device__ int intersectLeaf(Ray& ray, int leaf, float* d_vertice, unsigned int * d_face, KDLeafNode* d_leafNode, int* d_leafNodeIndex,Point& d_normal, Point& d_intersection)
{
	if(leaf==-1)
		return -1;
	int start = d_leafNode[leaf].begin;
	int end = d_leafNode[leaf].begin+d_leafNode[leaf].numObject;
	float minDis=10000;
	int minObjNum=-1;
	Point planeNormal;//store the normal of the intersected plane
	for(int m=start; m<end; m++)
	{
		Point PT1,PT2,PT3;

		PT1.x = d_vertice[d_face[d_leafNodeIndex[m]*3]*3];
		PT1.y = d_vertice[d_face[d_leafNodeIndex[m]*3]*3+1];
		PT1.z = d_vertice[d_face[d_leafNodeIndex[m]*3]*3+2];

		PT2.x = d_vertice[d_face[d_leafNodeIndex[m]*3+1]*3];
		PT2.y = d_vertice[d_face[d_leafNodeIndex[m]*3+1]*3+1];
		PT2.z = d_vertice[d_face[d_leafNodeIndex[m]*3+1]*3+2];

		PT3.x = d_vertice[d_face[d_leafNodeIndex[m]*3+2]*3];
		PT3.y = d_vertice[d_face[d_leafNodeIndex[m]*3+2]*3+1];
		PT3.z = d_vertice[d_face[d_leafNodeIndex[m]*3+2]*3+2];

		//get plane
		Point normal=getNormal(PT1,PT3,PT2);
		float d, denom, nom,t;
		Point intersectionP, intersection;
		d = 0-Dot(normal , PT1);
		//find intersection points
		denom=Dot(normal,ray.dir);
		if(fabs(denom)<EPSILON)
		{//parallel no intersection
			continue;
		}
		else
		{
			nom=Dot(normal, ray.pos)+d;
			t=-nom/denom;//distance
			if(t<=0)
			{//interseciton is on the back of the ray's start point 
				continue;
			}
			else
			{
				//whether in the triangle's plane
				//intersectionP=make_Point(ray.pos.x+ray.dir.x*t,ray.pos.y+ray.dir.y*t,ray.pos.z+ray.dir.z*t);

				intersectionP.x = ray.pos.x+ray.dir.x*t;
				intersectionP.y = ray.pos.y+ray.dir.y*t;
				intersectionP.z = ray.pos.z+ray.dir.z*t;

				if( t <minDis&&pointInTriangle(PT1,PT2,PT3,intersectionP))
				{
					minDis=t;	//min distance;
					minObjNum=m;
					planeNormal=normal;

					d_normal.x = normal.x ;
					d_normal.y = normal.y ;
					d_normal.z = normal.z ;

					d_intersection.x = intersectionP.x;
					d_intersection.y = intersectionP.y;
					d_intersection.z = intersectionP.z;



				}
			}
		}
	}
	return minObjNum;
}

////////////////////////////////////////////////////////////////////////////////
//intersectBox method
////////////////////////////////////////////////////////////////////////////////
__device__ bool intersectBox(Ray& r, Point& boxmin, Point& boxmax, float &tnear, float &tfar)
{
	//if(fabs(r.dir.x)<EPSILON || fabs(r.dir.y)<EPSILON || fabs(r.dir.z)<EPSILON)
	//	return false;
	// compute intersection of ray with all six bbox planes
	Point invR ;//= make_Point(1.0f) / r.dir;
	invR.x = 1.0/r.dir.x;
	invR.y = 1.0/r.dir.y;
	invR.z = 1.0/r.dir.z;
	Point tbot = multi(invR , minus(boxmin , r.pos));
	Point ttop = multi(invR , minus(boxmax , r.pos));

	// re-order intersections to find smallest and largest on each axis
	Point tmin = make(min(ttop.x, tbot.x),min(ttop.y, tbot.y),min(ttop.z, tbot.z));
	Point tmax = make(max(ttop.x, tbot.x),max(ttop.y, tbot.y),max(ttop.z, tbot.z));

	// find the largest tmin and the smallest tmax
	float largest_tmin = max(max(tmin.x, tmin.y), max(tmin.x, tmin.z));
	float smallest_tmax = min(min(tmax.x, tmax.y), min(tmax.x, tmax.z));

	tnear = largest_tmin;
	tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

__device__ int traverse(Ray& ray, KDInsideNode* d_insideNode, float* vertice, unsigned int* face, int* leafNodeIndex, KDLeafNode* leafNode, Point& normal, Point& intersection)
{

	int cur_index =0;
	float tnear,tfar;
	Point intersectionPoint;
	float intersectionvalue =0;

	Stack stack[200];
	int capacity = 0;
	int result;

	stack[capacity].index =0;
	stack[capacity].leafNode = false;
	capacity++;
	while(capacity>0)
	{
		 
		capacity--;

		while(stack[capacity].leafNode&&capacity>=0)
		{
			result = intersectLeaf(ray, stack[capacity].index, vertice, face, leafNode, leafNodeIndex,normal, intersection);
			if(result!=-1)
				return result;
			else
			{
				capacity--;
				continue;
			}
		}
		if(!stack[capacity].leafNode)
		{
			cur_index = stack[capacity].index;

			if(intersectBox(ray, d_insideNode[cur_index].aabb.minPoint, d_insideNode[cur_index].aabb.maxPoint, tnear, tfar))
			{
				intersectionPoint.x = ray.pos.x + tnear*ray.dir.x;
				intersectionPoint.y = ray.pos.y + tnear*ray.dir.y;
				intersectionPoint.z = ray.pos.z + tnear*ray.dir.z;

				switch(d_insideNode[cur_index].splitAxis)
				{
				case Axis_X:
					intersectionvalue = intersectionPoint.x;
					break;
				case Axis_Y:
					intersectionvalue = intersectionPoint.y;
					break;	
				case Axis_Z:
					intersectionvalue = intersectionPoint.z;
					break;
				}
				if(intersectionvalue < d_insideNode[cur_index].splitValue)
				{ // left part
					if(d_insideNode[cur_index].right!=-1)
					{
						if(d_insideNode[cur_index].RightLeaf)
							stack[capacity].leafNode = true;
						else
							stack[capacity].leafNode = false;

						stack[capacity].index = d_insideNode[cur_index].right;
						capacity++;
					}
					if(d_insideNode[cur_index].left!=-1)
					{
						if(d_insideNode[cur_index].LeftLeaf)
							stack[capacity].leafNode = true;
						else
							stack[capacity].leafNode = false;

						stack[capacity].index = d_insideNode[cur_index].left;
						capacity++;
					}
				}
				else
				{ // right part
					if(d_insideNode[cur_index].left!=-1)
					{
						if(d_insideNode[cur_index].LeftLeaf)
							stack[capacity].leafNode = true;
						else
							stack[capacity].leafNode = false;

						stack[capacity].index = d_insideNode[cur_index].left;
						capacity++;
					}
					if(d_insideNode[cur_index].right!=-1)
					{
						if(d_insideNode[cur_index].RightLeaf)
							stack[capacity].leafNode = true;
						else
							stack[capacity].leafNode = false;

						stack[capacity].index = d_insideNode[cur_index].right;
						capacity++;
					}

				}
			}
		}
	}
	if(capacity ==0)
		return -1;
}

__global__ void Render(Ray* rray, float* gpu_vertice, unsigned int* gpu_face, int* gpu_index, KDInsideNode* gpu_insideNode, KDLeafNode* gpu_leafNode,unsigned char* gpu_color)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	 
	
		 
	
	if(index<WIDTH*HEIGH)
	{
		Point diffuse;
		Point ambient;
		Point specular; 
		Point color;
		Point intersectionPoint;
		Point normal;
		float dif;
		float cos_data;
		Point reflect;
		Point rayFromEye;

		
		Ray ray;
		ray.pos.x = rray[index].pos.x;
		ray.pos.y = rray[index].pos.y;
		ray.pos.z = rray[index].pos.z;
		
		ray.dir.x = rray[index].dir.x;
		ray.dir.y = rray[index].dir.y;
		ray.dir.z = rray[index].dir.z;
 
		Light light;
		light.pos.x = 5;
		light.pos.y = 5;
		light.pos.z = 5;
		light.col.x =1;
		light.col.y =1;
		light.col.z =1;


		int kind = traverse(ray,gpu_insideNode, gpu_vertice, gpu_face, gpu_index,gpu_leafNode,normal,intersectionPoint);

		ambient.x = ambient.y = ambient.z =0.0;
		diffuse.x = diffuse.y = diffuse.z =0.0;
		specular.x = specular.y = specular.z =0.0;

		if(kind!=-1)
		{
			  
			ambient.x = AMBIENT*RGB_R;
			ambient.y = AMBIENT*RGB_G;
			ambient.z = AMBIENT*RGB_B;

			Point p;
			p.x = intersectionPoint.x;
			p.y = intersectionPoint.y;
			p.z = intersectionPoint.z;
			Normalize(p);

			dif  =Dot(p,normal);
			if(dif>0)
			{
				Ray temp;
				Point temp_point1,temp_point2;
				temp.pos = intersectionPoint;
				temp.dir = minus(light.pos , intersectionPoint);
				kind = traverse(temp,gpu_insideNode, gpu_vertice, gpu_face, gpu_index,gpu_leafNode,temp_point1,temp_point2);
				if(kind ==-1)
				{	
					color.x = ambient.x;
					color.y = ambient.y;
					color.z = ambient.z; 
					gpu_color[2+3*index] = (unsigned char)(color.x*255);
					gpu_color[1+3*index] = (unsigned char)(color.y*255);
					gpu_color[0+3*index] =   (unsigned char)(color.z*255);
					return;
				}
				else
				{
					diffuse.x = dif*DIFFUSE;
					diffuse.y = dif*DIFFUSE;
					diffuse.z = dif*DIFFUSE;
				}
			}
			else
			{
				color.x = ambient.x;
				color.y = ambient.y;
				color.z = ambient.z;
				gpu_color[2+3*index] = (unsigned char)(color.x*255);
				gpu_color[1+3*index] = (unsigned char)(color.y*255);
				gpu_color[0+3*index] =   (unsigned char)(color.z*255);
				return;
			}
			if(dif>0)
			{
				reflect.x = normal.x *dif*2-p.x;
				reflect.y = normal.y *dif*2-p.y;
				reflect.z = normal.z *dif*2-p.z;
				
				rayFromEye.x = ray.pos.x - intersectionPoint.x;
				rayFromEye.y = ray.pos.y - intersectionPoint.y;
				rayFromEye.z = ray.pos.z - intersectionPoint.z;
				Normalize(rayFromEye);

				cos_data = reflect.x*rayFromEye.x + reflect.y*rayFromEye.y + reflect.z*rayFromEye.z;
	 
				if(cos_data>0)
				{
					cos_data = pow(cos_data,PHON); 
					specular.x = light.col.x*cos_data*SPECULAR;
					specular.y = light.col.y*cos_data*SPECULAR;
					specular.z = light.col.z*cos_data*SPECULAR;
				}

			}
			
		}			
			color.x = diffuse.x  + ambient.x + specular.x;
			color.y = diffuse.y  + ambient.y + specular.y;
			color.z = diffuse.z  + ambient.z + specular.z;
			if(color.x >1.0)
				color.x =1.0;
			if(color.y >1.0)
				color.y =1.0;		
			if(color.z >1.0)
				color.z =1.0;

		 
		gpu_color[2+3*index] = (unsigned char)(color.x*255);
		gpu_color[1+3*index] = (unsigned char)(color.y*255);
		gpu_color[0+3*index] =   (unsigned char)(color.z*255);

		 
		

	}
	 
		
		 
	 
}

void processMouse(int button, int state, int x, int y) {

	if( button == GLUT_WHEEL_UP )
		scale -= 0.05;
	if( button == GLUT_WHEEL_DOWN )
		scale += 0.05;
	if ( state == GLUT_DOWN )
	{
		mouse_x = x; mouse_y = y;

		if( button == GLUT_LEFT_BUTTON )
			left_tag = true;
		if( button == GLUT_RIGHT_BUTTON )
			right_tag = true;
		//	cout << "left down!" << endl;
	}

	if ( state == GLUT_UP )
	{
		left_tag = false;
		right_tag = false;
		if( button == GLUT_LEFT_BUTTON )
		{
			anger_x += anger_move_x;
			anger_y += anger_move_y;
			anger_move_x = 0;
			anger_move_y = 0;
		}

		if( button == GLUT_RIGHT_BUTTON )
		{
			shift_x += shift_move_x;
			shift_y += shift_move_y;
			shift_move_x = 0;
			shift_move_y = 0;
		}
		//	cout << "left up!" << endl;
	}

}


void processMouseActiveMotion(int x, int y) {

	if ( left_tag )
	{
		anger_move_x = ( y - mouse_y ) * 1.0 / 800;
		anger_move_y = ( x - mouse_x ) * 1.0 / 800;
		//	cout << anger_x << endl;
	}
	if ( right_tag )
	{
		shift_move_x = x - mouse_x;
		shift_move_y = mouse_y - y;
	}
}

int   main(int   argc,   char   **argv) 
{ 

	if(mesh.loadFile("export/dolphins.obj"))
	{
		cout<<"successful";
	}
	else
		cout<<"failer";


	int* index = new int[mesh.m_numFace];
	for(int i=0;i<mesh.m_numFace;i++)
		index[i] =i;

	tree.ConstructKDTree(tree.root,mesh.m_numFace,index,mesh.face,mesh.vertice,mesh.aabb,0,true);
	tree.PrepareMemory();

//	kdnode = tree.root;





	///////////////////////////////////////////////////////////////////////////////////////////////////
	Ray* rays;
	rays  = new Ray[WIDTH *HEIGH];
	//pixelData = new  GLubyte[WIDTH *HEIGH*3];
	pixelData = (GLubyte*)malloc(WIDTH*HEIGH*3);
	if(pixelData == 0)
		return 0;

	//CUDA_Prepare(&mesh, &tree,pixelData);
	//unsigned char gpu_color[3];
	int pixel_y=0,pixel_x =0;//pixles
	float i_inc = 1.0/WIDTH;
	float j_inc = 1.0/HEIGH;
	float j=-0.5,i;
	int count_data =0;


	for(pixel_y=0; pixel_y < HEIGH; j+=j_inc, pixel_y++ )
	{
		for(pixel_x =0,i=-0.5*WIDTH/HEIGH;pixel_x < WIDTH; i+=i_inc*WIDTH/HEIGH, pixel_x++)
		{
		//	cout<<"pixel_y "<<pixel_y<<"  pixel_x "<<pixel_x<<endl;
			rays[count_data].dir.x= i; //+camera.look_at.x -camera.location.x;
			rays[count_data].dir.y= -1.5 ;//+camera.look_at.y -camera.location.y;
			rays[count_data].dir.z =-j;// camera->direction.z;//*(camera.look_at.z -camera.location.z);

			rays[count_data].pos.x = 0;//camera->location.x;
			rays[count_data].pos.y = 7;//camera->location.y;
			rays[count_data].pos.z = 0;//camera->location.z;

			/*
			Render(rays[count_data],mesh.vertice, mesh.face, tree.leafNodeIndex, tree.insideNode, tree.leafNode,gpu_color);

			pixelData[3*count_data] = gpu_color[0];
			pixelData[3*count_data+1] = gpu_color[1];
			pixelData[3*count_data+2] = gpu_color[2];
			*/
			count_data++;

		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////



	/*
	glutInit(&argc,   argv); 

	glutInitDisplayMode(GLUT_DEPTH | GLUT_RGBA | GLUT_DOUBLE); 

	glutInitWindowPosition(300, 30);
	glutInitWindowSize(WIDTH, HEIGH);
	glutCreateWindow( "kdtree"); 

	glViewport(0,   0,   800,  800); 
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 

	glOrtho(-200,   200,   -200,   200, 1, 10000); 
	// gluPerspective (0, 800 /600.0, 10.0, 100.0);
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity(); 


	gluLookAt(0,0,1000,0,0,0,0,1.0,0);

	glutDisplayFunc(draw); 
	glutIdleFunc(draw);

	glEnable(GL_DEPTH_TEST);

	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutKeyboardFunc(keyboard);

	glutMainLoop();
	*/


	Ray* gpu_rays;
	unsigned int* gpu_face;
	float* gpu_vertice;
	int* gpu_leafNodeIndex;
	KDLeafNode* gpu_leafNode;
	KDInsideNode* gpu_insideNode;
	GLubyte* gpu_color;

	cudaMalloc((void**)&gpu_rays, sizeof( Ray )*(WIDTH*HEIGH));
	cudaMalloc((void**)&gpu_face, sizeof( unsigned int )*(mesh.m_numFace)*3);
	cudaMalloc((void**)&gpu_vertice, sizeof( float )*(mesh.m_numVertice)*3);
	cudaMalloc((void**)&gpu_leafNodeIndex, sizeof( int )*(tree.SizeofLeafNode));
	cudaMalloc((void**)&gpu_leafNode, sizeof( KDLeafNode )*(tree.numLeafNode));
	cudaMalloc((void**)&gpu_insideNode, sizeof( KDInsideNode )*(tree.numInsideNode));
	cudaMalloc((void**)&gpu_color, sizeof( GLubyte )*(WIDTH*HEIGH*3));

	cudaMemcpy(gpu_rays, rays, sizeof(Ray)*WIDTH*HEIGH, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_face, mesh.face, sizeof( unsigned int )*(mesh.m_numFace)*3, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_vertice, mesh.vertice, sizeof( float )*(mesh.m_numVertice)*3, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_leafNodeIndex, tree.leafNodeIndex, sizeof( int )*(tree.SizeofLeafNode), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_leafNode, tree.leafNode, sizeof( KDLeafNode )*(tree.numLeafNode), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_insideNode, tree.insideNode, sizeof( KDInsideNode )*(tree.numInsideNode), cudaMemcpyHostToDevice);

	int blocks = WIDTH*HEIGH/512;
	if(WIDTH*HEIGH%512 !=0)
		blocks++;
	printf("\nblocks:%d\n",blocks);
	
	Render<<<blocks,512>>>(gpu_rays,gpu_vertice, gpu_face, gpu_leafNodeIndex, gpu_insideNode, gpu_leafNode,gpu_color);
	cutilDeviceSynchronize();
 
	cudaMemcpy(pixelData, gpu_color, sizeof( GLubyte )*(WIDTH*HEIGH*3), cudaMemcpyDeviceToHost);


	//for(int i=0;i<WIDTH*HEIGH*3;i++)
	//	pixelData[i]=120;

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(WIDTH, HEIGH);
	glutCreateWindow("RayTracing");
	glutDisplayFunc(&myDisplay);
	glutMainLoop();
	delete [] rays;
	delete [] index;
	free(pixelData);


	cudaFree(gpu_rays);
	cudaFree(gpu_vertice);
	cudaFree(gpu_face);
	cudaFree(gpu_leafNodeIndex);
	cudaFree(gpu_insideNode);
	cudaFree(gpu_leafNode);
	cudaFree(gpu_color);




	return   0;                           
} 

void draw( void )
{
	glPushMatrix();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glColor4f(1.0,1.0,1.0,1.0);

	glRotatef(10.0f, 1.0f, 0.0f, 0.0f);
	glRotatef(10.0f, 0.0f, 1.0f, 0.0f);

	glRotatef( ( anger_x + anger_move_x ) * 90, 1.0f, 0.0f, 0.0f);
	glRotatef( ( anger_y + anger_move_y ) * 90, 0.0f, 1.0f, 0.0f);

	glTranslatef( ( shift_x + shift_move_x ) * 0.5, 0.0f, 0.0f);
	glTranslatef( 0.0f, ( shift_y + shift_move_y ) * 0.5, 0.0f);

	glScalef( scale, scale, scale );

	drawGraph();


	glPopMatrix();
	glutSwapBuffers();

	//	glFlush();
}

void drawGraph()
{

	KDNode* cur;
	//KDLeafNode* cur;
	//cur = &(tree.leafNode[index_node]);
	cur = kdnode;
	drawBox(cur->aabb.maxPoint,cur->aabb.minPoint);

	glBegin(GL_LINES);
	//glBegin(GL_LINE_LOOP);
	for( int i=0;i<cur->numObject;i++)
	{
		int face = cur->object[i];
		//int face = tree.leafNodeIndex[(cur->begin)+i];
		glVertex3f(mesh.vertice[mesh.face[face*3]*3],mesh.vertice[mesh.face[face*3]*3+1],mesh.vertice[mesh.face[face*3]*3+2]);
		glVertex3f(mesh.vertice[mesh.face[face*3+1]*3],mesh.vertice[mesh.face[face*3+1]*3+1],mesh.vertice[mesh.face[face*3+1]*3+2]);
		glVertex3f(mesh.vertice[mesh.face[face*3+2]*3],mesh.vertice[mesh.face[face*3+2]*3+1],mesh.vertice[mesh.face[face*3+2]*3+2]);
	}
	glEnd();


}

void drawBox(Point& p1,Point& p2)
{
	// glColor3f(1.0f,0.0f,0.0f);
	glBegin(GL_LINES);
	glVertex3f(p1.x,p1.y,p1.z);
	glVertex3f(p1.x,p1.y,p2.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p1.x,p1.y,p1.z);
	glVertex3f(p1.x,p2.y,p1.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p1.x,p1.y,p1.z);
	glVertex3f(p2.x,p1.y,p1.z);
	glEnd();


	glBegin(GL_LINES);
	glVertex3f(p2.x,p2.y,p1.z);
	glVertex3f(p2.x,p2.y,p2.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p2.x,p1.y,p2.z);
	glVertex3f(p2.x,p2.y,p2.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p1.x,p2.y,p2.z);
	glVertex3f(p2.x,p2.y,p2.z);
	glEnd();

	///
	glBegin(GL_LINES);
	glVertex3f(p1.x,p1.y,p2.z);
	glVertex3f(p1.x,p2.y,p2.z);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(p1.x,p1.y,p2.z);
	glVertex3f(p2.x,p1.y,p2.z);
	glEnd();


	glBegin(GL_LINES);
	glVertex3f(p1.x,p2.y,p1.z);
	glVertex3f(p1.x,p2.y,p2.z);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(p1.x,p2.y,p1.z);
	glVertex3f(p2.x,p2.y,p1.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p2.x,p1.y,p1.z);
	glVertex3f(p2.x,p2.y,p1.z);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(p2.x,p1.y,p1.z);
	glVertex3f(p2.x,p1.y,p2.z);
	glEnd();

}

void keyboard(unsigned char key,int x, int y)
{
	switch(key)
	{
	case 's':
		scale +=0.5;
		break;
	case 'm':
		scale -=0.5;
		break;
	case 'k':
		index_node++;
		kdnode = kdnode->left;

	}
}