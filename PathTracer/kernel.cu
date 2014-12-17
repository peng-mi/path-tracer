#include <windows.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <time.h>
//#include <cutil_inline.h>
//#include "datamanager.h"
#include <math.h>
//#include <cutil_math.h> 
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
 
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/remove.h>
#include <thrust/copy.h>

//#include <cuda.h>
//#include <curand_kernel.h>
#include "Node.h"
#include <iostream>
#include "KDTree.h"
#include "glm.h"
using namespace std;
#define PARAMETER_SIZE 40

#define INF 100000			
#define EPSILON 0.0000001
#define PI 3.14159265
#define MAX_FLOAT 10000000.0f


#define PHON 30
#define AMBIENT 0.2
#define SPECULAR 0.4
#define DIFFUSE 0.4
#define RGB_R 1.0
#define RGB_G 1.0
#define RGB_B 1.0
 


unsigned int* gpu_face;
float* gpu_vertice;
int* gpu_leafNodeIndex;
KDLeafNode* gpu_leafNode;
KDInsideNode* gpu_insideNode;
GLubyte* gpu_color;
int d_numFace;
int d_numVertex;
int d_width;
int d_height;
int d_sampleSize;
int d_blockSize;
float3* gpu_sec_pos;
float3* gpu_sec_posDouble;
float3* gpu_sec_normal;
float3* gpu_sec_normalDouble;
int* gpu_sec_index;
int* gpu_sec_indexDouble;
int d_numActiveRay;
GLubyte* gpu_Scolor;
float* gpu_parameter; 
float* gpu_rand;
bool d_lightIn;
bool d_secondary;


struct is_filled
{
	__host__ __device__
	bool operator() (const int x)
	{
		return x != -1;
	}
};
 
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
	float v=sqrtf(vector.x*vector.x+vector.y*vector.y+vector.z*vector.z);
	v  =1.0f/v;
	vector.x *= v; 
	vector.y *= v; 
	vector.z *= v;
}
__device__ Point Point_minus(Point& p1, Point& p2)
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

__device__ bool intersectionTriangle(Ray& ray, Point& p1, Point& p2, Point& p3, Point& normal, float& t)
{
	float det, inv_det, u, v;
	Point e1 = Point_minus(p2, p1);
	Point e2 = Point_minus(p3, p1);
	Point pvec = CrossProduct(ray.dir, e2);

	normal = CrossProduct(e1, e2);
	Normalize(normal);

	det = Dot(e1, pvec);
	if(det < EPSILON && det > -EPSILON)
		return false;
	inv_det = 1.0f/det;
	Point tvec = Point_minus(ray.pos, p1);

	u = Dot(tvec, pvec)*inv_det;
	if( u <0.0f || u >1.0f)
		return false;

	Point qvec = CrossProduct(tvec, e1);
	v = Dot(ray.dir, qvec)*inv_det;
	if( v<0.0f || (u+v)>1.0f)
		return false;

	t = Dot(e2, qvec)*inv_det;
	if( t > 0.000001)
		return true;
	
	return false;
}

__device__ 
bool intersectLeaf(Ray& ray, int leaf, float* d_vertice, unsigned int * d_face, KDLeafNode* d_leafNode, int* d_leafNodeIndex,Point& d_normal, float& minDis)
{
	if(leaf==-1)
		return false;
	int start = d_leafNode[leaf].begin;
	int end = d_leafNode[leaf].begin+d_leafNode[leaf].numObject;
	 
	Point normal;
	int  minObj = -1;
	float t;

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

		if(intersectionTriangle(ray,PT1,PT2,PT3,normal,t))
		{
			if(t < minDis)
			{
				minObj = m;
				minDis = t;
				d_normal.x = normal.x;
				d_normal.y = normal.y;
				d_normal.z = normal.z;
			}
		}
	}

	if(minObj == -1)
		return false;
	else
		return true;
}

//intersectBox method
__device__ bool intersectBox(Ray& r, Point& boxmin, Point& boxmax, float &tnear, float &tfar)
{
	Point invR ;//= make_Point(1.0f) / r.dir;
	invR.x = 1.0/r.dir.x;
	invR.y = 1.0/r.dir.y;
	invR.z = 1.0/r.dir.z;
	Point tbot = multi(invR , Point_minus(boxmin , r.pos));
	Point ttop = multi(invR , Point_minus(boxmax , r.pos));

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

__device__ 
bool traverse(Ray& ray, KDInsideNode* d_insideNode, float* vertice, unsigned int* face, int* leafNodeIndex, KDLeafNode* leafNode, Point& normal, float& minDis)
{

	int cur_index =0;
	float tnear,tfar;
	Point intersectionPoint;
	float intersectionvalue =0;

	Normalize(ray.dir);
	Stack stack[12];
	int capacity = 0;
	int result = -1;

	float distance = MAX_FLOAT;
	Point tmp_normal;
	stack[capacity].index =0;
	stack[capacity].leafNode = false;
	capacity++;
	while(capacity>0)
	{
		capacity--;
		while(stack[capacity].leafNode&&capacity>=0)//leaf node intersection test
		{
			if(intersectLeaf(ray, stack[capacity].index, vertice, face, leafNode, leafNodeIndex,tmp_normal, distance))
			{
				result = capacity;
				minDis = distance;
				normal.x = tmp_normal.x;
				normal.y = tmp_normal.y;
				normal.z = tmp_normal.z;
			}
			capacity--;
			if(capacity == 0)
			{
				if(result != -1)
					return true;
				else
					return false;
			}
			//continue;
		}

		if(!stack[capacity].leafNode)// interal node
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
	if(result == -1)
		return false;
	else
		return true;
}

 
 
 

__device__ void getRotationMatrix(Point& pa, Point& pb, float* matrix)
{
	float cos_h, sin_h;
	Normalize(pa);
	Normalize(pb);

	cos_h = Dot(pa, pb);
	sin_h = sqrtf(1-cos_h*cos_h)*(-1.0);

	Point tmp;
	tmp.x = pa.y*pb.z - pa.z*pb.y;
	tmp.y = pa.z*pb.x - pa.x*pb.z;
	tmp.z = pa.x*pb.y - pa.y*pb.x;

	matrix[0] = tmp.x*tmp.x*(1.0f-cos_h) +cos_h;
	matrix[1] = tmp.x*tmp.y*(1.0f-cos_h) +sin_h*tmp.z;
	matrix[2] = tmp.x*tmp.z*(1.0f-cos_h) +sin_h*tmp.y;

	matrix[3] = tmp.x*tmp.y*(1.0f-cos_h) +sin_h*tmp.z;
	matrix[4] = tmp.y*tmp.y*(1.0f-cos_h) +cos_h;
	matrix[5] = tmp.y*tmp.z*(1.0f-cos_h) -sin_h*tmp.x;

	matrix[6] = tmp.x*tmp.z*(1.0f-cos_h) +sin_h*tmp.y;
	matrix[7] = tmp.z*tmp.y*(1.0f-cos_h) +sin_h*tmp.x;
	matrix[8] = tmp.z*tmp.z*(1.0f-cos_h) +cos_h;
}

__device__ void getSampleRayDir(float* gpu_rand, Point& mainDir, float* rayDir,int ray_index, int numActiveRay, int width, int height, int sampleSize)
{
	//random value generation
	float random_1  = gpu_rand[ray_index];
	float random_2  = gpu_rand[ray_index + width*height*sampleSize];
	random_1 = random_1*2.0f*PI;
	random_2 = acosf(sqrtf(random_2));

	Point random_ray;
	random_ray.x = sinf(random_1)*cosf(random_2);
	random_ray.y = sinf(random_1)*sinf(random_2);
	random_ray.z = cosf(random_2);

	float matrix[9];
	Point normalize_dir;
	normalize_dir.x = 0.0f;
	normalize_dir.y = 0.0f;
	normalize_dir.z = 1.0f;

	getRotationMatrix(mainDir, normalize_dir, matrix);

	 
	rayDir[0] = random_ray.x*matrix[0] + random_ray.y*matrix[1] + random_ray.z*matrix[2];
	rayDir[1] = random_ray.x*matrix[3] + random_ray.y*matrix[4] + random_ray.z*matrix[5];
	rayDir[2] = random_ray.x*matrix[6] + random_ray.y*matrix[7] + random_ray.z*matrix[8];
}

__global__ 
void GenerateSecondImage(bool secondary, int numActiveRay,int* gpu_sec_index, unsigned char* gpu_color, unsigned char* s_color, int sampleSize)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index >= numActiveRay)
		return;
	
	int node_index = gpu_sec_index[index];

	uint3 tmp_color;
	tmp_color.x = 0;
	tmp_color.y = 0;
	tmp_color.z = 0;
	
	#pragma unroll 32
	for(int i = 0; i< sampleSize; i++)
	{
		tmp_color.x += s_color[3*(i + index*sampleSize) + 0];
		tmp_color.y += s_color[3*(i + index*sampleSize) + 1];
		tmp_color.z += s_color[3*(i + index*sampleSize) + 2];
	}
	if(secondary == true)
	{
		gpu_color[node_index*3+2] =   (unsigned char)(tmp_color.x/sampleSize);
		gpu_color[node_index*3+1] = (unsigned char)(tmp_color.y/sampleSize);
		gpu_color[node_index*3+0] = (unsigned char)(tmp_color.z/sampleSize);
	}
	else
	{
		gpu_color[node_index*3+2] +=   (unsigned char)(tmp_color.x/sampleSize);
		gpu_color[node_index*3+1] += (unsigned char)(tmp_color.y/sampleSize);
		gpu_color[node_index*3+0] += (unsigned char)(tmp_color.z/sampleSize);
	}
}

__global__ 
void SecondRender(float* gpu_rand, int numActiveRay, int sampleSize, int width, int height, float* parameter, float* gpu_vertice, unsigned int* gpu_face, int* gpu_index, KDInsideNode* gpu_insideNode, KDLeafNode* gpu_leafNode, float3* gpu_sec_node, float3* gpu_sec_normal, int* gpu_sec_index, unsigned char* s_color)
{
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index >= sampleSize*numActiveRay)
		return;

	s_color[3*index +0] = 0;
	s_color[3*index +1] = 0;
	s_color[3*index +2] = 0;

	int idx = index/sampleSize;

	Point normal;
	normal.x = gpu_sec_normal[idx].x;
	normal.y = gpu_sec_normal[idx].y;
	normal.z = gpu_sec_normal[idx].z;

	Ray ray;
	ray.pos.x = gpu_sec_node[idx].x;
	ray.pos.y = gpu_sec_node[idx].y;
	ray.pos.z = gpu_sec_node[idx].z;

	Point eye;
	eye.x = parameter[9] -ray.pos.x ;
	eye.y = parameter[10] -ray.pos.y;
	eye.z = parameter[11] -ray.pos.z;

	Normalize(eye);

	float cos_ray = Dot(eye, normal);
	if(cos_ray <0.0f)
		return;
	Point main_dir;
	main_dir.x = 2.0f*cos_ray*normal.x - eye.x;
	main_dir.y = 2.0f*cos_ray*normal.y - eye.y;
	main_dir.z = 2.0f*cos_ray*normal.z - eye.z;

	
	float dir[3]; 
	getSampleRayDir(gpu_rand, main_dir, dir ,index, numActiveRay, width, height,sampleSize);
	ray.dir.x = dir[0];
	ray.dir.y = dir[1];
	ray.dir.z = dir[2];
	Normalize(ray.dir);

	float min_dis = MAX_FLOAT;
	if(traverse(ray, gpu_insideNode, gpu_vertice, gpu_face,gpu_index, gpu_leafNode, normal, min_dis))
	{
		s_color[3*index +0] = 128;
		s_color[3*index +1] = 128;
		s_color[3*index +2] = 128;
	}
}


__global__ 
void RenderTracer(bool lightIn, int numFace, int numNode, int gpu_width, int gpu_height, float* parameter, float* gpu_vertice, unsigned int* gpu_face, int* gpu_index, KDInsideNode* gpu_insideNode, KDLeafNode* gpu_leafNode,unsigned char* gpu_color, float3* gpu_sec_pos, float3* gpu_sec_normal, int* gpu_sec_index)
{
	//parameter [0-8]: rotation matrix, [9-17]: camera_location, direction, lookat, [18-23] lightPos light col, may be more lights
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index<gpu_width*gpu_height)
	{
		gpu_sec_index[index] = -1;
		gpu_sec_pos[index].x = MAX_FLOAT;
		gpu_sec_pos[index].y = MAX_FLOAT;
		gpu_sec_pos[index].z = MAX_FLOAT;
		gpu_sec_normal[index].x = MAX_FLOAT;
		gpu_sec_normal[index].y = MAX_FLOAT;
		gpu_sec_normal[index].z = MAX_FLOAT;

		Point diffuse;
		Point ambient;
		Point specular; 
		Point color;
		Point intersectionPoint;
		Point normal;
		float dif;
		Light light;

		light.pos.x = parameter[18];
		light.pos.y = parameter[19];
		light.pos.z = parameter[20];
		light.col.x = parameter[21];
		light.col.y = parameter[22];
		light.col.z = parameter[23];
		 

		Ray ray;
		ray.pos.x = parameter[9];
		ray.pos.y = parameter[10];
		ray.pos.z = parameter[11];
 

		float id_x, id_y;
		unsigned int _width, _height;
		_height = index/gpu_width;
		_width = index%gpu_width;

		id_x = (float)_width/gpu_height - 0.5f*(float)gpu_width/gpu_height;
		id_y = (float)_height/gpu_height -0.5f;

		 
		ambient.x = id_x;
		ambient.y = -id_y;
		ambient.z = -1.0f;
		Normalize(ambient);
		

		ray.dir.x = ambient.x*parameter[0] + ambient.y*parameter[1] + ambient.z*parameter[2];
		ray.dir.y = ambient.x*parameter[3] + ambient.y*parameter[4] + ambient.z*parameter[5];
		ray.dir.z =	ambient.x*parameter[6] + ambient.y*parameter[7] + ambient.z*parameter[8];
		Normalize(ray.dir);

		//////////////////////////////////
		//traverse this tree get the intersection point and normal
		ambient.x = ambient.y = ambient.z =0.0;
		diffuse.x = diffuse.y = diffuse.z =0.0;
		specular.x = specular.y = specular.z =0.0;

		float min_dis = MAX_FLOAT;
		if(traverse(ray, gpu_insideNode, gpu_vertice, gpu_face,gpu_index, gpu_leafNode, normal, min_dis))
		{
			

			intersectionPoint.x = ray.pos.x + min_dis*ray.dir.x;
			intersectionPoint.y = ray.pos.y + min_dis*ray.dir.y;
			intersectionPoint.z = ray.pos.z + min_dis*ray.dir.z;

			gpu_sec_index[index] = index;
			gpu_sec_pos[index].x = intersectionPoint.x;
			gpu_sec_pos[index].y = intersectionPoint.y;
			gpu_sec_pos[index].z = intersectionPoint.z;
			gpu_sec_normal[index].x = normal.x;
			gpu_sec_normal[index].y = normal.y;
			gpu_sec_normal[index].z = normal.z;

		 
		 
			Point p;
			p.x =  light.pos.x -intersectionPoint.x;
			p.y =  light.pos.y -intersectionPoint.y;
			p.z =  light.pos.z -intersectionPoint.z;
			Normalize(p);

			dif  =Dot(p,normal);
			if(dif>0)
			{

				Ray rayEye;
				rayEye.pos.x = intersectionPoint.x;
				rayEye.pos.y = intersectionPoint.y;
				rayEye.pos.z = intersectionPoint.z;
				rayEye.dir.x = light.pos.x -intersectionPoint.x;
				rayEye.dir.y = light.pos.y -intersectionPoint.y;
				rayEye.dir.z = light.pos.z -intersectionPoint.z ;
				Normalize(rayEye.dir);	 

				if(lightIn != false)
				{
					diffuse.x = DIFFUSE*dif*RGB_R;
					diffuse.y = DIFFUSE*dif*RGB_G;
					diffuse.z = DIFFUSE*dif*RGB_B;  
				}
				else if(traverse(rayEye, gpu_insideNode, gpu_vertice, gpu_face,gpu_index, gpu_leafNode, normal, min_dis) == false)//shadows and occluused
				{		
					diffuse.x = DIFFUSE*dif*RGB_R;
					diffuse.y = DIFFUSE*dif*RGB_G;
					diffuse.z = DIFFUSE*dif*RGB_B;  
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
		gpu_color[2+3*index] =(unsigned char)(color.x*255);
		gpu_color[1+3*index] =(unsigned char)(color.y*255);
		gpu_color[0+3*index] =(unsigned char)(color.z*255);
	}
}


__global__ 
void RenderKDTree(bool lightIn, int numFace, int numNode, int gpu_width, int gpu_height, float* parameter, float* gpu_vertice, unsigned int* gpu_face, int* gpu_index, KDInsideNode* gpu_insideNode, KDLeafNode* gpu_leafNode,unsigned char* gpu_color, float3* gpu_sec_pos, float3* gpu_sec_normal, int* gpu_sec_index)
{
	//parameter [0-8]: rotation matrix, [9-17]: camera_location, direction, lookat, [18-23] lightPos light col, may be more lights
	int index = blockIdx.x*blockDim.x + threadIdx.x;
	if(index<gpu_width*gpu_height)
	{
		gpu_sec_index[index] = -1;
		gpu_sec_pos[index].x = MAX_FLOAT;
		gpu_sec_pos[index].y = MAX_FLOAT;
		gpu_sec_pos[index].z = MAX_FLOAT;
		gpu_sec_normal[index].x = MAX_FLOAT;
		gpu_sec_normal[index].y = MAX_FLOAT;
		gpu_sec_normal[index].z = MAX_FLOAT;

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

		Light light;


		light.pos.x = parameter[18];
		light.pos.y = parameter[19];
		light.pos.z = parameter[20];
		light.col.x = parameter[21];
		light.col.y = parameter[22];
		light.col.z = parameter[23];
		 

		Ray ray;
		ray.pos.x = parameter[9];
		ray.pos.y = parameter[10];
		ray.pos.z = parameter[11];
 

		float id_x, id_y;
		unsigned int _width, _height;
		_height = index/gpu_width;
		_width = index%gpu_width;

		id_x = (float)_width/gpu_height - 0.5f*(float)gpu_width/gpu_height;
		id_y = (float)_height/gpu_height -0.5f;

		 
		ambient.x = id_x;
		ambient.y = -id_y;
		ambient.z = -1.0f;
		Normalize(ambient);
		

		ray.dir.x = ambient.x*parameter[0] + ambient.y*parameter[1] + ambient.z*parameter[2];
		ray.dir.y = ambient.x*parameter[3] + ambient.y*parameter[4] + ambient.z*parameter[5];
		ray.dir.z =	ambient.x*parameter[6] + ambient.y*parameter[7] + ambient.z*parameter[8];
		Normalize(ray.dir);

	

		//////////////////////////////////
		//traverse this tree get the intersection point and normal
		ambient.x = ambient.y = ambient.z =0.0;
		diffuse.x = diffuse.y = diffuse.z =0.0;
		specular.x = specular.y = specular.z =0.0;

		float min_dis = MAX_FLOAT;
		if(traverse(ray, gpu_insideNode, gpu_vertice, gpu_face,gpu_index, gpu_leafNode, normal, min_dis))
		{
			

			intersectionPoint.x = ray.pos.x + min_dis*ray.dir.x;
			intersectionPoint.y = ray.pos.y + min_dis*ray.dir.y;
			intersectionPoint.z = ray.pos.z + min_dis*ray.dir.z;

			gpu_sec_index[index] = index;
			gpu_sec_pos[index].x = intersectionPoint.x;
			gpu_sec_pos[index].y = intersectionPoint.y;
			gpu_sec_pos[index].z = intersectionPoint.z;
			gpu_sec_normal[index].x = normal.x;
			gpu_sec_normal[index].y = normal.y;
			gpu_sec_normal[index].z = normal.z;

		 
			ambient.x = AMBIENT*RGB_R;
			ambient.y = AMBIENT*RGB_G;
			ambient.z = AMBIENT*RGB_B;
	
			Point p;
			p.x =  light.pos.x -intersectionPoint.x;
			p.y =  light.pos.y -intersectionPoint.y;
			p.z =  light.pos.z -intersectionPoint.z;
			Normalize(p);

			dif  =Dot(p,normal);
			if(dif>0)
			{

				Ray rayEye;
				rayEye.pos.x = intersectionPoint.x;
				rayEye.pos.y = intersectionPoint.y;
				rayEye.pos.z = intersectionPoint.z;
				rayEye.dir.x = light.pos.x -intersectionPoint.x;
				rayEye.dir.y = light.pos.y -intersectionPoint.y;
				rayEye.dir.z = light.pos.z -intersectionPoint.z ;
				Normalize(rayEye.dir);

				//specular light calculation
				reflect.x = normal.x*2*dif-p.x;
				reflect.y = normal.y*2*dif-p.y;
				reflect.z = normal.z*2*dif-p.z;
				Normalize(reflect);

				rayFromEye.x = ray.pos.x - intersectionPoint.x;
				rayFromEye.y = ray.pos.y - intersectionPoint.y;
				rayFromEye.z = ray.pos.z - intersectionPoint.z;
				Normalize(rayFromEye);
				cos_data = Dot(reflect, rayFromEye);

				if(lightIn != false)//shadows and occluused
				{
					diffuse.x = DIFFUSE*dif*RGB_R;
					diffuse.y = DIFFUSE*dif*RGB_G;
					diffuse.z = DIFFUSE*dif*RGB_B;

					if(cos_data>0)
					{ 
						cos_data = powf(cos_data,PHON); 
						specular.x = light.col.x*cos_data*SPECULAR;
						specular.y = light.col.y*cos_data*SPECULAR;
						specular.z = light.col.z*cos_data*SPECULAR;
					}	  
				}
				else if(traverse(rayEye, gpu_insideNode, gpu_vertice, gpu_face,gpu_index, gpu_leafNode, normal, min_dis) == false)
				{
					diffuse.x = DIFFUSE*dif*RGB_R;
					diffuse.y = DIFFUSE*dif*RGB_G;
					diffuse.z = DIFFUSE*dif*RGB_B;

					if(cos_data>0)
					{ 
						cos_data = powf(cos_data,PHON); 
						specular.x = light.col.x*cos_data*SPECULAR;
						specular.y = light.col.y*cos_data*SPECULAR;
						specular.z = light.col.z*cos_data*SPECULAR;
					}	

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
		gpu_color[2+3*index] =(unsigned char)(color.x*255);
		gpu_color[1+3*index] =(unsigned char)(color.y*255);
		gpu_color[0+3*index] =(unsigned char)(color.z*255);
	}
}

extern "C"
void RayCasting(float* matrix, GLubyte* pixelData)
{
	cudaMemcpy(gpu_parameter, matrix, sizeof(float)*PARAMETER_SIZE, cudaMemcpyHostToDevice);

	int blockSize = d_blockSize;
	int blocks = d_width*d_height/blockSize;
	if( (d_width*d_height)%blockSize !=0)
		blocks++;
	RenderKDTree<<<blocks,blockSize>>>(d_lightIn, d_numFace, d_numVertex, d_width, d_height, gpu_parameter,gpu_vertice, gpu_face, gpu_leafNodeIndex, gpu_insideNode, gpu_leafNode,gpu_color, gpu_sec_pos, gpu_sec_normal, gpu_sec_index);
	//cutilDeviceSynchronize();
	cudaDeviceSynchronize();
	cudaMemcpy(pixelData, gpu_color, sizeof( GLubyte )*(d_width*d_height*3), cudaMemcpyDeviceToHost);
}

extern "C"
void RayTracing(float* matrix, GLubyte* pixelData)
{
	cudaMemcpy(gpu_parameter, matrix, sizeof(float)*PARAMETER_SIZE, cudaMemcpyHostToDevice);

	int blockSize = d_blockSize;
	int blocks = d_width*d_height/blockSize;
	if( (d_width*d_height)%blockSize !=0)
		blocks++;
	RenderTracer<<<blocks,blockSize>>>(d_lightIn, d_numFace, d_numVertex, d_width, d_height, gpu_parameter,gpu_vertice, gpu_face, gpu_leafNodeIndex, gpu_insideNode, gpu_leafNode,gpu_color, gpu_sec_pos, gpu_sec_normal, gpu_sec_index);
	//cutilDeviceSynchronize();
	cudaDeviceSynchronize();
	//compact the secondary ray
	thrust::device_ptr<int> index_ptr(gpu_sec_index);
	thrust::device_ptr<int> indexDouble_ptr(gpu_sec_indexDouble);
	cudaMemcpy(gpu_sec_indexDouble, gpu_sec_index,sizeof(int)*d_width*d_height,cudaMemcpyDeviceToDevice);
	thrust::device_ptr<float3> pos_ptr(gpu_sec_pos);
	thrust::device_ptr<float3> posDou_ptr(gpu_sec_posDouble);
	thrust::device_ptr<float3> normal_ptr(gpu_sec_normal);
	thrust::device_ptr<float3> normalDou_ptr(gpu_sec_normalDouble);
	thrust::copy_if(pos_ptr, pos_ptr+d_width*d_height,index_ptr, posDou_ptr, is_filled()); 
	thrust::copy_if(normal_ptr, normal_ptr+d_width*d_height,index_ptr, normalDou_ptr, is_filled()); 
	
	d_numActiveRay = thrust::count(index_ptr, index_ptr+d_width*d_height, -1);
	d_numActiveRay = d_width*d_height - d_numActiveRay; 
	printf("result:%d\n",d_numActiveRay);
	thrust::remove_copy(indexDouble_ptr, indexDouble_ptr+d_width*d_height, index_ptr, -1);

	blocks = (d_numActiveRay*d_sampleSize)/blockSize;
	if( (d_numActiveRay*d_sampleSize)%blockSize != 0)
		blocks++;
	printf("blocks:%d, threads:%d\n", blocks, d_numActiveRay*d_sampleSize);

	cudaDeviceSynchronize();
	//cutilDeviceSynchronize();
	SecondRender<<<blocks, blockSize>>>(gpu_rand, d_numActiveRay, d_sampleSize, d_width, d_height ,gpu_parameter, gpu_vertice, gpu_face, gpu_leafNodeIndex, gpu_insideNode, gpu_leafNode, gpu_sec_posDouble, gpu_sec_normalDouble, gpu_sec_index, gpu_Scolor);
	//cutilDeviceSynchronize();
	cudaDeviceSynchronize();

	blocks = (d_numActiveRay)/blockSize;
	if( (d_numActiveRay)%blockSize != 0)
		blocks++;
	GenerateSecondImage<<<blocks, blockSize>>>(d_secondary, d_numActiveRay,gpu_sec_index, gpu_color, gpu_Scolor,d_sampleSize);
	
	cudaDeviceSynchronize();
	//cutilDeviceSynchronize();
	cudaMemcpy(pixelData, gpu_color, sizeof( GLubyte )*(d_width*d_height*3), cudaMemcpyDeviceToHost);
	//cutilDeviceSynchronize();
	cudaDeviceSynchronize();
}
 
extern "C"
void Prepare_Data(bool secondary, bool lightIn, int blockSize, int sampleSize, int width, int height, GLMmodel* glm_model, KDTree& tree)
{
	d_numFace = glm_model->numtriangles;
	d_numVertex = glm_model->numvertices;
	d_width = width;
	d_height = height;
	d_sampleSize = sampleSize;
	d_blockSize = blockSize;
	d_lightIn = lightIn;
	d_secondary = secondary;

	cudaMalloc((void**)&gpu_rand, sizeof(float)*2*d_width*d_height*d_sampleSize);
	cudaMalloc((void**)&gpu_Scolor, sizeof(GLubyte)*3*d_width*d_height*d_sampleSize);
	cudaMalloc((void**)&gpu_parameter, sizeof(float)*PARAMETER_SIZE);
	cudaMalloc((void**)&gpu_color, sizeof( GLubyte )*(d_width*d_height*3));
	cudaMalloc((void**)&gpu_face, sizeof( unsigned int )*(glm_model->numtriangles)*3);
	cudaMalloc((void**)&gpu_vertice, sizeof( float )*(glm_model->numvertices+1)*3);
	cudaMalloc((void**)&gpu_leafNode, sizeof(KDLeafNode)*(tree.numLeafNode));
	cudaMalloc((void**)&gpu_insideNode, sizeof(KDInsideNode)*(tree.numInsideNode));
	cudaMalloc((void**)&gpu_leafNodeIndex, sizeof(int)*(tree.SizeofLeafNode));

	cudaMalloc((void**)&gpu_sec_pos, sizeof(float3)*(d_width*d_height));
	cudaMalloc((void**)&gpu_sec_posDouble, sizeof(float3)*(d_width*d_height));
	cudaMalloc((void**)&gpu_sec_normal, sizeof(float3)*(d_width*d_height));
	cudaMalloc((void**)&gpu_sec_normalDouble, sizeof(float3)*(d_width*d_height));
	cudaMalloc((void**)&gpu_sec_index, sizeof(int)*(d_width*d_height));
	cudaMalloc((void**)&gpu_sec_indexDouble, sizeof(int)*(d_width*d_height));
	
	cudaMemcpy(gpu_face, glm_model->vindices, sizeof( unsigned int )*(glm_model->numtriangles)*3, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_vertice, glm_model->vertices, sizeof( float )*(glm_model->numvertices+1)*3, cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_leafNode, tree.leafNode, sizeof(KDLeafNode)*(tree.numLeafNode), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_insideNode, tree.insideNode, sizeof(KDInsideNode)*(tree.numInsideNode), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_leafNodeIndex, tree.leafNodeIndex, sizeof(int)*(tree.SizeofLeafNode), cudaMemcpyHostToDevice);

	float* tmp_random;
	tmp_random = (float*)malloc(sizeof(float)*2*d_width*d_height*d_sampleSize);
	for(int i = 0; i<d_width*d_height*d_sampleSize*2; i++)
		tmp_random[i] = ((float)rand()/RAND_MAX);

	cudaMemcpy(gpu_rand, tmp_random, sizeof(float)*2*d_width*d_height*d_sampleSize, cudaMemcpyHostToDevice);
	free(tmp_random);
}

extern "C"
void Clean_Data()
{
	cudaFree(gpu_parameter);
	cudaFree(gpu_vertice);
	cudaFree(gpu_face);
	cudaFree(gpu_leafNodeIndex);
	cudaFree(gpu_insideNode);
	cudaFree(gpu_leafNode);
	cudaFree(gpu_color);

	cudaFree(gpu_sec_pos);
	cudaFree(gpu_sec_posDouble);
	cudaFree(gpu_sec_normal);
	cudaFree(gpu_sec_normalDouble);
	cudaFree(gpu_sec_index);
	cudaFree(gpu_sec_indexDouble);
	cudaFree(gpu_Scolor);
	cudaFree(gpu_rand);
}


