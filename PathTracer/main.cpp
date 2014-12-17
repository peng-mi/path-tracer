#include <windows.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <stdlib.h>
#include <stdio.h> 
#include <math.h>
#include <time.h>
 
#include "Node.h"
#include <iostream>
#include "KDTree.h"
#include "glm.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream> 

using namespace std;

#define MINNUM		1.0e-6
#define RAD(angle)  (angle*3.14159265359/180.0)
#define  INVRAD(rad) (rad * 180.0/3.14159265359) 
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
double g_theta = 0.0f;
double g_phi = 0.0f;
double g_rho =0.0f;
double scale = 50;
bool left_tag = false;
bool right_tag = false;
bool middle_tag = false;
int mouse_x = 0;
int mouse_y = 0;
double scale_move_x =0;
Point orig_cam;
Point orig_view;
Point  tmpX,tmpY, tmpZ;
void processMouseActiveMotion(int x, int y);
void processMouse(int button, int state, int x, int y);
void TranslateView(float xoffset, float yoffset);
void keyboard(unsigned char key,int x, int y);
void Rotation();
void ForWard(float ratio);

extern "C" void Clean_Data();
extern "C" void Prepare_Data(bool second,bool g_lightIn, int blockSize, int sampleSize, int width, int height, GLMmodel* glm_model, KDTree& tree);
extern "C" void RayTracing(float* matrix, GLubyte* pixelData);
extern "C" void RayCasting(float* matrix, GLubyte* pixelData);

void Render();
void GetMatrix(Point& pa, Point& pb, float* matrix);
void parseSetting(const char* filename);


KDTree tree;
Camera g_camera;
Light g_light;
int g_width, g_height;
float g_matrix[40];
GLubyte* pixelData;
KDNode* kdnode;
int* index;
GLMmodel* glm_model;

bool g_normalized = true;
bool g_raycasting = true;
int  g_leafNodeSize = 0;
char* g_filename;
int g_sample_size = 0;
int g_blockSize = 0;
bool g_light_in = false;
bool g_secondary = false;

void Point_Normal(Point& vec)
{
	float sum = sqrtf(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
	sum = 1.0f/sum;
	vec.x *= sum;
	vec.y *= sum;
	vec.z *= sum;
}

void Point_CrossProduct(Point& pa, Point& pb, Point& out)
{
	out.x = pa.y*pb.z - pa.z*pb.y;
	out.y = pa.z*pb.x - pa.x*pb.z;
	out.z = pa.x*pb.y - pa.y*pb.x;
}

float Point_Dot(Point& pa, Point& pb)
{
	return pa.x*pb.x + pa.y*pb.y + pa.z*pb.z;
}


void Init_Camera()
{
	orig_view.x = g_camera.lookat.x;
	orig_view.y = g_camera.lookat.y;
	orig_view.z = g_camera.lookat.z;

	g_rho = sqrtf(g_camera.location.x*g_camera.location.x + g_camera.location.y*g_camera.location.y + g_camera.location.z*g_camera.location.z);

	g_phi = INVRAD(acos(g_camera.location.y/g_rho));

	if(fabsf(g_camera.location.x) > 0.0000001)
		g_theta = INVRAD(atan(g_camera.location.z/g_camera.location.x));
	else
		g_theta = 0.0f;

	g_theta += 90.0f;

	if(g_phi >=0.0)
	{
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = 1.0f;
		g_camera.direction.z = 0.0f;
	}
	else
	{
		
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = -1.0f;
		g_camera.direction.x = 0.0f;
	}

	tmpY.x = g_camera.direction.x;
	tmpY.y = g_camera.direction.y;
	tmpY.z = g_camera.direction.z;

	tmpZ.x =(float)( sin(RAD(g_phi))*cos(RAD(g_theta)) );
	tmpZ.y =(float)( cos(RAD(g_phi)) );
	tmpZ.z =(float)( sin(RAD(g_phi))*sin(RAD(g_theta)) );

	Point_Normal(tmpZ);
	Point_CrossProduct(tmpY, tmpZ, tmpX);
	
	float tmpNorm = sqrtf(tmpX.x*tmpX.x + tmpX.y*tmpX.y + tmpX.z*tmpX.z);
	if(tmpNorm < MINNUM)
	{
		tmpX.x =1.0f; tmpX.y =0.0f; tmpX.z = 0.0f;
	}

	Point_CrossProduct(tmpZ, tmpX, tmpY);
	tmpNorm = sqrtf(tmpY.x*tmpY.x + tmpY.y*tmpY.y + tmpY.z*tmpY.z);
	if(tmpNorm < MINNUM)
	{
		tmpY.x = 0.0f; tmpY.y =1.0f; tmpY.z = 0.0f;
	}

	Point_Normal(tmpY);
	Point_Normal(tmpX);

	orig_cam.x = tmpX.x*g_camera.location.x + tmpY.x*g_camera.location.y + tmpZ.x*g_camera.location.z;
	orig_cam.y = tmpX.y*g_camera.location.x + tmpY.y*g_camera.location.y + tmpZ.y*g_camera.location.z;
	orig_cam.z = tmpX.z*g_camera.location.x + tmpY.z*g_camera.location.y + tmpZ.z*g_camera.location.z;
}


void ForWard(float ratio)
{
	orig_cam.z = orig_cam.z*ratio;

	g_camera.location.x = tmpX.x* orig_cam.x + tmpY.x*orig_cam.y + tmpZ.x*orig_cam.z;
	g_camera.location.y = tmpX.y* orig_cam.x + tmpY.y*orig_cam.y + tmpZ.y*orig_cam.z;
	g_camera.location.z = tmpX.z* orig_cam.x + tmpY.z*orig_cam.y + tmpZ.z*orig_cam.z;

	g_camera.lookat.x = tmpX.x*orig_view.x + tmpY.x*orig_view.y + tmpZ.x*orig_view.z;
	g_camera.lookat.y = tmpX.y*orig_view.x + tmpY.y*orig_view.y + tmpZ.y*orig_view.z;
	g_camera.lookat.z = tmpX.z*orig_view.x + tmpY.z*orig_view.y + tmpZ.z*orig_view.z;
}

void TranslateView(float xoffset, float yoffset)
{
	orig_cam.x += xoffset;
	orig_cam.y -= yoffset;

	orig_view.x += xoffset;
	orig_view.y -= yoffset;

	tmpZ.x =(float)( sin(RAD(g_phi))*cos(RAD(g_theta)) );
	tmpZ.y =(float)( cos(RAD(g_phi)) );
	tmpZ.z =(float)( sin(RAD(g_phi))*sin(RAD(g_theta)) );

	if(g_phi >=0.0)
	{
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = 1.0f;
		g_camera.direction.z = 0.0f;
	}
	else
	{
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = -1.0f;
		g_camera.direction.z = 0.0f;
	}

	Point_Normal(tmpZ);

	Point_CrossProduct(tmpY, tmpZ, tmpX);
	
	float tmpNorm = sqrtf(tmpX.x*tmpX.x + tmpX.y*tmpX.y + tmpX.z*tmpX.z);
	if(tmpNorm < MINNUM)
	{
		tmpX.x =1.0f; tmpX.y =0.0f; tmpX.z = 0.0f;
	}

	Point_CrossProduct(tmpZ, tmpX, tmpY);
	tmpNorm = sqrtf(tmpY.x*tmpY.x + tmpY.y*tmpY.y + tmpY.z*tmpY.z);
	if(tmpNorm < MINNUM)
	{
		tmpY.x = 0.0f; tmpY.y =1.0f; tmpY.z = 0.0f;
	}

	Point_Normal(tmpX);
	Point_Normal(tmpY);

	g_camera.location.x = tmpX.x*orig_cam.x + tmpY.x*orig_cam.y + tmpZ.x*orig_cam.z;
	g_camera.location.y = tmpX.y*orig_cam.x + tmpY.y*orig_cam.y + tmpZ.y*orig_cam.z;
	g_camera.location.z = tmpX.z*orig_cam.x + tmpY.z*orig_cam.y + tmpZ.z*orig_cam.z;

	g_camera.lookat.x = tmpX.x*orig_view.x + tmpY.x*orig_view.y + tmpZ.x*orig_view.z;
	g_camera.lookat.y = tmpX.y*orig_view.x + tmpY.y*orig_view.y + tmpZ.y*orig_view.z;
	g_camera.lookat.z = tmpX.z*orig_view.x + tmpY.z*orig_view.y + tmpZ.z*orig_view.z;
}

void Rotation()
{
	float xoffset, yoffset;
	xoffset = anger_move_x;
	yoffset = anger_move_y;
	
	g_theta += 0.1f*xoffset;
	if(g_theta >= 360.0)
		g_theta -= 360.0;
	else if(g_theta <=-360.0f)
		g_theta += 360.0f;

	g_phi += 0.1f*yoffset;
	if(g_phi > 180.0f)
		g_phi -=  360.0f;
	else if(g_phi <-180.0f)
		g_phi += 360.0f;

	if(g_phi >=0.0)
	{
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = 1.0f;
		g_camera.direction.z = 0.0f;
	}
	else
	{
		g_camera.direction.x = 0.0f;
		g_camera.direction.y = -1.0f;
		g_camera.direction.z = 0.0f;
	}

	tmpY.x = g_camera.direction.x;
	tmpY.y = g_camera.direction.y;
	tmpY.z = g_camera.direction.z;

	tmpZ.x =(float)( sin(RAD(g_phi))*cos(RAD(g_theta)) );
	tmpZ.y =(float)( cos(RAD(g_phi)) );
	tmpZ.z =(float)( sin(RAD(g_phi))*sin(RAD(g_theta)) );

	Point_Normal(tmpZ);
	Point_CrossProduct(tmpY, tmpZ, tmpX);
	
	float tmpNorm = sqrtf(tmpX.x*tmpX.x + tmpX.y*tmpX.y + tmpX.z*tmpX.z);
	if(tmpNorm < MINNUM)
	{
		tmpX.x =1.0f; tmpX.y =0.0f; tmpX.z = 0.0f;
	}

	Point_CrossProduct(tmpZ, tmpX, tmpY);
	tmpNorm = sqrtf(tmpY.x*tmpY.x + tmpY.y*tmpY.y + tmpY.z*tmpY.z);
	if(tmpNorm < MINNUM)
	{
		tmpY.x = 0.0f; tmpY.y =1.0f; tmpY.z = 0.0f;
	}

	Point_Normal(tmpX);
	Point_Normal(tmpY);

	g_camera.location.x = tmpX.x*orig_cam.x + tmpY.x*orig_cam.y + tmpZ.x*orig_cam.z;
	g_camera.location.y = tmpX.y*orig_cam.x + tmpY.y*orig_cam.y + tmpZ.y*orig_cam.z;
	g_camera.location.z = tmpX.z*orig_cam.x + tmpY.z*orig_cam.y + tmpZ.z*orig_cam.z;

	g_camera.lookat.x = tmpX.x*orig_view.x + tmpY.x*orig_view.y + tmpZ.x*orig_view.z;
	g_camera.lookat.y = tmpX.y*orig_view.x + tmpY.y*orig_view.y + tmpZ.y*orig_view.z;
	g_camera.lookat.z = tmpX.z*orig_view.x + tmpY.z*orig_view.y + tmpZ.z*orig_view.z;
}

void myDisplay(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(g_width,g_height,GL_BGR_EXT,GL_UNSIGNED_BYTE,pixelData);

	glutSwapBuffers();
}

void processMouse(int button, int state, int x, int y) {

	 
	if ( state == GLUT_DOWN )
	{
		mouse_x = x; mouse_y = y;

		if( button == GLUT_LEFT_BUTTON )
			left_tag = true;
		if( button == GLUT_RIGHT_BUTTON )
			right_tag = true;
		if( button == GLUT_MIDDLE_BUTTON)
			middle_tag = true;

	}

	if ( state == GLUT_UP )
	{
		left_tag = false;
		right_tag = false;
		middle_tag = false;
		if( button == GLUT_LEFT_BUTTON )
		{
			//anger_move_x = 
			/*anger_x += anger_move_x;
			anger_y += anger_move_y;
			anger_move_x = 0;
			anger_move_y = 0;*/
		}

		if( button == GLUT_RIGHT_BUTTON )
		{
			shift_x += shift_move_x;
			shift_y += shift_move_y;
			shift_move_x = 0;
			shift_move_y = 0;
		}
	}
}

void processMouseActiveMotion(int x, int y) {

	if ( left_tag )
	{
		anger_move_x = (x - mouse_x );
		anger_move_y = (y - mouse_y );

		mouse_x = x;
		mouse_y = y;
	}
	if ( right_tag )
	{
		shift_move_x = x - mouse_x;
		shift_move_y = mouse_y - y;
	}
	if(middle_tag)
	{
		if(x > mouse_x)
			ForWard(0.99);
		else
			ForWard(1.01);

		mouse_x = x;


	}
}

void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'p':
		printf("camera_location:%f, %f, %f\n",g_camera.location.x, g_camera.location.y, g_camera.location.z);
		printf("camera_lookat:%f, %f, %f\n",g_camera.lookat.x, g_camera.lookat.y, g_camera.lookat.z);
		printf("camera_direction:%f, %f, %f\n",g_camera.direction.x, g_camera.direction.y, g_camera.direction.z);
		break;
	}
}

void Render()
{
	
	TranslateView(shift_move_x/g_width, shift_move_y/g_height);
	Rotation();	

	Point tmp_a, tmp_b;
	tmp_a.x = g_camera.lookat.x - g_camera.location.x;
	tmp_a.y = g_camera.lookat.y - g_camera.location.y;
	tmp_a.z = g_camera.lookat.z - g_camera.location.z;
	tmp_b.x = 0;
	tmp_b.y = 0;
	tmp_b.z = -1;
	GetMatrix(tmp_a, tmp_b, g_matrix);
	
	g_matrix[9]  = g_camera.location.x;
	g_matrix[10] = g_camera.location.y;
	g_matrix[11] = g_camera.location.z;
	
	g_matrix[12] = g_camera.direction.x;
	g_matrix[13] = g_camera.direction.y;
	g_matrix[14] = g_camera.direction.z;

	g_matrix[15] = g_camera.lookat.x;
	g_matrix[16] = g_camera.lookat.y;
	g_matrix[17] = g_camera.lookat.z;

	g_matrix[18] = g_light.pos.x;
	g_matrix[19] = g_light.pos.y;
	g_matrix[20] = g_light.pos.z;

	g_matrix[21] = g_light.col.x;
	g_matrix[22] = g_light.col.y;
	g_matrix[23] = g_light.col.z;

	clock_t t;
	t = clock();
	
	if(g_raycasting)
		RayCasting( g_matrix,  pixelData);
	else
		RayTracing( g_matrix,  pixelData);

	t =  clock() -t;
	printf("Time:%f\n", (float)t/CLOCKS_PER_SEC);
}

void parseSetting(const char* filename)
{
	
	ifstream iFile;
	iFile.open(filename, ifstream::in);
	string buf;
	int id =0;

	while(iFile.good())
	{
		getline(iFile, buf, ':');
		if(buf.find("resolution") != string::npos)
		{
			getline(iFile, buf, ',');
			g_width = atoi(buf.c_str());
			getline(iFile, buf, '\n');
			g_height = atoi(buf.c_str());
		}
		else if(buf.find("camera_location")!= string::npos)
		{
			getline(iFile, buf, ',');
			g_camera.location.x = (float)atof(buf.c_str());
			getline(iFile, buf, ',');
			g_camera.location.y = (float)atof(buf.c_str());
			getline(iFile, buf, '\n');
			g_camera.location.z = (float)atof(buf.c_str());
		}
		else if(buf.find("camera_lookat")!= string::npos)
		{
			getline(iFile, buf, ',');
			g_camera.lookat.x = (float)atof(buf.c_str());
			getline(iFile, buf, ',');
			g_camera.lookat.y = (float)atof(buf.c_str());
			getline(iFile, buf, '\n');
			g_camera.lookat.z = (float)atof(buf.c_str());
		}
		else if(buf.find("camera_direction")!= string::npos)
		{
			getline(iFile, buf, ',');
			g_camera.direction.x = (float)atof(buf.c_str());
			getline(iFile, buf, ',');
			g_camera.direction.y = (float)atof(buf.c_str());
			getline(iFile, buf, '\n');
			g_camera.direction.z = (float)atof(buf.c_str());
		}
		else if(buf.find("light_pos")!= string::npos)
		{
			getline(iFile, buf, ',');
			g_light.pos.x = (float)atof(buf.c_str());
			getline(iFile, buf, ',');
			g_light.pos.y = (float)atof(buf.c_str());
			getline(iFile, buf, '\n');
			g_light.pos.z = (float)atof(buf.c_str());
		}
		else if(buf.find("light_color")!= string::npos)
		{
			getline(iFile, buf, ',');
			g_light.col.x = (float)atof(buf.c_str());
			getline(iFile, buf, ',');
			g_light.col.y = (float)atof(buf.c_str());
			getline(iFile, buf, '\n');
			g_light.col.z = (float)atof(buf.c_str());
		}
		else if(buf.find("normalized")!= string::npos)
		{
			getline(iFile, buf, '\n');
			int tmp = atoi(buf.c_str());
			if(tmp == 1)
				g_normalized = true;
			else
				g_normalized = false;
		}
		else if(buf.find("ray_cast")!= string::npos)
		{
			getline(iFile, buf, '\n');
			int tmp = atoi(buf.c_str());
			if(tmp == 1)
				g_raycasting = true;
			else
				g_raycasting = false;
		}
		else if(buf.find("file")!= string::npos)
		{
			getline(iFile, buf, '\n');
			g_filename =strdup(buf.c_str());
		}
		else if(buf.find("sample_size") != string::npos)
		{
			getline(iFile, buf, '\n');
			g_sample_size = atoi(buf.c_str());
		}
		else if(buf.find("leafNodeSize") != string::npos)
		{
			getline(iFile, buf, '\n');
			g_leafNodeSize = atoi(buf.c_str());
		}
		else if(buf.find("blockSize") != string::npos)
		{
			getline(iFile, buf, '\n');
			g_blockSize = atoi(buf.c_str());
		}
		else if(buf.find("light_in") != string::npos)
		{
			getline(iFile, buf, '\n');
			int tmp = atoi(buf.c_str());
			if(tmp == 1)
				g_light_in = true;
			else
				g_light_in = false;
			 
		}
		else if(buf.find("only_second") != string::npos)
		{
			getline(iFile, buf, '\n');
			int tmp = atoi(buf.c_str());
			if(tmp == 1)
				g_secondary = true;
			else
				g_secondary = false;
		}
	}
	Init_Camera();
}

void GetMatrix(Point& pa, Point& pb, float* matrix)
{
	float sum;
	float cos_h;

	sum = sqrt(pa.x*pa.x + pa.y*pa.y + pa.z*pa.z);
	sum = 1.0f/sum;
	pa.x *= sum;
	pa.y *= sum;
	pa.z *= sum;

	sum = sqrt(pb.x*pb.x + pb.y*pb.y + pb.z*pb.z);
	sum = 1.0f/sum;
	pb.x *= sum;
	pb.y *= sum;
	pb.z *= sum;

	cos_h = pa.x*pb.x + pa.y*pb.y + pa.z*pb.z;
	float sin_h  = sqrt(1-cos_h*cos_h)*(-1.0);

	Point temp;
	temp.x = temp.y = temp.z =0.0;
	temp.x = pa.y*pb.z - pa.z*pb.y;
	temp.y = pa.z*pb.x - pa.x*pb.z;
	temp.z = pa.x*pb.y - pa.y*pb.x;

	matrix[0] = temp.x*temp.x*(1.0-cos_h) + cos_h;
	matrix[1] = temp.x*temp.y*(1.0-cos_h) - sin_h*temp.z;
	matrix[2] = temp.x*temp.z*(1.0-cos_h) + sin_h*temp.y;

	matrix[3] = temp.x*temp.y*(1.0-cos_h) + sin_h*temp.z;
	matrix[4] = temp.y*temp.y*(1.0-cos_h) + cos_h;
	matrix[5] = temp.y*temp.z*(1.0-cos_h) - sin_h*temp.x;

	matrix[6] = temp.x*temp.z*(1.0-cos_h) - sin_h*temp.y;
	matrix[7] = temp.z*temp.y*(1.0-cos_h) + sin_h*temp.x;
	matrix[8] = temp.z*temp.z*(1.0-cos_h) + cos_h;
}

int   main(int   argc,   char   **argv) 
{ 
	

	parseSetting(argv[1]);
	string file_name(argv[1]);
	size_t found = file_name.find_last_of('/');
	if(found != string::npos)
		file_name = file_name.substr(0, found+1);
	file_name.append(g_filename);

	

	pixelData = (GLubyte*)malloc(g_width*g_height*3);
	if(pixelData == 0)
		return 0;

	float max_dimension;
	float dimension[3];
	AABB aabb;

	glm_model= glm_model->glmReadOBJ((char*)file_name.c_str());
	glm_model->glmDimensions(glm_model,dimension);
	
	index = new int[glm_model->numtriangles];
	for(int i=0;i<glm_model->numtriangles;i++)
		index[i] =i;

	if(g_normalized)
	{
		max_dimension = dimension[0];
		if(max_dimension < dimension[1])
			max_dimension = dimension[1];
		if(max_dimension < dimension[2])
			max_dimension = dimension[2];
		glm_model->glmScale(glm_model,1.0f/max_dimension);
	
		max_dimension = 1.0f/max_dimension;

	
		aabb.minPoint.x = glm_model->m_min[0]*max_dimension;
		aabb.minPoint.y = glm_model->m_min[1]*max_dimension;
		aabb.minPoint.z = glm_model->m_min[2]*max_dimension;
	
		aabb.maxPoint.x = glm_model->m_max[0]*max_dimension;
		aabb.maxPoint.y = glm_model->m_max[1]*max_dimension;
		aabb.maxPoint.z = glm_model->m_max[2]*max_dimension;
	}
	else
	{
		aabb.minPoint.x = glm_model->m_min[0];
		aabb.minPoint.y = glm_model->m_min[1];
		aabb.minPoint.z = glm_model->m_min[2];
	
		aabb.maxPoint.x = glm_model->m_max[0];
		aabb.maxPoint.y = glm_model->m_max[1];
		aabb.maxPoint.z = glm_model->m_max[2];
	}

	printf("number of triangels: %d\n",glm_model->numtriangles);

	tree.ConstructKDTree(g_leafNodeSize, tree.root,glm_model->numtriangles,index,glm_model->vindices,glm_model->vertices,aabb,0,true);
	tree.PrepareMemory();
	printf("depth:%d\n",tree.m_floor);


	Prepare_Data(g_secondary, g_light_in, g_blockSize, g_sample_size,g_width, g_height, glm_model,tree);
	
	

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(800, 500);
	glutInitWindowSize(g_width, g_height);
	glutCreateWindow("RayTracing");
	glutIdleFunc(Render);
	glutDisplayFunc(&myDisplay);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutMainLoop();


	//clean data
	atexit(Clean_Data);

	free(pixelData);
	free(tree.leafNode);
	free(tree.insideNode);
	free(tree.leafNodeIndex);
	delete [] index;
	glm_model->glmDelete(glm_model);
	
	return   0;                           
} 
