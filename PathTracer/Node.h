#pragma once




struct Point
{
	float x;
	float y;
	float z;
}; 

struct AABB
{
	Point minPoint;
	Point maxPoint;
};

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

struct Camera
{
	Point location;
	Point lookat;
	Point direction;
};