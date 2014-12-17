#include "Node.h"


enum SplitAxis {Axis_X,Axis_Y,Axis_Z,LEAF};

struct KDInsideNode
{
	SplitAxis		splitAxis;
	float splitValue;
	int left;
	int right;
	bool LeftLeaf;
	bool RightLeaf;
	AABB aabb;
};

struct KDLeafNode
{
	int begin;
	int numObject;
	AABB aabb;
};

struct KDNode
{
	SplitAxis		splitAxis;
	int* object;
	int numObject;

	KDNode* left;
	KDNode* right;
	KDNode* parent;
	float splitValue;
	short floor;
	AABB aabb; 
};

class KDTree
{
public:
	KDTree();
	KDTree(AABB aabb);
	~KDTree();

	unsigned int m_floor;
	AABB bounding;
	KDNode* root;

	KDLeafNode* leafNode;
	KDInsideNode* insideNode;
	int* leafNodeIndex;

	void ConstructKDTree(int leafNodeSize, KDNode* parent,int num,int*index, unsigned int *face, float* vertice, AABB& aabb,int floor,bool isLeft);
	void SetIndex(int parentID, KDNode* current, int& InterSize, int& LeafSize, bool isLeft);
	void PrepareMemory();
	SplitAxis GetSplit(AABB &aabb,AABB& left, AABB& right,float& splitValue);
	void Split(int* index,unsigned int* face, float* vertice, int num, int* left,  int* right, int& numLeft, int& numRight, float& splitValue, SplitAxis axis);
	void Split_New(int* index,unsigned int* face, float* vertice, int num, int* left,int& numLeft, float& splitValue, SplitAxis axis, bool isLeft);

public:
	int numLeafNode;
	int numInsideNode;
	int SizeofLeafNode;

};
