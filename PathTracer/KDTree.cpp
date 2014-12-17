#include "KDTree.h"
#include<stdio.h>
#include<iostream>
#include<fstream>
#include <stdlib.h>

 
//delete this line if you want to evenly divide the space
#define MEDIAN

int compare (const void * a, const void * b)
{
	int index =-1;
	if((*(float*)a - *(float*)b) >= 0.0)
		index =1;
	return index;
}

KDTree::KDTree(void)
{
	root =0;
	numLeafNode =0;
	numInsideNode =0;
	m_floor =0;
	SizeofLeafNode =0;
}

KDTree::KDTree(AABB aabb)
{
	root =0; 
	bounding.maxPoint.x = aabb.maxPoint.x;
	bounding.maxPoint.y = aabb.maxPoint.y;
	bounding.maxPoint.z = aabb.maxPoint.z;

	bounding.minPoint.x = aabb.minPoint.x;
	bounding.minPoint.y = aabb.minPoint.y;
	bounding.minPoint.z = aabb.minPoint.z;

	numLeafNode =0;
	numInsideNode =0;
	m_floor =0;
	SizeofLeafNode =0;
}
KDTree::~KDTree(void)
{
	delete [] leafNode;
	delete [] insideNode;
	delete [] leafNodeIndex;
}

/*
object		prepare GPU memory, transfer KDTree to array
*/
void KDTree::PrepareMemory()
{ 

	leafNode = new KDLeafNode[numLeafNode];
	insideNode = new KDInsideNode[numInsideNode];

	KDNode* cur = root;

	int indexNode =0;
	int inter_size =1;
	int leaf_size =0;
	leafNodeIndex = new int[SizeofLeafNode];
	SizeofLeafNode =0;

	insideNode[0].splitValue = root->splitValue;
	insideNode[0].splitAxis = root->splitAxis;

	insideNode[0].aabb.maxPoint.x=root->aabb.maxPoint.x;
	insideNode[0].aabb.maxPoint.y=root->aabb.maxPoint.y;
	insideNode[0].aabb.maxPoint.z=root->aabb.maxPoint.z;
	insideNode[0].aabb.minPoint.x=root->aabb.minPoint.x;
	insideNode[0].aabb.minPoint.y=root->aabb.minPoint.y;
	insideNode[0].aabb.minPoint.z=root->aabb.minPoint.z;
	insideNode[0].LeftLeaf = false;
	insideNode[0].RightLeaf = false;
	SetIndex(indexNode,cur->left, inter_size,leaf_size,true);
	SetIndex(indexNode,cur->right, inter_size,leaf_size,false);
}

/*
object		set index for each item of array to help it to find its children nodes.
parameter	parentID	parent ID
parameter	current		current node
parameter	InterSize	the internode's positon in InterNode array
parameter	LeafSize	the leafnode's position in LeafNode array
parameter	isLeaf		the indicate the current node is the left or the right child of its parent
*/
void KDTree::SetIndex(int parentID, KDNode* current, int& InterSize, int& LeafSize, bool isLeft)
{
	if(current)
	{
		int id = InterSize;
		if(current->splitAxis != LEAF)
		{

			insideNode[id].splitAxis = current->splitAxis;
			insideNode[id].splitValue = current->splitValue;

			insideNode[id].aabb.maxPoint.x = current->aabb.maxPoint.x;
			insideNode[id].aabb.maxPoint.y = current->aabb.maxPoint.y;
			insideNode[id].aabb.maxPoint.z = current->aabb.maxPoint.z;
			insideNode[id].aabb.minPoint.x = current->aabb.minPoint.x;
			insideNode[id].aabb.minPoint.y = current->aabb.minPoint.y;
			insideNode[id].aabb.minPoint.z = current->aabb.minPoint.z;

			if(isLeft)
			{
				insideNode[parentID].left = id;
				insideNode[parentID].LeftLeaf =  false;
			}
			else
			{
				insideNode[parentID].right = id;
				insideNode[parentID].RightLeaf =  false;
			}

			InterSize++;

			SetIndex(id,current->left, InterSize,LeafSize,true);
			SetIndex(id,current->right, InterSize,LeafSize,false);
		}
		else
		{

			if(isLeft)
			{
				insideNode[parentID].LeftLeaf =  true; 
				insideNode[parentID].left = LeafSize;//index
			}
			else
			{
				insideNode[parentID].RightLeaf =  true; 
				insideNode[parentID].right = LeafSize;
			}
			leafNode[LeafSize].begin = SizeofLeafNode;
			leafNode[LeafSize].numObject = current->numObject;
			LeafSize++;
			for(int i=0;i<current->numObject;i++)
			{
				leafNodeIndex[SizeofLeafNode] = current->object[i];
				SizeofLeafNode++;
			}
			delete [] current->object;
			current->object =0;
		}
	}
	else
	{
		if(isLeft)
			insideNode[parentID].left = -1;
		else
			insideNode[parentID].right = -1;
	}
}

/*
object						  construct a kdtree; this method is recursive 
parameter	parent			  the parent node
parameter	num 			  the length of index 
parameter	index			  a array to store the index of face we want to split
parameter	face			  face array
parameter	vertice			  vertice array
parameter	aabb			  the AABB of the sub geometry
parameter   floor			  the floor for each node
parameter	isLeft			  indicate the current node is the left or the right of its parent
*/

void KDTree::ConstructKDTree(int leafNodeSize, KDNode* parent,int num, int* index, unsigned int *face, float* vertice, AABB& aabb,int floor, bool isLeft)
{
	if(root ==0)
	{

		root = new KDNode();
		root->numObject = num;//face num
		root->object = index;//int[] 1,2,3,4...
		root->floor = floor;//0
		//assign mesh aabb to the root node
		root->aabb.maxPoint.x = aabb.maxPoint.x;
		root->aabb.maxPoint.y = aabb.maxPoint.y;
		root->aabb.maxPoint.z = aabb.maxPoint.z;
		root->aabb.minPoint.x = aabb.minPoint.x;
		root->aabb.minPoint.y = aabb.minPoint.y;
		root->aabb.minPoint.z = aabb.minPoint.z;

		AABB leftAABB, rightAABB;
		int* leftIndex = new int[num];
		int* rightIndex =new int[num];
		root->splitAxis = GetSplit(aabb,leftAABB,rightAABB, root->splitValue); 

		int numLeft, numRight;

		Split(index,face,vertice,num,leftIndex,rightIndex,numLeft,numRight,root->splitValue,root->splitAxis);
		//Split_New(index,face,vertice,num,leftIndex,numLeft,root->splitValue,root->splitAxis,true);

#ifdef MEDIAN
		switch(root->splitAxis)
		{
		case Axis_X:
			leftAABB.maxPoint.x = root->splitValue;
			rightAABB.minPoint.x = root->splitValue;
			break;
		case Axis_Y:
			leftAABB.maxPoint.y = root->splitValue;
			rightAABB.minPoint.y = root->splitValue;
			break;
		case Axis_Z:
			leftAABB.maxPoint.z = root->splitValue;
			rightAABB.minPoint.z = root->splitValue;
			break;
		}
#endif
		m_floor++;
		ConstructKDTree(leafNodeSize, root,numLeft, leftIndex, face,vertice, leftAABB,floor+1,true);
		delete [] leftIndex;
		//free(leftIndex) ;

		//int* rightIndex = (int*)malloc(sizeof(int)*num);
		//Split_New(index,face,vertice,num,rightIndex,numRight,root->splitValue,root->splitAxis,false);
		ConstructKDTree(leafNodeSize, root,numRight, rightIndex,face, vertice, rightAABB,floor+1,false);	
		numInsideNode++;
		//free(rightIndex) ;
		delete [] rightIndex;
	}
	else
	{
		if(num<=leafNodeSize)
		{//the node is small enough to do ray tracing
			if(num == 0)
			{
				if(isLeft)
					parent->left =0;
				else
					parent->right =0;
				return;
			}
			else
			{
				KDNode* node = new KDNode();
				if(isLeft)
					parent->left = node;
				else
					parent->right = node;
				node->aabb.maxPoint.x = aabb.maxPoint.x;
				node->aabb.maxPoint.y = aabb.maxPoint.y;
				node->aabb.maxPoint.z = aabb.maxPoint.z;
				node->aabb.minPoint.x = aabb.minPoint.x;
				node->aabb.minPoint.y = aabb.minPoint.y;
				node->aabb.minPoint.z = aabb.minPoint.z;

				node->splitAxis = LEAF;
				node->numObject = num;
				node->floor = floor;
				node->parent = parent;
				//parent->left = parent->right = NULL;
				node->object = index;
				numLeafNode++;
				SizeofLeafNode +=num;
				node->left =0;
				node->right=0;
				return;
			}
		}
		else
		{//the need still need to be splited
			KDNode* node = new KDNode();
			if(isLeft)
				parent->left = node;
			else
				parent->right = node;
			node->numObject = num;
			node->object = index;
			node->floor = floor;
			node->parent = parent;

			node->aabb.maxPoint.x = aabb.maxPoint.x;
			node->aabb.maxPoint.y = aabb.maxPoint.y;
			node->aabb.maxPoint.z = aabb.maxPoint.z;
			node->aabb.minPoint.x = aabb.minPoint.x;
			node->aabb.minPoint.y = aabb.minPoint.y;
			node->aabb.minPoint.z = aabb.minPoint.z;

			AABB leftAABB, rightAABB;

			//int* leftIndex = (int*)malloc(sizeof(int)*num);
			int* leftIndex = new int[num];
			int* rightIndex =new int[num];
			
			node->splitAxis = GetSplit(aabb,leftAABB,rightAABB, node->splitValue);

			int numLeft, numRight;
			//Split_New(index,face,vertice,num,leftIndex,numLeft,node->splitValue,node->splitAxis, true);
			Split(index,face,vertice,num,leftIndex,rightIndex,numLeft,numRight,node->splitValue,node->splitAxis);

#ifdef MEDIAN
			switch(node->splitAxis)
			{
			case Axis_X:
				leftAABB.maxPoint.x = node->splitValue;
				rightAABB.minPoint.x = node->splitValue;
				break;
			case Axis_Y:
				leftAABB.maxPoint.y = node->splitValue;
				rightAABB.minPoint.y = node->splitValue;
				break;
			case Axis_Z:
				leftAABB.maxPoint.z = node->splitValue;
				rightAABB.minPoint.z = node->splitValue;
				break;
			}
#endif
			if(m_floor < floor+1)
				m_floor = floor+1;
			
			ConstructKDTree(leafNodeSize, node,numLeft, leftIndex,face, vertice, leftAABB,floor+1,true);
			if(numLeft>leafNodeSize)
			//	free(leftIndex) ;
				delete [] leftIndex ;

			//int* rightIndex = (int*)malloc(sizeof(int)*num);
			//Split_New(index,face,vertice,num,rightIndex,numRight,node->splitValue,node->splitAxis, false);
			ConstructKDTree(leafNodeSize, node,numRight, rightIndex,face, vertice, rightAABB,floor+1,false);

			numInsideNode++;
			if(numRight>leafNodeSize)
			//	free(rightIndex);
				delete [] rightIndex ;
		}
	}

}


void KDTree::Split_New(int* index,unsigned int* face, float* vertice, 
				   int num, int* left, int& numLeft, float& splitValue, SplitAxis axis, bool isLeft)
{
	numLeft = 0;
	int offset;

	switch(axis)
	{
	case Axis_X:
		offset =0;
		break;
	case Axis_Y:
		offset =1;
		break;
	case Axis_Z:
		offset =2;
		break;
	default:
		offset =-1;
	}

#ifdef MEDIAN
	float* temp;
	temp = new float[num*3];
	for(int i=0;i<num;i++)
	{	
		temp[3*i+0] = vertice[face[index[i]*3]*3 + offset];
		temp[3*i+1] = vertice[face[index[i]*3+1]*3 + offset];
		temp[3*i+2] = vertice[face[index[i]*3+2]*3 + offset];
	}

	qsort (temp, num*3, sizeof(float), compare);

	splitValue = temp[num*3/2];
#endif
	for(int i=0;i<num;i++)
	{
		if(isLeft)
		{
			if ( (vertice[face[index[i]*3]*3 + offset]<=splitValue) 
			|| (vertice[face[index[i]*3+1]*3 + offset]<=splitValue) 
			|| (vertice[face[index[i]*3+2]*3 + offset]<=splitValue) )
			left[numLeft++]=index[i];

		}
		else
		{
			if ( (vertice[face[index[i]*3]*3 + offset]>=splitValue) 
			|| (vertice[face[index[i]*3+1]*3 + offset]>=splitValue) 
			|| (vertice[face[index[i]*3+2]*3 + offset]>=splitValue) )
			left[numLeft++]=index[i];
		}
	} 
#ifdef MEDIAN
	delete [] temp;
#endif
}


///////////////////////////////////////////////////////////////////////////////////
////parameters:		in: index face vertice num splitValue, axis
////				out: left, right, numLeft, numRight, 
/*
object						   Split geometry objects
parameter	index			   a array to indicate which face we need to split 
parameter	face			   face information
parameter	num				   the length of index
parameter	left			   a array to store the left part of split plane
parameter	right			   a array to store the right part of split plane
parameter	numLeft			   the length of left array
parameter	numRight		   the length of right array
parameter	splitValue		   split plane's value
parameter	axis			   the split plane's axis
*/
void KDTree::Split(int* index,unsigned int* face, float* vertice, 
				   int num, int* left,  int* right, int& numLeft, int& numRight, float& splitValue, SplitAxis axis)
{
	numLeft = 0;
	numRight = 0;
	int offset;

	switch(axis)
	{
	case Axis_X:
		offset =0;
		break;
	case Axis_Y:
		offset =1;
		break;
	case Axis_Z:
		offset =2;
		break;
	default:
		offset =-1;
	}

#ifdef MEDIAN
	float* temp;
	temp = new float[num*3];
	for(int i=0;i<num;i++)
	{	
		temp[3*i+0] = vertice[face[index[i]*3]*3 + offset];
		temp[3*i+1] = vertice[face[index[i]*3+1]*3 + offset];
		temp[3*i+2] = vertice[face[index[i]*3+2]*3 + offset];
	}

	qsort (temp, num*3, sizeof(float), compare);

	splitValue = temp[num*3/2];
#endif

	for(int i=0;i<num;i++)
	{
		if ( (vertice[face[index[i]*3]*3 + offset]<=splitValue) 
			&& (vertice[face[index[i]*3+1]*3 + offset]<=splitValue) 
			&& (vertice[face[index[i]*3+2]*3 + offset]<=splitValue) )
			left[numLeft++]=index[i];
		else if ( (vertice[face[index[i]*3]*3 + offset]>=splitValue) 
			&& (vertice[face[index[i]*3+1]*3 + offset]>=splitValue) 
			&& (vertice[face[index[i]*3+2]*3 + offset]>=splitValue) )
			right[numRight++]=index[i];
		else
		{
			left[numLeft++]=index[i]; 
			right[numRight++]=index[i];
		}
	} 
#ifdef MEDIAN
	delete [] temp;
#endif
}


/////////////////////////////////////////////////////////////////////////////
///get spilitaxis for current aabb and split the aabb into left and right////
///parameter:	in: aabb
///				out: left, right splitValue
/*
object							Get the split plane
parameter	the bounding box   the size of father node
parameter	splitValue		   get the splited value
return		SplitAxis		   the split plane
*/
SplitAxis KDTree::GetSplit(AABB &aabb,AABB& left, AABB& right,float& splitValue)
{
	float rangeX = aabb.maxPoint.x - aabb.minPoint.x;
	float rangeY = aabb.maxPoint.y - aabb.minPoint.y;
	float rangeZ = aabb.maxPoint.z - aabb.minPoint.z;

	if(rangeX>=rangeY)
	{
		if(rangeX>=rangeZ)
		{
			splitValue = (aabb.maxPoint.x + aabb.minPoint.x)*0.5;

			left.minPoint.x = aabb.minPoint.x;
			left.minPoint.y = aabb.minPoint.y;
			left.minPoint.z = aabb.minPoint.z;

			left.maxPoint.x = splitValue;
			left.maxPoint.y = aabb.maxPoint.y;
			left.maxPoint.z = aabb.maxPoint.z;

			right.minPoint.x = splitValue;
			right.minPoint.y = aabb.minPoint.y;
			right.minPoint.z = aabb.minPoint.z;

			right.maxPoint.x = aabb.maxPoint.x;
			right.maxPoint.y = aabb.maxPoint.y;
			right.maxPoint.z = aabb.maxPoint.z;
			return Axis_X;
		}
		else
		{
			splitValue = (aabb.maxPoint.z + aabb.minPoint.z)*0.5;

			left.minPoint.x = aabb.minPoint.x;
			left.minPoint.y = aabb.minPoint.y;
			left.minPoint.z = aabb.minPoint.z;

			left.maxPoint.x = aabb.maxPoint.x;
			left.maxPoint.y = aabb.maxPoint.y;
			left.maxPoint.z = splitValue;

			right.minPoint.x = aabb.minPoint.x;
			right.minPoint.y = aabb.minPoint.y;
			right.minPoint.z = splitValue;

			right.maxPoint.x = aabb.maxPoint.x;
			right.maxPoint.y = aabb.maxPoint.y;
			right.maxPoint.z = aabb.maxPoint.z;
			return Axis_Z;
		}
	}
	else
	{
		if(rangeY>=rangeZ)
		{
			splitValue = (aabb.maxPoint.y + aabb.minPoint.y)*0.5;

			left.minPoint.x = aabb.minPoint.x;
			left.minPoint.y = aabb.minPoint.y;
			left.minPoint.z = aabb.minPoint.z;

			left.maxPoint.x = aabb.maxPoint.x;
			left.maxPoint.y = splitValue ;
			left.maxPoint.z = aabb.maxPoint.z;

			right.minPoint.x = aabb.minPoint.x;
			right.minPoint.y = splitValue ;
			right.minPoint.z = aabb.minPoint.z;

			right.maxPoint.x = aabb.maxPoint.x;
			right.maxPoint.y = aabb.maxPoint.y;
			right.maxPoint.z = aabb.maxPoint.z;

			return Axis_Y;
		}
		else
		{
			splitValue = (aabb.maxPoint.z + aabb.minPoint.z)*0.5;

			left.minPoint.x = aabb.minPoint.x;
			left.minPoint.y = aabb.minPoint.y;
			left.minPoint.z = aabb.minPoint.z;

			left.maxPoint.x = aabb.maxPoint.x;
			left.maxPoint.y = aabb.maxPoint.y;
			left.maxPoint.z = splitValue;

			right.minPoint.x = aabb.minPoint.x;
			right.minPoint.y = aabb.minPoint.y;
			right.minPoint.z = splitValue;

			right.maxPoint.x = aabb.maxPoint.x;
			right.maxPoint.y = aabb.maxPoint.y;
			right.maxPoint.z = aabb.maxPoint.z;
			return Axis_Z;
		}
	}
}

