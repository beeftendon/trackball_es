// Animation.h: interface for the CAnimation class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ANIMATION_H__E481D60A_F0B0_4FED_B47D_7CA923708947__INCLUDED_)
#define AFX_ANIMATION_H__E481D60A_F0B0_4FED_B47D_7CA923708947__INCLUDED_

//#if _MSC_VER > 1000
#pragma once
//#endif // _MSC_VER > 1000

#include <windows.h>		// Header File For Windows
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
#include <stdio.h>
#include <math.h>

class CAnimation  
{

public:

    float tau; // time scale in seconds
    float blobSize;
    int blobCount;
	float speed; // speed of motion, in weird units or deg/s

	float Zheight;
	float windowWidth; //=blobMaxX*blobSize;
    float windowHeight; //=blobMaxY*blobSize;

private:
	const float PI;
    //double blobMaxX; // = 160; //windowWidth/blobSize;
    //double blobMaxY; // = 120; //windowHeight/blobSize;

    struct SBlob 
    {
	    double x, y; //defines the top left right x,y coordinate for each blob
    };

    SBlob *blobArray, *pBlob;


    

    const int maxNumBlobs; // = 1000;


    float *pVertex;
    //double vertexArray[4000]; //8*maxNumBlobs];  up to 12 if 3 points per vertex
	float vertexArray[40000];  // 12*maxNumBlobs -- 4 vertices of 3 points each


	

public:
	CAnimation();
	virtual ~CAnimation();
    void CreateBlobs(void);
    void RenderScene(float dt);
	void RenderScene_rotation(float dt);
	void RenderScene_rotation_halves(float dt);
	void RenderScene_rotation_grad(float dt);

};

#endif // !defined(AFX_ANIMATION_H__E481D60A_F0B0_4FED_B47D_7CA923708947__INCLUDED_)
