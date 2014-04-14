// Animation.cpp: implementation of the CAnimation class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Animation.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CAnimation::CAnimation() : 
    
    maxNumBlobs(5000),
	PI(3.14159)

{
    // set width and height
    //windowWidth=1.2;   //blobMaxX*blobSize;
    //windowHeight=1;  //blobMaxY*blobSize;

	//blobMaxX=30; //windowWidth/blobSize;
    //blobMaxY=20; //windowHeight/blobSize;
}

CAnimation::~CAnimation()  // delete any explicitly allocated memory ("new")
{
	delete blobArray;
}


void CAnimation::CreateBlobs(void) // this is called in trackstim declaration
{
    blobCount=0;
    blobArray = new SBlob[maxNumBlobs];

    for (blobCount=0; blobCount<maxNumBlobs; blobCount++) 
	{
        blobArray[blobCount].x=((rand()*windowWidth/RAND_MAX));
        blobArray[blobCount].y=((rand()*windowHeight/RAND_MAX));
    }

}


// Alternative rendering function called to draw next frame
// this rendering function uses a vertex array to plot all squares in one go
// should be more efficient but seems to be not much faster than RenderScene()
void CAnimation::RenderScene(float dt) 
{
    int i;
	double Zoff=Zheight;
	double Xoff=0; // to cener the damn thing
	double Yoff = 0;
	static float cumxsteps=0;
    
    pBlob = blobArray;  // array defined by pointer to first element; increment this to go through array
	float xstep = speed*dt;

	cumxsteps += xstep; // keep track of this to keep them on a grid


    for (i=0, pVertex = vertexArray; i<blobCount; i++, pBlob++) {
        pBlob->x += xstep;  //translate the blobs horizontally
        
        //according to coherence check whether we move this blob to a new random location
        if (rand () < dt / tau * RAND_MAX) 
		{
            //randomize blob postion
            pBlob->x = cumxsteps - blobSize*(int)(cumxsteps/blobSize) + blobSize*(int)((((float)rand()/(float)RAND_MAX))*windowWidth/blobSize); //+ xStepCumulative;
            pBlob->y = blobSize*(int)(((float)rand()*windowHeight/(float)RAND_MAX)/blobSize); //+ xStepCumulative;
        }
        else if (pBlob->x > windowWidth) 
		{ //else if scrolled off the edge, or ADD IN LESS THAN 0!
            //wrap around
            pBlob->x = pBlob->x - blobSize*(int)(windowWidth/blobSize); // put whole thing inside the mod! CHANGE THIS!
        }
		else if (pBlob->x < 0)
		{
			pBlob->x = pBlob->x + blobSize*(int)(windowWidth/blobSize);
		};

        //write new blob

        *pVertex++ = Xoff + pBlob->x;
        *pVertex++ = Yoff + pBlob->y;
		*pVertex++ = Zoff;

        *pVertex++ = Xoff + pBlob->x;
        *pVertex++ = Yoff + pBlob->y + blobSize;
		*pVertex++ = Zoff;

        *pVertex++ = Xoff + pBlob->x + blobSize;
        *pVertex++ = Yoff + pBlob->y + blobSize;
		*pVertex++ = Zoff;

        *pVertex++ = Xoff + pBlob->x + blobSize;
        *pVertex++ = Yoff + pBlob->y;
		*pVertex++ = Zoff;
    }
        
    
	glTranslatef(-windowWidth/2,-windowHeight/2,0.0f);
    
	glEnableClientState(GL_VERTEX_ARRAY);

    glVertexPointer(3, GL_FLOAT, 0, &vertexArray);  // do my stuff this way -- may be faster!
	
	glDrawArrays(GL_QUADS,0, blobCount*4);

    glDisableClientState(GL_VERTEX_ARRAY);

	// Flush drawing commands, move into graphics buffer from memory (?)
    glFlush();

    //terminate if we have run maxloops

}

// renders scene in polar coordinates
void CAnimation::RenderScene_rotation(float dt) 
{
    int i;

	static float cumxsteps=0;

	float xstep=speed*dt;  // x becomes theta
	float blobTheta=blobSize*PI/180; // blob extent in theta // make a parameter?
	float radius=30;
	float height, theta;
	float Xoff=0.0f;
	float Yoff=0.0f;
	float Zoff=0.0f;
	
	if (abs(speed) > .0001)
		cumxsteps += xstep; // keep track of this to keep them on a grid
	else
		cumxsteps = 0;
    
    pBlob = blobArray;  // array defined by pointer to first element; increment this to go through array


    for (i=0, pVertex = vertexArray; i<blobCount; i++, pBlob++) {
        //not the first frame so move each dot //DO THIS BY TIME!
		pBlob->x += xstep;
        
        //according to coherence check whether we move this blob to a new random location
        //if (rand() > coherence * RAND_MAX) {
		if (rand () < dt / tau * RAND_MAX) 
		{    // do it with dt/tau, not coherence
            //randomize blob postion
            //pBlob->x = ((rand()*2*PI/RAND_MAX)); //+ xStepCumulative;
            pBlob->x = cumxsteps - blobTheta*(int)(cumxsteps/blobTheta) + blobTheta*(int)((((float)rand()/(float)RAND_MAX))*2*PI/blobTheta); //create on grid...;
            //pBlob->x = blobTheta*(int)((((float)rand()/(float)RAND_MAX))*2*PI/blobTheta); //create on grid...;
            pBlob->y = (((float)rand()/(float)RAND_MAX)); //+ xStepCumulative;
        }
        /*if (pBlob->x > 2*PI) 
		{ //else if scrolled off the edge, or ADD IN LESS THAN 0!
            //wrap around
            pBlob->x = pBlob->x - 2*PI;
			//pBlob->y = ((rand()*blobMaxY/RAND_MAX))*blobSize;
		}
		if (pBlob->x < 0)
		{
			pBlob->x = 2*PI - pBlob->x;
			//pBlob->y = ((rand()*blobMaxY/RAND_MAX))*blobSize;            
        }*/ // WHY WAS I EVEN WORRIED ABOUT THIS? IT'S PERIODIC!

        //write new blob

		theta=pBlob->x;
		//height=pBlob->y*Zheight/windowHeight; // goes to new height
		height=radius*blobTheta*((int)(pBlob->y*(Zheight)/(radius*blobTheta))); // goes to new height


        *pVertex++ = Xoff + radius*cos(theta);  // fields of pointer
        *pVertex++ = Yoff + radius*sin(theta);
		*pVertex++ = Zoff + height;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta);
        *pVertex++ = Yoff + radius*sin(theta);
		*pVertex++ = Zoff + height + radius*blobTheta;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta);
		*pVertex++ = Zoff + height + radius*blobTheta;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta);
		*pVertex++ = Zoff + height;  // z value of 0
    }
        
    glEnableClientState(GL_VERTEX_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertexArray);  // do my stuff this way -- may be faster!

    glDrawArrays(GL_QUADS,0, blobCount*4);

    glDisableClientState(GL_VERTEX_ARRAY);

	//glLoadIdentity();
	
	// Flush drawing commands, move into graphics buffer from memory (?)
    glFlush();

};

void CAnimation::RenderScene_rotation_halves(float dt) 
{
    int i;

	float xstep=speed*dt;  // x becomes theta
	float blobTheta=5*PI/180; // blob extent in theta // make a parameter?
	float radius=10;
	float height, theta;
	float Xoff=0.0f;
	float Yoff=0.0f;
	float Zoff=-5.0f;
	bool top;
    
    pBlob = blobArray;  // array defined by pointer to first element; increment this to go through array


    for (i=0, pVertex = vertexArray; i<blobCount; i++, pBlob++) {
		top = (sin(pBlob->x) > 0); // it's on the top half
        //not the first frame so move each dot //DO THIS BY TIME!
		if (top)
			pBlob->x += xstep;
		else
			pBlob->x += - xstep;
        
        //according to coherence check whether we move this blob to a new random location
        //if (rand() > coherence * RAND_MAX) {
		if (rand () < dt / tau * RAND_MAX) 
		{    // do it with dt/tau, not coherence
            //randomize blob postion
            pBlob->x = ((rand()*2*PI/RAND_MAX)); //+ xStepCumulative;
            pBlob->y = ((rand()*windowHeight/RAND_MAX)); //+ xStepCumulative;
        }
        else if ((pBlob->x > PI) && top) 
		{ //else if scrolled off the edge, or ADD IN LESS THAN 0!
            //wrap around
            pBlob->x = pBlob->x - PI;
		}
		else if ((pBlob->x < 0) && top)
		{
			pBlob->x = PI + pBlob->x;
        }
		else if ((pBlob->x > 2*PI) && !top)
			pBlob->x = pBlob->x - PI;
		else if ((pBlob->x < PI) && !top)
			pBlob->x = pBlob->x + PI;


        //write new blob

		theta=pBlob->x;
		height=pBlob->y*Zheight/windowHeight; // goes to new height

        *pVertex++ = Xoff + radius*cos(theta - blobTheta/2);  // fields of pointer
        *pVertex++ = Yoff + radius*sin(theta - blobTheta/2);
		*pVertex++ = Zoff + height;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta - blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta - blobTheta/2);
		*pVertex++ = Zoff + height + blobSize;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta/2);
		*pVertex++ = Zoff + height + blobSize;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta/2);
		*pVertex++ = Zoff + height;  // z value of 0
    }
        
    glEnableClientState(GL_VERTEX_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertexArray);  // do my stuff this way -- may be faster!

    glDrawArrays(GL_QUADS,0, maxNumBlobs*4);

    glDisableClientState(GL_VERTEX_ARRAY);

	// Flush drawing commands, move into graphics buffer from memory (?)
    glFlush();

};

void CAnimation::RenderScene_rotation_grad(float dt) 
{
    int i;

	float xstep=speed*dt;  // x becomes theta
	float blobTheta=5*3.14/180; // blob extent in theta // make a parameter?
	float radius=10;
	float height, theta;
	float Xoff=0.0f;
	float Yoff=0.0f;
	float Zoff=-5.0f;
	bool top;
    
    pBlob = blobArray;  // array defined by pointer to first element; increment this to go through array


    for (i=0, pVertex = vertexArray; i<blobCount; i++, pBlob++) {
		top = (sin(pBlob->x) > 0); // it's on the top half
        //not the first frame so move each dot //DO THIS BY TIME!
		pBlob->x += xstep*sin(pBlob->x);
        
        //according to coherence check whether we move this blob to a new random location
		if (rand () < dt / tau * RAND_MAX) 
		{    // do it with dt/tau, not coherence
            //randomize blob postion
            pBlob->x = ((rand()*2*PI/RAND_MAX)); //+ xStepCumulative;
            pBlob->y = ((rand()*windowHeight/RAND_MAX)); //+ xStepCumulative;
        }
        else if ((pBlob->x > PI) && top) 
		{ //else if scrolled off the edge, or ADD IN LESS THAN 0!
            //wrap around
            pBlob->x = pBlob->x - PI;
		}
		else if ((pBlob->x < 0) && top)
		{
			pBlob->x = PI + pBlob->x;
        }
		else if ((pBlob->x > 2*PI) && !top)
			pBlob->x = pBlob->x - PI;
		else if ((pBlob->x < PI) && !top)
			pBlob->x = pBlob->x + PI;


        //write new blob


		theta=pBlob->x;
		height=pBlob->y*Zheight/windowHeight; // goes to new height

        *pVertex++ = Xoff + radius*cos(theta - blobTheta/2);  // fields of pointer
        *pVertex++ = Yoff + radius*sin(theta - blobTheta/2);
		*pVertex++ = Zoff + height;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta - blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta - blobTheta/2);
		*pVertex++ = Zoff + height + blobSize;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta/2);
		*pVertex++ = Zoff + height + blobSize;  // z value of 0

        *pVertex++ = Xoff + radius*cos(theta + blobTheta/2);
        *pVertex++ = Yoff + radius*sin(theta + blobTheta/2);
		*pVertex++ = Zoff + height;  // z value of 0
    }
        
    glEnableClientState(GL_VERTEX_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertexArray);  // do my stuff this way -- may be faster!

    glDrawArrays(GL_QUADS,0, maxNumBlobs*4);

    glDisableClientState(GL_VERTEX_ARRAY);

	//glLoadIdentity();
	
	// Flush drawing commands, move into graphics buffer from memory (?)
    glFlush();

};

