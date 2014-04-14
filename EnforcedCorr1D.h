#pragma once
#include <cmath>
#include <windows.h>		// Header File For Windows
#include <math.h>			// Math Library Header File
#include <stdio.h>			// Header File For Standard Input/Output


class EnforcedCorr1D
	/*
	Class for creating all the stimuli on the 200 Hz screen (can change to DLP)...
	*/
{


private:
	float Pwhite; // probability of intensity of 1
	int numberofpixels;
	int parity;
	int parityNUM;
	int parityvals[10][360]; // -1 <--> 0, and 1 <--> 1; saves 8 MAX!
	int pixvals[360];
	int pixvalsnew[360];
	int updatetype; // this will be a switch in the update to determine what kind of update rule to use
	int direction;
	bool circularBC; // circular boundary conditions?
	void copypixvalstopixvalsnew();
	void copypixvalsnewtopixvals();
	void update1(int x1t1,int x2t2); // functions to do basic correlation across 1 pixel groups
	void update1corr(int x1t1,int x2t2, float c);

	// OBSOLETE -- UNUSED code!
	void update2(int parity); // parity 0 or 1 of update
	bool test2(int v1,int v2,int position);
	// add here update2(); etc.
	// and update3(); etc.

	void updateparityvals();
	void updateparitypixvals(int P,int NUM);
	int computeparity(int NUM, int loc);

	void updateparitypixvals_arbitrarypattern(int P,int NUM, int (&toff)[10], int (&xoff)[10]);
	int computeparity_arbitrarypattern(int NUM, int loc, int (&toff)[10], int (&xoff)[10]);

	int select_frames;
public:
	
	//////////basic
	EnforcedCorr1D();  // **
	virtual ~EnforcedCorr1D(); // **

	////////////// FLAG FOR UPDATES
	bool PIXUPDATED;

	//////////////perform basic settings
	void setnumberofpixels(int numin);
	void setPwhite(float Pin);
	void setupdatetype(int updatein);
	void setBC(bool circBC);
	void setparity(int parityin);
	void setparityNUM(int parityNUMin);

	/////////performreal operations
	void updatepixels();
	void randomizepixels();
	int getpixelvalue(int pixnum);
	void setDir(int dir);
	void setBackFrames(int inframes);

} ;