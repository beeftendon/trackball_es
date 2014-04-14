#pragma once
#include <cmath>
#include <windows.h>		// Header File For Windows
#include <math.h>			// Math Library Header File
#include <stdio.h>			// Header File For Standard Input/Output


class RandomC
	/*
	Class for creating all the stimuli on the 200 Hz screen (can change to DLP)...
	*/
{


private:
	float dt, tau, PI, x, last;
	int count;

public:
	
	//////////basic
	RandomC();  // **
	virtual ~RandomC(); // **

	//////////////set various things in the stimulus
	void setTau(float tau_in);
	void setdt(float dt_in);
	void resetSeed(int seed_in);
	float Rrand();
	float Rrandn();
	float Rrandnt();
	float readLastValue();
	int readCount();
	void setCount(int count_in);
	
} ;