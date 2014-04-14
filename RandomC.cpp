#include "stdafx.h"
#include "RandomC.h"

/////////////// functions to set and write stimulus stuff from fly measurements //////////////////////////////

RandomC::RandomC()
{

	tau = 0.2;
	dt = 1/240; // throw up at 240 Hz
	PI = 3.141592654; // pretty close
	x=0;
	last = 0;
	count = 0;

} ;

RandomC::~RandomC()
{
};

void RandomC::setTau(float tau_in)
{
	tau = tau_in;
};

void RandomC::setdt(float dt_in)
{
	dt = dt_in;
};

void RandomC::resetSeed(int seed_in)
{
	srand(seed_in);
};

float RandomC::Rrand()
{
	// flat distribution 0 to 1
	last = ((float)rand())/((float)RAND_MAX+1.0f);
	return last;
};

float RandomC::Rrandn()
{
	// gaussian, mean = 0, std = 1
	float rn1 = ((float)rand())/((float)RAND_MAX+1.0);            //between 0 and <1
    float rn2 = ((float)rand())/((float)RAND_MAX+1.0);
	last = sqrt( -2.0*log(1.0 - rn1) ) * cos( 2 * PI* rn2 );
    return last;
};

float RandomC::Rrandnt()
{
	// time series, updates at interval dt with correlation time tau
	
	float r = RandomC::Rrandn();            //between 0 and <1

	x += dt/tau*(r*sqrt(2.0f)/sqrt(dt/tau) - x);
	
	last = x;

	return x;
};

float RandomC::readLastValue()
{
	return last;
};

int RandomC::readCount()
{
	return count;
};

void RandomC::setCount(int count_in)
{
	count = count_in;
};