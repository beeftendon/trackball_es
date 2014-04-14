#include "stdafx.h"
#include "TrackStim.h"
#include <math.h>

/////////////// functions to set and write stimulus stuff from fly measurements //////////////////////////////

TrackStim::TrackStim()
{

	PI=3.14159f;
	ScreenUpdateRate=120; // throw up at 200 Hz
	

	monitor_position_file = NULL;   // file to write to for positions
	monitor_stimulus_file = NULL;   // file to write for stimulus info
	parameter_file = NULL;          // parameters out, for a check
	error_file = NULL; 
	paramfilename[0] = NULL;

	TrackStim::setZeroTime();
	
	// initialize FMeas
	Fmeas.x = Fmeas.y = Fmeas.theta = Fmeas.xold = Fmeas.yold = Fmeas.thetaold = 0.0f;
	Fmeas.vx = Fmeas.vy = Fmeas.vtheta = 0.0f;
	Fmeas.tcurr = Fmeas.told = Fmeas.dt = 0.0f;
	Fmeas.dx = Fmeas.dy = Fmeas.dtheta = 0.0f;

	Fcurr=Fmeas;
	Scurr=Fcurr;
	Sfull=Scurr;
	Saccumulate=Sfull;

	// initialize Scurr, Sfull
	//Scurr.x = Scurr.y = Scurr.theta = 0.0f;    // all we keep track of here
	//Sfull.x = Sfull.y = Sfull.theta = 0.0f;    // and here

	// initialize stimulus stuff in case not read in
	numepochs=1;
	epochchoose=0;
	TrackStim::Stimulus=TrackStim::initializeStimulus();

	TrackStim::ViewPorts=TrackStim::initializeViewPositions();
	for (int ii=0; ii<500; ii++)
	{
		Map.x[ii] = 0; Map.y[ii] = 0;
	};
	Map.currpos = 0;
	
	// set up random dot stuff...
	Rdots.tau = Stimulus.tau;
	Rdots.blobSize = Stimulus.spacing;  // in degrees
	Rdots.windowHeight = Stimulus.arenasize*3;
	Rdots.windowWidth = Stimulus.arenasize*3;
	Rdots.Zheight = Stimulus.arenaheight;
	Rdots.blobCount = (int)(400*Stimulus.density);  // assign based on density
	Rdots.CreateBlobs();

	GammaCorrect = 2.2f;
	DistanceAccumulated = 0;
	FwdAccumulated = 0;
	framenumber = 0;

	VIEWCUBECORNER = FALSE;
	USE_FRUSTUM = TRUE; // no real need to do this yet -- looks same with and without this shift, in cylindrical world...

} ;

TrackStim::~TrackStim()
{
	// close files in use
	if (!(monitor_stimulus_file == NULL))
		fclose(monitor_stimulus_file); 
	if (!(monitor_position_file == NULL)) 
		fclose(monitor_position_file); 
	if (!(parameter_file == NULL))
		fclose(parameter_file);
	if (!(error_file == NULL))
		fclose(error_file);

};

TrackStim::STIMPARAMS TrackStim::initializeStimulus()
{
	TrackStim::STIMPARAMS Stim;   // declare it here...
	Stim.stimtype = 0;
	Stim.lum = .5;
	Stim.contrast = 1;
	Stim.duration = 1000;   // in seconds
	Stim.transgain = 1;
	Stim.rotgain = 1;
	Stim.trans.mean = 0;
	Stim.trans.amp = 0;
	Stim.trans.per = 1;
	Stim.trans.phase = 0;
	Stim.rot.mean = 0;
	Stim.rot.amp = 0;
	Stim.rot.per = 1; 
	Stim.rot.phase = 0;
	Stim.stimtrans.mean = 360;
	Stim.stimtrans.amp = 0;
	Stim.stimtrans.per = 1;
	Stim.stimtrans.phase = 0;
	Stim.stimrot.mean = 0;
	Stim.stimrot.amp = 0;
	Stim.stimrot.per = 1; 
	Stim.stimrot.phase = 30;
	Stim.spacing = 30; //.010; //.1 or 10
	Stim.density = .5f;
	Stim.tau = .1f;  // in seconds
	Stim.arenasize = 100; // in mm
	Stim.arenaheight = 30; // in mm
	Stim.tau2 = 0.3f; // in seconds...
	Stim.spacing2 = 0;
	Stim.right.rotgain = 1;
	Stim.right.transgain = 1;
	Stim.left.rotgain = 1;
	Stim.left.transgain = 1;
	Stim.changemethod = 0;

	Stim.aperture.location = 0;
	Stim.aperture.width = 0;
	
	return Stim;   // pass it back.
}

TrackStim::VIEWPOSITIONS TrackStim::initializeViewPositions()
{
	TrackStim::VIEWPOSITIONS temp;
	temp.x[0] = 0; temp.x[1] = 250; temp.x[2] = 500; temp.x[3] = 500;
	temp.y[0] = 0; temp.y[1] = 300; temp.y[2] = 0; temp.y[3] = 300;
	temp.w[0] = 200; temp.w[1] = 200; temp.w[2] = 200; temp.w[3] = 200;
	temp.h[0] = 200; temp.h[1] = 200; temp.h[2] = 200; temp.h[3] = 200;

	return temp;
};
void TrackStim::incrementFrameNumber()
{
	framenumber++;
};

void TrackStim::setZeroTime()
{
	QueryPerformanceCounter(&CounterBeginTime);  // time == 0
	QueryPerformanceFrequency(&CounterFreq);     // frequency of counter, needed to calculate time
};

float TrackStim::queryCurrTime()
{
	LARGE_INTEGER TimeNow;
	QueryPerformanceCounter(&TimeNow);
	return (float)(TimeNow.QuadPart - CounterBeginTime.QuadPart)/(float)CounterFreq.QuadPart;
};


float TrackStim::sineWave(float period)  // freq in Hz
{
	//DON'T USE THIS! PUT IN A SCREEN PRESENTATION TIME, FROM FCURR, FOR INSTANCE

	//float tc = TrackStim::queryCurrTime() - epochtimezero;
	float tc = Fcurr.tcurr - epochtimezero; 
	if (period > 0)   // ignore all if period is 0
		return sin(2*PI*tc/period);
	else
		return 0;
};

void TrackStim::setFlyMeasuredPos(float inx, float iny, float intheta, float intime)
{	

	Fmeas.told = Fmeas.tcurr;
	Fmeas.tcurr = intime; // set time of the newest measurements -- could add as argument, i think.
	Fmeas.dt = Fmeas.tcurr - Fmeas.told; 
	Fmeas.xold = Fmeas.x; Fmeas.yold = Fmeas.y; Fmeas.thetaold = Fmeas.theta; // store old values
	Fmeas.x = inx; Fmeas.y = iny; Fmeas.theta = intheta; // add in new values
	// update velocities
	Fmeas.vx = (inx - Fmeas.xold) / Fmeas.dt;
	Fmeas.vy = (iny - Fmeas.yold) / Fmeas.dt;
	Fmeas.vtheta = (intheta - Fmeas.thetaold) / Fmeas.dt; 

};

void TrackStim::setFlyMeasuredDelta(float dx, float dy, float dtheta, float intime)
{	

	Fmeas.told = Fmeas.tcurr;
	Fmeas.tcurr = intime; // set time of the newest measurements -- could add as argument, i think.
	Fmeas.dt = Fmeas.tcurr - Fmeas.told; 
	Fmeas.xold = Fmeas.x; Fmeas.yold = Fmeas.y; Fmeas.thetaold = Fmeas.theta; // store old values
	Fmeas.x += dx; Fmeas.y += dy; Fmeas.theta += dtheta; // add in new values
	// update velocities
	Fmeas.vx = (dx) / Fmeas.dt;
	Fmeas.vy = (dy) / Fmeas.dt;
	Fmeas.vtheta = (dtheta) / Fmeas.dt; 

};

void TrackStim::setFlyCurrentPos(float currtime)    // NOT THE MEASUREMENT
{
	float dt, dx, dy, dtheta; // to store this time bin's movement...
	
	// store old values...
	Fcurr.told = Fcurr.tcurr;
	Fcurr.xold = Fcurr.x;
	Fcurr.yold = Fcurr.y;
	Fcurr.thetaold = Fcurr.theta;
	
	// update with new values from Fmeas, interpolating if not the same time taken
	Fcurr.tcurr = currtime;
	dt = Fcurr.tcurr - Fmeas.tcurr; // time between last measurement and now
	Fcurr.dt = dt; 
	dx = dt * Fmeas.vx;
	Fcurr.x = Fmeas.x + (dx); // linearly interpolate them all, but keeps up with measured position -- corrects back each time...
	dy = dt * Fmeas.vy;
	Fcurr.y = Fmeas.y + (dy); 
	dtheta = dt * Fmeas.vtheta;
	Fcurr.theta = Fmeas.theta + (dtheta); 

	// these are what should be used in any integration with gain, as they take into account interpolating over and undershoots accurately
	Fcurr.dx = Fcurr.x - Fcurr.xold; // note: these are not the same as above, b/c Fmeas is used to recalibrate with new measurements
	Fcurr.dy = Fcurr.y - Fcurr.yold;
	Fcurr.dtheta = Fcurr.theta - Fcurr.thetaold;
	
};


void TrackStim::writePos()
{
	if (monitor_position_file == NULL)
	{
		monitor_position_file=fopen("_position_output.txt","a");
		fprintf(monitor_position_file,TrackStim::paramfilename);
		fprintf(monitor_position_file,"\n");
		fprintf(monitor_position_file,TrackStim::worldfilename);
		fprintf(monitor_position_file,"\n");
		
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf (monitor_position_file, asctime (timeinfo) );
		
	};
	{	
		// fprintf(monitor_position_file,"%8d %4.3f %5.3f %5.3f %4.3f \n",framenumber,Fcurr.tcurr,Fcurr.x,Fcurr.y,Fcurr.theta);

		fprintf(monitor_position_file,"%8d %4.3f %5.3f %5.3f %4.3f \n",framenumber,Fmeas.tcurr,Fmeas.x,Fmeas.y,Fmeas.theta); // replace Fcurr with Fmeas, 080930 to kill overshoots in fwd reconstruction
		// Fcurr was putting in 1/240 s extra to give a small overshoot on all measurements -- this makes the data a tad noisier, which i don't need.
		// KEEP IT THIS WAY UNLESS THERE'S A GOOD REASON; IN THIS PROGRAM, YOU MEASURE EACH CYCLE, WRITE THIS EACH CYCLE, NO REASON TO USE FCURR
	};
};

void TrackStim::writeStim()
{
	if (monitor_stimulus_file == NULL)
	{
		monitor_stimulus_file=fopen("_stimulus_output.txt","a");
	};
	{	
		fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %4.3f \n",framenumber,Sfull.tcurr,epochchoose,Sfull.x,Sfull.y,Sfull.theta);
		//fprintf(monitor_stimulus_file,"%8d %8.3f %2d %5.4f %5.4f %3.2f \n",framenumber,Sfull.tcurr,epochchoose,Scurr.x,Scurr.dx,Scurr.dt);
	};
};

void TrackStim::writeError()
{
	// bool error = (Fmeas.dt > 0.01); // skipped frame, can expand this definition of errors, make error a global variable
	if (Fmeas.dt > 0.01)
	{
		if (error_file == NULL)
		{
			error_file=fopen("_error_out.txt","a");
		};
		fprintf(error_file,"EXCEEDED dt @ %8d %2d %3.3f\n",framenumber,epochchoose,Fmeas.dt);
	};
};


void TrackStim::updateStimCoords() // something wrong here. not sure what. figure out. seems not to be updating properly with dx from to Scurr.
{
	float totdist = sqrt(Fcurr.dx*Fcurr.dx + Fcurr.dy*Fcurr.dy);
	
	Scurr.tcurr = Fcurr.tcurr;
	if (totdist > 0)
	{
		Scurr.dx = Fcurr.dx * Stimulus.transgain + (Fcurr.dx / totdist * Stimulus.trans.mean * Fcurr.dt);
		Scurr.dy = Fcurr.dy * Stimulus.transgain + (Fcurr.dy / totdist * Stimulus.trans.mean * Fcurr.dt);
	}
	else
	{
		Scurr.dx = 0 ;
		Scurr.dy = 0 ;
	};
	Scurr.dtheta = Fcurr.dtheta * Stimulus.rotgain * -1 + Stimulus.rot.mean * Fcurr.dt; // dac 090226 -- added -1

	Scurr.x += Scurr.dx;
	Scurr.y += Scurr.dy;
	Scurr.theta += Scurr.dtheta;
	Scurr.dt = Fcurr.dt; // same dt as in update, not in measure (ostensibly screen update)
	DistanceAccumulated += sqrt(Scurr.dx*Scurr.dx + Scurr.dy*Scurr.dy); // keep track

	//090218 -- add in to keep track of changes in fwd, lat... -- UNDO following code, from the mouse conversion routine
	/* 
	*ldx = fwd * cos(currang) - lat * sin(currang); // currang is integrated angle
	*ldy = fwd * sin(currang) + lat * cos(currang);
	*/
	Scurr.dfwd = Scurr.dx * cos(Scurr.theta - Scurr.dtheta) - Scurr.dy * sin(-Scurr.theta + Scurr.dtheta);
	Scurr.dlat = Scurr.dx * sin(-Scurr.theta +Scurr.dtheta) + Scurr.dy * cos(Scurr.theta - Scurr.dtheta);
	FwdAccumulated += Scurr.dfwd;

	if (totdist > 0)
	{   // as it stands, this leaves out trans.mean here b/c it is done for Scurr...
		Sfull.x = Scurr.x + (Fcurr.dx / totdist * Stimulus.trans.amp * TrackStim::sineWave(Stimulus.trans.per)); // adds in time dependent stuff here, so it's not integrating it
		Sfull.y = Scurr.y + (Fcurr.dy / totdist * Stimulus.trans.amp * TrackStim::sineWave(Stimulus.trans.per));
	}
	else
	{
		Sfull.x = Scurr.x ; // adds in time dependent stuff here, so it's not integrating it
		Sfull.y = Scurr.y ;
	};
	Sfull.theta = Scurr.theta + (Stimulus.rot.amp * TrackStim::sineWave(Stimulus.rot.per));
	Sfull.theta += Saccumulate.theta;
	Sfull.tcurr = Scurr.tcurr;
	Sfull.dt = Scurr.dt;

	TrackStim::checkBoundaries();  // do this last after each update... ?? reason to keep separate?
};


float TrackStim::mod(float x,float y) // better mod than normal C++ mod? deals with signs appropriately.
{
	if (x>=0)
		return x - y*(int)(x/y);
	else
		return x - y*(int)((x/y)-1);
}


void TrackStim::checkBoundaries()
{
	Sfull.x = mod(Sfull.x,Stimulus.arenasize);
	Sfull.y = mod(Sfull.y,Stimulus.arenasize);
};
//////////////////making real stimuli//////////////////////////

float TrackStim::getDLPColor(float DLPintensity, char channel)
{
	// returns the digital value to get that intensity in that channel; single gamma for both
	//const static float gamma = 2.9;  // check these numbers...
	//const static float lowsat = 30/255;
	//const static float highsat = 200/255;
	//const static float GBratio = 0.83; 


	// changed 081119, after recalibration of projector on this computer
	const static float gamma = 2.15f;  // check these numbers...
	const static float lowsat = 0/255;
	const static float highsat = 247/255;
	const static float GBratio = 0.81f; 
	float temp = 0;
	
	if (DLPintensity > 0)
	{
		if (channel == 'G')
		{
			temp = lowsat + pow(DLPintensity * GBratio, 1/gamma);
		};

		if (channel == 'B')
		{
			temp = lowsat + pow(DLPintensity, 1/gamma);
		};
	};

	if (temp < 1)
		return temp;
	else
	{

	}
		return 1.0f;
};

void TrackStim::setColorTriples()
{
	float bgcol = Stimulus.lum * (1 - Stimulus.contrast);
	float fgcol = Stimulus.lum * (1 + Stimulus.contrast);
	float mncol = Stimulus.lum;
	backgroundtriple[0]=0.0f;
	backgroundtriple[1]=getDLPColor(bgcol,'G');
	backgroundtriple[2]=getDLPColor(bgcol,'B');
	foregroundtriple[0]=0.0f;
	foregroundtriple[1]=getDLPColor(fgcol,'G');
	foregroundtriple[2]=getDLPColor(fgcol,'B');
	meangroundtriple[0]=0.0f;
	meangroundtriple[1]=getDLPColor(mncol,'G');
	meangroundtriple[2]=getDLPColor(mncol,'B');
	blacktriple[0]=blacktriple[1]=blacktriple[2]=0.0f;
};


void TrackStim::setBackground() // sets background color
{
	// float color = pow((Stimulus.lum * (1 - Stimulus.contrast/2)) , (1/GammaCorrect));
	// glClearColor(color, color, color, 0.0f);		
	//
	//float bgcol = Stimulus.lum * (1 - Stimulus.contrast/2);
	//glClearColor(0.0f,getDLPColor(bgcol,'G'),getDLPColor(bgcol,'B'),0.0f);
	//glClearColor(backgroundtriple[0],backgroundtriple[1],backgroundtriple[2],0.0f);
	glClearColor(0.1, 0.1, 0.1, 0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

};  


void TrackStim::setBackgroundColor()
{
	// float color = pow((Stimulus.lum * (1 - Stimulus.contrast/2)) , (1/GammaCorrect));
	// glColor3f(color, color, color);				// 
	
	//float bgcol = Stimulus.lum * (1 - Stimulus.contrast/2);
	//glColor3f(0.0f,getDLPColor(bgcol,'G'),getDLPColor(bgcol,'B'));
	glColor3f(backgroundtriple[0],backgroundtriple[1],backgroundtriple[2]);
};

void TrackStim::setForeground() // sets foreground color
{
	// float color = pow((Stimulus.lum * (1 + Stimulus.contrast/2)) , (1/GammaCorrect));
	// glColor3f(color, color, color);
	
	//float col = Stimulus.lum * (1 + Stimulus.contrast/2);
	//glColor3f(0.0f,getDLPColor(col,'G'),getDLPColor(col,'B'));
	glColor3f(foregroundtriple[0],foregroundtriple[1],foregroundtriple[2]);};  


void TrackStim::readstring(FILE *f,char *string)
{
	do
	{
		fgets(string, 255, f);
	} while ((string[0] == '/') || (string[0] == '\n'));
	return;
}

void TrackStim::readXYZ(char * szFile)   // arena size from elsewhere -- spacing in stimulus
{
	/* to get this to work, use code:
	paramfilename = GetFileName() ;
	for (int ii=0; ii<260; ii++) paramfilekeep[ii]=paramfilename[ii]; // seems very klugey, and below
	 -- and pass this paramfilekeep to this function
	*/

	float x, y, h;
	int numposts;
	FILE *filein;
	char oneline[255];
	
	// add in dialog box here, use global variable to keep world name
	filein = fopen(szFile, "rt");				// File To Load World Data From
	for (int ii=0; ii<260; ii++) TrackStim::worldfilename[ii]=szFile[ii];

	readstring(filein,oneline);
	sscanf(oneline, "NUMPOSTS %d\n", &numposts);
	TrackStim::NumObs = numposts; // update global variable, to be used often

	for (int loop = 0; loop < min(numposts,1000); loop++)
	{
		readstring(filein,oneline);
		sscanf(oneline, "%f %f %f", &x, &y, &h);
		World[loop].x = x*Stimulus.arenasize;   // update world
		World[loop].y = y*Stimulus.arenasize;
		World[loop].z = h*Stimulus.arenaheight;
	}
	fclose(filein);

};   


void TrackStim::readNoise()   // arena size from elsewhere -- spacing in stimulus
{
	/* to get this to work, use code:
	paramfilename = GetFileName() ;
	for (int ii=0; ii<260; ii++) paramfilekeep[ii]=paramfilename[ii]; // seems very klugey, and below
	 -- and pass this paramfilekeep to this function
	*/

	float x;
	int numpoints;
	FILE *filein;
	char oneline[255];
	
	// add in dialog box here, use global variable to keep world name
	filein = fopen("E:/Code/trackball_final/Trackball/Debug/Data/noise.txt", "rt");				// File To Load World Data From

	readstring(filein,oneline);
	sscanf(oneline, "NUMPOINTS %d\n", &numpoints);

	for (int loop = 0; loop < min(numpoints,14400); loop++)
	{
		readstring(filein,oneline);
		sscanf(oneline, "%f", &x);
		noiseValues[loop]=x;
	}
	fclose(filein);

};   

float TrackStim::retrieveNextNoiseValue()
{
	static int count=-1;
	count++;
	return noiseValues[count/2 % 14400]; // noise is sampled at 120Hz, 2 minutes
};



float TrackStim::calculateNextNoiseValue()
{
	// produces gaussian white noise (std = 1) with correlation time Stimulus.tau (in seconds)
	// since this keeps track of only 1 value, you can only produce 1 time series with it -- need to make this into its own package/class
	static float dt=1/ScreenUpdateRate/2;  // updates for each color on each frame... though only written once
	float tau=Stimulus.tau;
	static float x=0; // start at 0
	
	double rn1 = ((double)rand())/((double)RAND_MAX+1.0);            //between 0 and <1
    double rn2 = ((double)rand())/((double)RAND_MAX+1.0);
    double r = sqrt( -2.0*log(1.0 - rn1) ) * cos( 2 * PI* rn2 ); // r is gaussian distributed with std = 1

	x += dt/tau*(r*sqrt(2.0f)/sqrt(dt/tau) - x);

	return x;

};

void TrackStim::readViewPositions()
{
	FILE *filein;
	char oneline[255];
	
	// add in dialog box here, use global variable to keep world name
	if (GetFileAttributes("viewpositions.txt") != INVALID_FILE_ATTRIBUTES)
	{
		filein = fopen("viewpositions.txt", "rt");				// File To Load World Data From

		readstring(filein,oneline);
		sscanf(oneline, "%i %i %i %i", &ViewPorts.x[0], &ViewPorts.x[1], &ViewPorts.x[2], &ViewPorts.x[3]);
		readstring(filein,oneline);
		sscanf(oneline, "%i %i %i %i", &ViewPorts.y[0], &ViewPorts.y[1], &ViewPorts.y[2], &ViewPorts.y[3]);
		readstring(filein,oneline);
		sscanf(oneline, "%i %i %i %i", &ViewPorts.w[0], &ViewPorts.w[1], &ViewPorts.w[2], &ViewPorts.w[3]);
		readstring(filein,oneline);
		sscanf(oneline, "%i %i %i %i", &ViewPorts.h[0], &ViewPorts.h[1], &ViewPorts.h[2], &ViewPorts.h[3]);
		fclose(filein);
	}
	else
		ViewPorts = initializeViewPositions();
};

void TrackStim::writeViewPositions()
{
	FILE *fileout;
	
	fileout = fopen("viewpositions.txt", "w");
	fprintf ( fileout, "%i %i %i %i \n", ViewPorts.x[0], ViewPorts.x[1], ViewPorts.x[2], ViewPorts.x[3] );
	fprintf ( fileout, "%i %i %i %i \n", ViewPorts.y[0], ViewPorts.y[1], ViewPorts.y[2], ViewPorts.y[3] );
	fprintf ( fileout, "%i %i %i %i \n", ViewPorts.w[0], ViewPorts.w[1], ViewPorts.w[2], ViewPorts.w[3] );
	fprintf ( fileout, "%i %i %i %i \n", ViewPorts.h[0], ViewPorts.h[1], ViewPorts.h[2], ViewPorts.h[3] );
	fclose ( fileout );
};

void TrackStim::readScreenPosition()
{
	FILE *filein;
	char oneline[255];
	
	// add in dialog box here, use global variable to keep world name
	if (GetFileAttributes("screenpositions.txt") != INVALID_FILE_ATTRIBUTES)
	{
		filein = fopen("screenpositions.txt", "rt");				// File To Load World Data From

		//readstring(filein,oneline);

		readstring(filein,oneline);
		sscanf(oneline, "%f %f", &ScreenPosition.x[0], &ScreenPosition.x[1]);
		readstring(filein,oneline);
		sscanf(oneline, "%f %f", &ScreenPosition.y[0], &ScreenPosition.y[1]);
		readstring(filein,oneline);
		sscanf(oneline, "%f %f", &ScreenPosition.z[0], &ScreenPosition.z[1]);
		fclose(filein);
	}
	else
	{
		ScreenPosition.x[0]=0; ScreenPosition.x[1]=0;
		ScreenPosition.y[0]=0; ScreenPosition.y[1]=0;
		ScreenPosition.z[0]=0; ScreenPosition.z[1]=0;
	};

};

void TrackStim::writeScreenPosition()
{
	FILE *fileout;
	
	fileout = fopen("screenpositions.txt", "w");
	fprintf ( fileout, "%f %f \n", ScreenPosition.x[0], ScreenPosition.x[1] );
	fprintf ( fileout, "%f %f \n", ScreenPosition.y[0], ScreenPosition.y[1] );
	fprintf ( fileout, "%f %f \n", ScreenPosition.z[0], ScreenPosition.z[1] );
	fclose ( fileout );
};


void TrackStim::distancePulse()
{
	float time = Fcurr.tcurr - epochtimezero;
	if (time < Stimulus.tau)
		DistanceAccumulated += Stimulus.stimtrans.mean/ScreenUpdateRate/2;
};


void TrackStim::rotationPulse()
{
	float time = Fcurr.tcurr - epochtimezero;
	if (time < Stimulus.tau)
		Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2; // Sfull is reset each time stim is updated, so this is NOT adding cumulatively...
	
};

void TrackStim::rotationAccumulate()
{
	float time = Fcurr.tcurr - epochtimezero;
	Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2; // Sfull is reset each time stim is updated, so this is NOT adding cumulatively...
};

void TrackStim::rotationAccumulateWait()
{
	float time = Fcurr.tcurr - epochtimezero;
	if (time > Stimulus.tau2)
		Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2; // Sfull is reset each time stim is updated, so this is NOT adding cumulatively...
};


void TrackStim::rotationNoiseUpdate()
{
	// noise acts like velocity here, not position
	Saccumulate.theta += retrieveNextNoiseValue()*PI/180*Stimulus.stimrot.mean/240; // assumes these come at 240 Hz
};

void TrackStim::rotationNoiseCalcUpdate()
{
	// Saccumulate.theta += calculateNextNoiseValue()*Stimulus.stimrot.mean*PI/180/240; // just sets the theta value...
	// Saccumulate.theta = calculateNextNoiseValue(); // just sets the theta value...
	r1.setdt(1/ScreenUpdateRate/2);
	r1.setTau(Stimulus.tau);
	Saccumulate.theta += r1.Rrandnt() * Stimulus.stimrot.mean*PI/180/240; // sets theta value... STILL needs testing, 081118
	// Saccumulate.theta = r1.Rrandnt();
}

void TrackStim::rotationXYBarCalcUpdate()
{
	static int counter = 0;
	if (counter == 0)
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
		r1.setdt(1/ScreenUpdateRate/2);
		r1.setTau(Stimulus.tau);
		r2.setdt(1/ScreenUpdateRate/2);
		r2.setTau(Stimulus.tau);
	};
	r1.Rrandnt(); // amplitude does it...
	r2.Rrandnt();
	
	counter ++;

};

void TrackStim::rotationDoublePulsePlus()
{
	float time = Fcurr.tcurr - epochtimezero;
	if (time < Stimulus.tau)
		Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2; // Sfull is reset each time stim is updated, so this is NOT adding cumulatively...
	if (time > Stimulus.tau2)
		Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2;

};

void TrackStim::rotationDoublePulseMinus()
{
	float time = Fcurr.tcurr - epochtimezero;
	if (time < Stimulus.tau)
		Saccumulate.theta += PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2; // Sfull is reset each time stim is updated, so this is NOT adding cumulatively...
	if (time > Stimulus.tau2)
		Saccumulate.theta -= PI/180*Stimulus.stimrot.mean/ScreenUpdateRate/2;
};
bool TrackStim::drawScene() // executes the drawing set by setStimType
{
	setColorTriples(); // do once for each draw scene
	
	glEnable(GL_SCISSOR_TEST);
	glViewport(0,0,800,600); // whole window
	glScissor(0,0,800,600);
	glClearColor(0.0f,0.0f,0.0f,0.0f); // set whole screen to be black; scissor testing makes backgrounds good in viewports
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
	glLoadIdentity();									// Reset The View
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // still on whole window

	for (int jj=0; jj<2; jj++)
	{
		
		// MUST UPDATE ONE VIEWPORT AT A TIME? 
		// update the model to match the current position, 2 times within frame... can do this outside loop
		TrackStim::setFlyCurrentPos(Fmeas.tcurr + jj/TrackStim::ScreenUpdateRate/2); // approximation to timing here... roughly correct. better timing from brian wilt?
		TrackStim::updateStimCoords(); // update every time here...		

		switch (Stimulus.stimtype) // dac 090302 moved position of this to after Fcurr update command...
		{
			case 14:
				// nothing
				break;
			case 36:
			case 15:
				rotationPulse();
				break;
			case 16:
				rotationDoublePulsePlus();
				break;
			case 17:
				rotationDoublePulseMinus();
				break;
			case 18:
				distancePulse();
				break;
			case 19:
				rotationNoiseCalcUpdate();
				break;
			case 21:
				updateReichardtIntensity();
				break;
			case 24:
			case 25:
				rotationXYBarCalcUpdate();
				break;
			case 37:
				rotationAccumulate();
				break;
			case 38:
				rotationAccumulateWait();
				break;
			case 43:
				updateReichardtIntensity();
				break;
			case 44:
				updateReichardtIntensity_binary();
				break;
			case 45:
				updateReichardtIntensity_binary();
				break;
			case 46:
				updateReichardtIntensity();
				break;
			case 47:
				updateReichardtIntensity_binary_3bars();
				break;
			case 51:
				updateReichardtIntensity();
				break;
			case 52:
				updateLumIntensity();
				break;
			case 53:
				updateReichardtIntensity();
				break;
			case 54:
				updateReichardtIntensity();
				break;
			case 55:
				updateReichardtIntensity();
				break;
			case 56:
				updateReichardtIntensity_discrete_gaussian();
				break;
			case 57:
				updateReichardtIntensity();
				break;
			case 61:
				updateReichCorr_everytime();
				break;
			case 62:
				updateReichCorr();
				break;
			case 63:
				updateReichCorr_everytime();
				break;
			case 64:
				updateReichardtIntensity_binary();
				break;
			case 65:
				updateReichardtIntensity_binary();
				break;
			case 66:
				updateReichardtIntensity_binary();
				break;
			case 67:
				updateReichCorrGrid_everytime();
				break;
			case 68:
				updateReichardtIntensity_binary();
				break;
			case 69:
				updateReichardtIntensity_binary();
				break;
			case 70:
				updateReichardtIntensity_binary();
				break;
			case 76:
				updateReichardtIntensity_binary();
				break;
			case 84:
				updateReichardtIntensity();
				break;
			case 85:
				updateReichardtIntensity_binary_isolating();
				break;


		};

		for (int ii=0; ii<3; ii++) // go through 3 different window for one time...
		{
			glViewport(ViewPorts.x[ii],ViewPorts.y[ii],ViewPorts.w[ii],ViewPorts.h[ii]);
			glScissor(ViewPorts.x[ii],ViewPorts.y[ii],ViewPorts.w[ii],ViewPorts.h[ii]);

			// new
			glClear(GL_DEPTH_BUFFER_BIT);
			glLoadIdentity();									// Reset The View
			// start required

			switch (jj) // have to do this inside loop
			{
				case 0:  // t1 in blue --- blue first, then green in order of DLP
					glColorMask(GL_FALSE,GL_FALSE,GL_TRUE,GL_TRUE);
					break;
				case 1:  // t2 in green
					glColorMask(GL_FALSE,GL_TRUE,GL_FALSE,GL_TRUE);
					break;
			}

			setBackground();
			setForeground();

			glClear(GL_COLOR_BUFFER_BIT);
			glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
			glLoadIdentity();
			if (VIEWCUBECORNER) // THIS IS NO LONGER USED
			{
				if (ii == 1)
					gluPerspective(90.0f,(float)ViewPorts.w[ii]/(float)ViewPorts.h[ii],0.1f,Stimulus.arenasize);
				else
					gluPerspective(90.0f,(float)ViewPorts.w[ii]/(float)ViewPorts.h[ii],0.1f,Stimulus.arenasize/sqrt(2.0f)); // peculiarity of the planes @ 45 degrees

				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();  // start at origin
				switch (ii) // cycle through view ports: theta is measurement from x axis -- look off into +/- pi/2; theta field in radians, supposedly
				{
					case 2:  // left window
						gluLookAt(Sfull.x + ScreenPosition.x[0]*cos(Sfull.theta) - ScreenPosition.y[0]*sin(Sfull.theta), // i think this is correct direction for rotation
							Sfull.y + ScreenPosition.x[0]*sin(Sfull.theta) + ScreenPosition.y[0]*cos(Sfull.theta), 
							0 + ScreenPosition.z[0],   
							Sfull.x+cos(Sfull.theta+PI/2)   + ScreenPosition.x[0]*cos(Sfull.theta) - ScreenPosition.y[0]*sin(Sfull.theta),
							Sfull.y+sin(Sfull.theta+PI/2)    + ScreenPosition.x[0]*sin(Sfull.theta) + ScreenPosition.y[0]*cos(Sfull.theta),
							1             + ScreenPosition.z[0],    
							0,0,1); // look up!
						break;
					case 1:  // center window
						gluLookAt(Sfull.x,                   Sfull.y,                                          0,   
							Sfull.x+cos(Sfull.theta),        Sfull.y+sin(Sfull.theta),                         0,     
							cos(Sfull.theta + PI/2),         sin(Sfull.theta + PI/2),                          1); // need to change this to 45 degrees later...
						break;
					case 0:  // right window
						gluLookAt(Sfull.x + ScreenPosition.x[1]*cos(Sfull.theta) - ScreenPosition.y[1]*sin(Sfull.theta), 
							Sfull.y + ScreenPosition.x[1]*sin(Sfull.theta) + ScreenPosition.y[1]*cos(Sfull.theta), 
							0 + ScreenPosition.z[1],   
							Sfull.x+cos(Sfull.theta-PI/2)     + ScreenPosition.x[1]*cos(Sfull.theta) - ScreenPosition.y[1]*sin(Sfull.theta),
							Sfull.y+sin(Sfull.theta-PI/2)      + ScreenPosition.x[1]*sin(Sfull.theta) + ScreenPosition.y[1]*cos(Sfull.theta),
							1             + ScreenPosition.z[1],    
							0,0,1); // look UP! -- change back!
						break;
				};
			}
			else // THIS IS EXECUTED
			{
				if (!USE_FRUSTUM) // initial way -- just a 90 degree view	
				{
					gluPerspective(90.0f,(float)ViewPorts.w[ii]/(float)ViewPorts.h[ii],0.1f,Stimulus.arenasize);
				}
				else
				{
					// remember, all screens are rotated, so viewports.h is really the width...
					// take screen height to be 4mm always, so that works, then centered in width,
					// which is proportional to the viewports.h
					//float adjustedwidth = 4.0f/ViewPorts.h[1]*ViewPorts.h[ii]; // commented out 110531, now projecting to same size on all sides
					float adjustedwidth = 4.0f;
					if (TRUE)
					{
						switch (ii)
						{
							case 1: // Center

							case 0:

							case 2:
								// changed this to symmetrical on 111027
								glFrustum(-1.5f, 1.5f, -1.7f, 1.7f, 1.5f, Stimulus.arenasize); // is this working?
								//glFrustum(l, r, b, t, n, f)
								//glFrustum(-5.0f, 5.0f, -5.0f, 0.0f, 1.f, 100); // is this working?
								break;
							//case 2:
								//glFrustum(-3.0f, 1.0f, -2.0f, 2.0f, 2.0f, Stimulus.arenasize); // is this working?
								//glFrustum(-1.0f, 3.0f, -2.0f, 2.0f, 2.0f, Stimulus.arenasize); // is this working?
								//glFrustum(-2.0f, 2.0f, -2.0f, 2.0f, 0.0f, Stimulus.arenasize); // is this working?
								break;
						};
					}
					else // DON'T USE THIS -- DISTORTS MORE, DOESN'T FIX PIXEL ROUNDING ERRORS
					{
						switch (ii)
						{
							case 1:
							case 0:
								glFrustum(-1.0f, 3.0f, -adjustedwidth/2, adjustedwidth/2, 2.0f, Stimulus.arenasize); // is this working?
								break;
							case 2:
								glFrustum(-3.0f, 1.0f, -adjustedwidth/2, adjustedwidth/2, 2.0f, Stimulus.arenasize); // is this working?
								break;
						};
					};
				};

				glMatrixMode(GL_MODELVIEW);
				glLoadIdentity();  // start at origin
				/////////***********NOTE -- CHECK THIS ALL OUT IF DOING FEEDBACK!
				float tempoff=0;
				switch (ii) // cycle through view ports: theta is measurement from x axis -- look off into +/- pi/2; theta field in radians, supposedly
				{
					case 2:  // left window
						/*
						gluLookAt(Sfull.x + ScreenPosition.x[0]*cos(Sfull.theta) - ScreenPosition.y[0]*sin(Sfull.theta), // i think this is correct direction for rotation
							Sfull.y + ScreenPosition.x[0]*sin(Sfull.theta) + ScreenPosition.y[0]*cos(Sfull.theta), 
							ScreenPosition.z[0] + Stimulus.arenaheight*5,   

							Sfull.x+cos(Sfull.theta-PI/2)   + ScreenPosition.x[0]*cos(Sfull.theta) - ScreenPosition.y[0]*sin(Sfull.theta),
							Sfull.y+sin(Sfull.theta-PI/2)   + ScreenPosition.x[0]*sin(Sfull.theta) + ScreenPosition.y[0]*cos(Sfull.theta),
							ScreenPosition.z[0] + Stimulus.arenaheight*5,    // look straight out... // HOIST IT UP to middle?

							cos(Sfull.theta), //+PI here?
							sin(Sfull.theta),
							0); // orient for fiber orientations
						break;
						*/
						gluLookAt(0, 0, 0,
							-1, 0, 0,
							0, 0, 1);
						break;
					case 1:  // center window
						// changed 100320
						//gluLookAt(Sfull.x,                   Sfull.y,                                           Stimulus.arenaheight*5,   
						//	Sfull.x+cos(Sfull.theta),        Sfull.y+sin(Sfull.theta),                          Stimulus.arenaheight*5,     
						//	cos(Sfull.theta-PI/2),        sin(Sfull.theta-PI/2),                          0); // need to change this to 45 degrees later...
						/*gluLookAt(Sfull.x + tempoff*cos(Sfull.theta),                   Sfull.y + tempoff*sin(Sfull.theta),                                           Stimulus.arenaheight*5,   
							Sfull.x+cos(Sfull.theta) + tempoff*cos(Sfull.theta),        Sfull.y+sin(Sfull.theta) + tempoff*sin(Sfull.theta),                          Stimulus.arenaheight*5,     
							-cos(Sfull.theta-PI/2),        -sin(Sfull.theta-PI/2),                          0);*/ // need to change this to 45 degrees later...
						gluLookAt(0, 0, 0,
							0, 0, -1,
							-1, 0, 0);
						break;
					case 0:  // right window
						/*gluLookAt(Sfull.x + ScreenPosition.x[1]*cos(Sfull.theta) - ScreenPosition.y[1]*sin(Sfull.theta), 
							Sfull.y + ScreenPosition.x[1]*sin(Sfull.theta) + ScreenPosition.y[1]*cos(Sfull.theta), 
							ScreenPosition.z[1] + Stimulus.arenaheight*5,   

							Sfull.x+cos(Sfull.theta+PI/2)     + ScreenPosition.x[1]*cos(Sfull.theta) - ScreenPosition.y[1]*sin(Sfull.theta),
							Sfull.y+sin(Sfull.theta+PI/2)      + ScreenPosition.x[1]*sin(Sfull.theta) + ScreenPosition.y[1]*cos(Sfull.theta),
							ScreenPosition.z[1] + Stimulus.arenaheight*5,    

							-cos(Sfull.theta),   // - PI here?
							-sin(Sfull.theta),
							0); // orient for fiber orientations
							*/

						gluLookAt(0, 0, 0,
							1, 0, 0,
							0, 0, -1);
						break;
				};
			};

			//drawWindow(0,PI/4,.2); // for potential use later...
			
			// all set, now draw the stimulus
			switch (Stimulus.stimtype)
			{
			case 0:
				drawPosts(); 
				break;
			case 1:
				drawCeiling();
				break;
			case 2:   // just a bar in the distance
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				//glRotatef(Sfull.theta*180/PI,0,0,1);
				drawBar(); 
				break;
			case 3:   // pure open loop bar
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawBar(); 
				break;
			case 4:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder(Fcurr.tcurr - epochtimezero);  // cylinder OL
				break;
			case 5:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHalves(Fcurr.tcurr - epochtimezero);  // expending and contracting poles
				break;
			case 6:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawBarCeiling(Fcurr.tcurr - epochtimezero);  // bars moving at an angle, above, only, OL
				break;
			case 7:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawInfCorridor(Fcurr.tcurr); 
				break;
			case 8:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				if (ii==0) // only update for one of 3 windows...
					drawLinearRDots(1/TrackStim::ScreenUpdateRate/2);
				else
					drawLinearRDots(0);
				break;
			case 9:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				if (ii==0)
					drawRotRDots(1/TrackStim::ScreenUpdateRate/2); 
				else
					drawRotRDots(0);
				break;
			case 10:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				if (ii==0)
					drawRotRDotsHalves(1/TrackStim::ScreenUpdateRate/2);
				else
					drawRotRDotsHalves(0);
				break;
			case 11:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				if (ii==0)
					drawRotRDotsGrads(1/TrackStim::ScreenUpdateRate/2); 
				else
					drawRotRDotsGrads(0);
				break;
			case 12:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawColorTest();
				break;
			case 13:
				glTranslatef(Sfull.x,Sfull.y,Stimulus.arenaheight*5);
				//glRotatef(Sfull.theta*180/PI,0,0,1);
				drawAlignmentGrid();
				break;
			case 14:
			case 15:
			case 16:
			case 17:
			case 19:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1); // Saccumulate is added into Sfull.theta
				drawCylinder(Fcurr.tcurr - epochtimezero);
				break;
			case 18:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHalves(Fcurr.tcurr - epochtimezero);
				break;
			case 20:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawCheckerboard(Fcurr.tcurr - epochtimezero);
				break;
			case 21:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple();
				break;
			case 22:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_deterministic(Fcurr.tcurr - epochtimezero);
				break;
			case 23:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawDoubleBars(Fcurr.tcurr - epochtimezero);
				break;
			case 24:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawDoubleBars_rand(Fcurr.tcurr - epochtimezero);
				break;
			case 25:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawDoubleBars_rand_indep(Fcurr.tcurr - epochtimezero);
				break;
			case 26:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawBrightCylinder(Fcurr.tcurr - epochtimezero);  // cylinder OL
				break;
			case 27:
				drawCylinder(Fcurr.tcurr - epochtimezero);  // cylinder CL
				break;
			case 28:
				drawSineCylinder(Fcurr.tcurr - epochtimezero);
				break;
			case 29:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawDoubleBars_int(Fcurr.tcurr - epochtimezero);
				break;
			case 30:
				drawScintillate(Fcurr.tcurr - epochtimezero);
				break;
			case 31:
				drawScintillate_dtheta(Fcurr.tcurr - epochtimezero);
				break;
			case 32:
				drawScintillate_complex_dtheta(Fcurr.tcurr - epochtimezero);
				break;			
			case 33:
				drawScintillate_complex_rand_dtheta(Fcurr.tcurr - epochtimezero);
				break;			
			case 34:
				drawScintillate_complex_dtheta_wait(Fcurr.tcurr - epochtimezero);
				break;		
			case 35:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHalves_fwd(Fcurr.tcurr - epochtimezero);  // expending and contracting poles
				break;
			case 36:
				drawCylinder(Fcurr.tcurr - epochtimezero);  // cylinder CL
				break;
			case 37:
				drawRasterRDotRot(Fcurr.tcurr - epochtimezero);
				break;
			case 38:
				drawCylinder(Fcurr.tcurr - epochtimezero);
				break;
			case 39:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawCylinderHalves_fwd(Fcurr.tcurr - epochtimezero);  // expanding and contracting poles
				break;
			case 40:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawRasterRDotHemispheres(Fcurr.tcurr - epochtimezero);
				break;
			case 41:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHemispheres(Fcurr.tcurr - epochtimezero);
				break;
			case 42:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHemispheresGain(Fcurr.tcurr - epochtimezero);
				break;
			case 43:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_translate();
				break;
			case 44:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary();
				break;
			case 45:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_translate_binary();
				break;
			case 46:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_spacing();
				break;
			case 47:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawReichardtSlits_multiple_binary_3bars();
				break;
			case 48:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder_Sawtooth(Fcurr.tcurr - epochtimezero);
				break;
			case 49:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHalves_fwd_Sawtooth(Fcurr.tcurr - epochtimezero);
				break;
			case 50:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder_new(Fcurr.tcurr - epochtimezero);
				break;
			case 51:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_widespace();
				break;
			case 52:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder_Sawtooth_randlum(Fcurr.tcurr - epochtimezero);
				break;
			case 53:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_dots(Fcurr.tcurr - epochtimezero);
				break;
			case 54:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_oneeye();
				break;
			case 55:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_nocenter();
				break;
			case 56:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_nocenter();
				break;
			case 57:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_nocenter_balanced();
				break;
			case 58:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_deterministic_steps(Fcurr.tcurr - epochtimezero);
				break;
			case 59:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawSingleBar(Fcurr.tcurr - epochtimezero);
				break;
			case 60:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder_apertured(Fcurr.tcurr - epochtimezero);  // cylinder OL
				break;
			case 61:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichCorr2Bars_nocenter();
				break;
			case 62:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichCorrFullBars();
				break;
			case 63:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichCorrFullBars();
				break;
			case 64:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary_nocenter();
				break;
			case 65:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_multiple_binary_4corners();
				break;
			case 66:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_multiple_binary_nocenter_spacing();
				break;
			case 67:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichCorrFullGrid();
				break;
			case 68:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary_nocenter_1eye();
				break;
			case 69:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary_nocenter_1eye_quiet();
				break;
			case 70:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary_nocenter_1eye_equalspace();
				break;
			case 71:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars(Fcurr.tcurr - epochtimezero);
				break;
			case 72:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderBars_wait(Fcurr.tcurr-epochtimezero);
				break;
			case 73:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_fwd(Fcurr.tcurr - epochtimezero);
				break;
			case 74:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderBars_wait_fwd(Fcurr.tcurr - epochtimezero);
				break;
			case 75:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_equal(Fcurr.tcurr - epochtimezero);
				break;
			case 76:
				glTranslatef(Sfull.x,Sfull.y,0.0f); // need to translate, as storing values in x and y, so moving!
				drawReichardtSlits_multiple_binary_nocenter_no45();
				break;
			case 77:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinder_continuous(Fcurr.tcurr - epochtimezero);
				break;
			case 78:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_monocular(Fcurr.tcurr - epochtimezero);
				break;
			case 79:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderBars_wait_monocular(Fcurr.tcurr - epochtimezero);
				break;
			case 80:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_competing_edges(Fcurr.tcurr - epochtimezero);
				break;
			case 81:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawScintillate_complex_continuous(Fcurr.tcurr - epochtimezero);
				break;
			case 82:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster(Fcurr.tcurr - epochtimezero);
				break;
			case 83:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_spaced(Fcurr.tcurr - epochtimezero);
				break;
			case 84:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_multiple_nocenter_1static(Fcurr.tcurr - epochtimezero);
				break;
			case 85:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_multiple_binary_nocenter_isolating();
				break;
			case 86:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_adapting_competing(Fcurr.tcurr - epochtimezero);
				break;
			case 87:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_adapting_competing_edgelast(Fcurr.tcurr - epochtimezero);
				break;
			case 88:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_4pt(Fcurr.tcurr - epochtimezero);
				break;
			case 89:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_6pt(Fcurr.tcurr - epochtimezero);
				break;
			case 90:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_6ptlong(Fcurr.tcurr - epochtimezero);
				break;
			case 91:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_staircase(Fcurr.tcurr - epochtimezero);
				break;
			case 92:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_staircase0(Fcurr.tcurr - epochtimezero);
				break;
			case 93:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_offcenter(Fcurr.tcurr - epochtimezero);
				break;
			case 94:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_randomphase(Fcurr.tcurr - epochtimezero);
				break;
			case 95:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawGrayWithVerticalNoise(Fcurr.tcurr - epochtimezero);
				break;
			case 96:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_original_gangster_stagger(Fcurr.tcurr - epochtimezero);
				break;
			case 97:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawSineCylinderWait(Fcurr.tcurr - epochtimezero);
				break;
			case 98:
				drawScintillate_random_posneg(Fcurr.tcurr - epochtimezero);
				break;
			case 99:
				drawScintillate_gradient_posneg(Fcurr.tcurr - epochtimezero);
				break;
			case 100:
				drawScintillate_N2parity_posneg(Fcurr.tcurr - epochtimezero);
				break;
			case 101:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawSineCylinderWait_tanh(Fcurr.tcurr - epochtimezero);
				break;
			case 102:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawScintillate_3pt_allchoice(Fcurr.tcurr - epochtimezero);
				break;
			case 103:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawBarSteps(Fcurr.tcurr - epochtimezero);
				break;
			case 104:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawSymmetricBars(Fcurr.tcurr - epochtimezero);
				break;
			case 105:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawScintillate_N2parity_posneg_waittime(Fcurr.tcurr - epochtimezero);
				break;
			case 106:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawScintillate_nograd_posneg(Fcurr.tcurr - epochtimezero);
				break;
			case 107:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_latinh_test(Fcurr.tcurr - epochtimezero);
				break;
			case 108:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderBars_wait_bw(Fcurr.tcurr-epochtimezero);
				break;
			case 109:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_half(Fcurr.tcurr - epochtimezero);
				break;
			case 110:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawScint_equilum(Fcurr.tcurr - epochtimezero);
				break;
			case 111:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHemispheres_revised(Fcurr.tcurr - epochtimezero);
				break;
			case 112:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_competing_edges_translation(Fcurr.tcurr - epochtimezero);
				break;
			case 113:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHemispheres_revised_graycenter(Fcurr.tcurr - epochtimezero);
				break;
			case 114:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderHemispheres_revised_graycenter_wait(Fcurr.tcurr - epochtimezero);
				break;
			case 115:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawExpandingBars_adapting_competing_edgelast_gray(Fcurr.tcurr - epochtimezero);
				break;
			case 116:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1);
				drawCylinderBars_wait_eye_orientation(Fcurr.tcurr-epochtimezero);
				break;
			case 117:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawSineCylinderWait_2D(Fcurr.tcurr - epochtimezero);
				break;
			case 118:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawLoomAzimuthal(Fcurr.tcurr - epochtimezero);
				break;
			case 119:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawSineCylinderWait_gray(Fcurr.tcurr - epochtimezero);
				break;
			case 120:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawMovingSquare(Fcurr.tcurr - epochtimezero);
				break;
			case 121:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				drawReichardtSlits_TRANSLATE_original_gangster_stagger(Fcurr.tcurr - epochtimezero);
				break;
			case 122:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawUniform(Fcurr.tcurr - epochtimezero);
				break;
			case 123:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawMovingBar_je(Fcurr.tcurr - epochtimezero);
				break;
			case 124:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawCenteringSquare(Fcurr.tcurr - epochtimezero);
				break;
			case 125:
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawMovingSquare_elevation_je(Fcurr.tcurr - epochtimezero);
				break;
			case 126: // added 03152013 by Jennifer Esch
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta,0.0f,0.0f,1.0f);
				drawLoomAzimuthal_je(Fcurr.tcurr - epochtimezero);
				break;

			case 127: // added 131002 by Yvette and Eugene
				glTranslatef(Sfull.x,Sfull.y,0.0f);
				glRotatef(Sfull.theta*180/PI,0,0,1); 
				drawStarField(Fcurr.tcurr - epochtimezero);
				break;
			case 128: // added 131028 by Yvette and Eugene
				//glTranslatef(Sfull.x, Sfull.y, 0.0f);
				//glRotatef(Sfull.theta * 180 / PI, 0, 0, 1);
				drawCylinderBarsES(Fcurr.tcurr - epochtimezero);
				break;
			case 129: // added 131029 by Yvette and Eugene
				//glTranslatef(Sfull.x, Sfull.y, 0.0f);
				//glRotatef(Sfull.theta * 180 / PI, 0, 0, 1);
				drawSimpleBarCW(Fcurr.tcurr - epochtimezero);
				break;
			case 130: // added 131029 by Yvette and Eugene
				//glTranslatef(Sfull.x, Sfull.y, 0.0f);
				//glRotatef(Sfull.theta * 180 / PI, 0, 0, 1);
				drawSimpleBarCCW(Fcurr.tcurr - epochtimezero);
				break;
			case 131: // added 131029 by Yvette and Eugene
				//glTranslatef(Sfull.x, Sfull.y, 0.0f);
				//glRotatef(Sfull.theta * 180 / PI, 0, 0, 1);
				drawSimpleBarCW_rphi(Fcurr.tcurr - epochtimezero);
				break;
			case 132: // added 131029 by Yvette and Eugene
				//glTranslatef(Sfull.x, Sfull.y, 0.0f);
				//glRotatef(Sfull.theta * 180 / PI, 0, 0, 1);
				drawSimpleBarCCW_rphi(Fcurr.tcurr - epochtimezero);
				break;
			};

		};

	};
	
	glClear(GL_DEPTH_BUFFER_BIT);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	// 120510 -- this drawmap is flickering for unknown reasons. having to do with 2d texture? seems fine otherwise.
	// though it's intended to draw in white, but it draws only in B,G channels, in spite of color mask command here
	TrackStim::drawMap();   // draw this last...

	glDisable(GL_SCISSOR_TEST); // for some reason, want to do this before swap buffer operation

	return TRUE;	
};  

void TrackStim::drawBackgroundCube()
{
	// obsolete, now that we can get glclear to work within each viewport, black outside...
	float d = Stimulus.arenasize*.999f;
	TrackStim::setBackgroundColor();  // set this up to draw in background color only...
	
	//glTranslatef(Sfull.x,Sfull.y,0.0f); // translate outside the function
	GLUquadricObj *cyl = gluNewQuadric();
	glPushMatrix();
	gluSphere(cyl,d,18,10);
	glPopMatrix();

	//glPushMatrix();
	//gluCylinder(cyl,d*.998,d*.998,Stimulus.arenaheight*2,18,1);
	//glPopMatrix();


};

void TrackStim::drawWindow(float azim, float zenith, float solid)
{
	//TrackStim::setBackground();
	glColor3f(blacktriple[0],blacktriple[1],blacktriple[2]);
	float d = 10.0f;
	glPushMatrix();
	glTranslatef(Sfull.x,Sfull.y,0.0f);
	glRotatef(Sfull.theta*180/PI,0,0,1);
	glRotatef(azim*180/PI,0,0,1);
	glRotatef(zenith*180/PI,0,1,0); // rotate so that z is in correct direction

	GLUquadricObj *cyl = gluNewQuadric();

	glPushMatrix();
	glTranslatef(0.0f,0.0f,-d);
	gluCylinder(cyl,d,d,2*d,10,1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0f,0.0f,-d);
	gluDisk(cyl,0,d,10,1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0.0f,0.0f,d);
	gluDisk(cyl,solid*d,d,10,1);
	glPopMatrix();

	glPopMatrix();
	TrackStim::setForeground();
};

void TrackStim::drawMap() // this is drawn in the upper right hand portion of teh screen. can ignore.
{
	static int maxpos = 500;  // don't update, just keep once
	float Zpos = -1000;

	// update map
	Map.x[Map.currpos]=Sfull.x;
	Map.y[Map.currpos]=Sfull.y;
	Map.currpos++;
	Map.currpos = Map.currpos % maxpos;  // need to mod it here so doesn't run over

	// set up view port and location
	glLoadIdentity();									// Reset The View
	glViewport(ViewPorts.x[3],ViewPorts.y[3],ViewPorts.w[3],ViewPorts.h[3]);
	glScissor(ViewPorts.x[3],ViewPorts.y[3],ViewPorts.w[3],ViewPorts.h[3]);
	glClearColor(0.0f,0.0f,0.0f,0.0f);
	glClear(GL_COLOR_BUFFER_BIT); 
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();
	gluPerspective(90.0f,(float)ViewPorts.w[3]/(float)ViewPorts.h[3],0.1f,100000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();  // start at origin
	gluLookAt(Stimulus.arenasize/2,Stimulus.arenasize/2,Zpos+Stimulus.arenasize/sqrt(2.0f),Stimulus.arenasize/2,Stimulus.arenasize/2,Zpos,0,-1,0);

	// draw that shit -- should be pretty quick
	glColor3f(1.0f,1.0f,1.0f); // draw map always in white
	glBegin(GL_LINE_STRIP);
	for (int ii=0; ii<Map.currpos-1; ii++)
	{
		glVertex3f(Map.x[ii], Map.y[ii], Zpos); // origin of the FIRST line
	};
	glEnd();

	glBegin(GL_LINE_STRIP);
	for (int ii=Map.currpos; ii<maxpos-1; ii++)
	{
		glVertex3f(Map.x[ii], Map.y[ii], Zpos); // origin of the FIRST line
	};
	glEnd();
	
	int choosecurrent=Map.currpos; choosecurrent--; choosecurrent %= maxpos;
	glBegin(GL_LINES);
		glVertex3f(Map.x[choosecurrent]+5*cos(Sfull.theta),Map.y[choosecurrent]+5*sin(Sfull.theta),Zpos);
		glVertex3f(Map.x[choosecurrent]-5*cos(Sfull.theta),Map.y[choosecurrent]-5*sin(Sfull.theta),Zpos);
	glEnd();
	
	glBegin(GL_LINE_LOOP); // draw out square to define arena
	glVertex3f(0.0f,0.0f,Zpos);
	glVertex3f(Stimulus.arenasize,0.0f,Zpos);
	glVertex3f(Stimulus.arenasize,Stimulus.arenasize,Zpos);
	glVertex3f(0.0f,Stimulus.arenasize,Zpos);
	glEnd();
	
};
void TrackStim::drawPosts() // draws posts if post file read in. not that useful.
{
	GLUquadricObj *cyl = gluNewQuadric();
	for (int ii=-1; ii<2; ii++)               // 2 outer loops to draw extra worlds around boundaries...
	{
		for (int jj=-1; jj<2; jj++)
		{	
			glPushMatrix();
			glTranslatef(ii*Stimulus.arenasize,jj*Stimulus.arenasize,0.0f);
			for (int kk=0; kk<NumObs; kk++)
			{
				glPushMatrix(); // save (0,0,0);
				glTranslatef(World[kk].x,World[kk].y, 2.0f);   // offset in z?
				gluCylinder(cyl, .2, 1.5, World[kk].z, 3, 1);
				glPopMatrix(); // return to (0,0,0);
			};
			glPopMatrix();
		};
	};

};

void TrackStim::drawPosts2()
{
	// almost works -- periodic boundaries are causing problems, though.
	float angle;
	for (int ii=-1; ii<2; ii++)               // 2 outer loops to draw extra worlds around boundaries...
	{
		for (int jj=-1; jj<2; jj++)
		{	
			glPushMatrix();
			glTranslatef(ii*Stimulus.arenasize,jj*Stimulus.arenasize,0.0f);
			for (int kk=0; kk<NumObs; kk++)
			{
				angle = atan((Sfull.y - World[kk].y) / (Sfull.x - World[kk].x)) + PI/2;
				glPushMatrix(); // save (0,0,0);
				glTranslatef(World[kk].x,World[kk].y, 2.0f);   // offset in z?
				glRotatef(angle*180/PI,0.0f,0.0f,1.0f); // rotate this much, just to be sure
				glBegin(GL_TRIANGLES);
				glVertex3f(0.0f,0.0f,0.0f);
				glVertex3f(0.75f,0.0f,World[kk].z);
				glVertex3f(-0.75f,0.0f,World[kk].z);
				glEnd();
				glPopMatrix(); // return to (0,0,0);
			};
			glPopMatrix();
		};
	};

};

void TrackStim::drawAlignmentGrid() // USEFUL
{
	glColor3f(0.0f,getDLPColor(1,'G'),getDLPColor(1,'B'));
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);  
	GLUquadricObj *cyl = gluNewQuadric();
	glPushMatrix();
	gluSphere(cyl,30,36,20);
	glPopMatrix();

};
void TrackStim::drawCeiling() // not useful with current screen positions
{
	int ii,jj,kk,ll;
	
	glBegin(GL_QUADS);
	for (ii=-1; ii<2; ii++)               // 2 outer loops to draw extra worlds around boundaries...
	{
		for (jj=-1; jj<2; jj++)
		{	
			for (kk=0; kk<Stimulus.arenasize/Stimulus.spacing/2; kk++)
			{
				for (ll=0; ll<Stimulus.arenasize/Stimulus.spacing/2; ll++)
				{
					glPushMatrix(); // save (0,0,0);
					glTranslatef(ii*Stimulus.arenasize + kk*Stimulus.spacing*2,
						jj*Stimulus.arenasize + ll*Stimulus.spacing*2,10.0f);
					glVertex3f(0.0f,0.0f,0.0f);
					glVertex3f(0.0f,Stimulus.spacing,0.0f);
					glVertex3f(Stimulus.spacing,Stimulus.spacing,0.0f);
					glVertex3f(Stimulus.spacing,0.0f,0.0f);
					glVertex3f(Stimulus.spacing,Stimulus.spacing,0.0f);
					glVertex3f(Stimulus.spacing,2*Stimulus.spacing,0.0f);
					glVertex3f(2*Stimulus.spacing,2*Stimulus.spacing,0.0f);
					glVertex3f(2*Stimulus.spacing,Stimulus.spacing,0.0f);
					glPopMatrix(); // return to (0,0,0);
				};
			};
		};
	};
	glEnd();
	
};

void TrackStim::drawUniform(float time)
{
	// for testing gamma correction
	int period = Stimulus.tau;
	float increment = 0.1; // hard code for ease here
	int choosechannels = int(Stimulus.spacing); 
	// 1 is for blue only
	// 2 is for green only
	// 3 is for blue and green together
	bool gammacorrect = bool(int(Stimulus.spacing2));
	// TRUE does current gamma correction
	// false does increments straight, without gamma correcting

	int epoch = int(time/period) % 11;
	float currlum = epoch*0.1;

	float uselum[3]={0.0f,0.0f,0.0f};
	
	if (gammacorrect)
	{
		uselum[0] = 0.0f;
		uselum[1] = getDLPColor(currlum,'G');
		uselum[2] = getDLPColor(currlum,'B');
	}
	else
	{
		uselum[0] = 0.0f;
		uselum[1] = currlum;
		uselum[2] = currlum;
	};

	if (choosechannels == 1)
	{
		uselum[1]=0.0f;
	};
	if (choosechannels == 2)
	{
		uselum[2]=0.0f;
	};

	glColor3f(uselum[0],uselum[1],uselum[2]);
	for (int ii=0; ii<36; ii++)
	{
		drawWedge(ii*36.0f,36.0f);
	};


};

void TrackStim::drawWedge(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

};

void TrackStim::drawWedge_mirror(GLfloat angle, GLfloat width)  // ** mirror symmetric version of above
{
	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();
	
	// draw mirror symmetric one. should just work. great!
	angle *= -1;
	width *= -1;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();


};

void TrackStim::drawSquare(float sqedge, float distance, float angle) // for jenn esch
{
	glPushMatrix();
	
	glRotatef(angle,0,0,1);

	GLfloat heightfactor = 10.0f;

	glBegin(GL_QUADS);
		glVertex3f(distance,-sqedge/2, -sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,-sqedge/2, sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, -sqedge/2 + Stimulus.arenaheight*heightfactor/2);
	glEnd();

	glPopMatrix();

};

void TrackStim::drawSquare_elevation(float sqedge, float elevation, float distance, float angle) // for jenn esch
{
	glPushMatrix();
	
	glRotatef(angle,0,0,1);

	GLfloat heightfactor = 10.0f;

	glBegin(GL_QUADS);
		glVertex3f(distance,-sqedge/2, elevation + -sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,-sqedge/2, elevation + sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, elevation + sqedge/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, elevation + -sqedge/2 + Stimulus.arenaheight*heightfactor/2);
	glEnd();

	glPopMatrix();

};

void TrackStim::drawBar_je(float sqedge, float ledge, float distance, float angle) // another one for JE
{
	glPushMatrix();
	
	glRotatef(angle,0,0,1);

	GLfloat heightfactor = 10.0f;

	
	glBegin(GL_QUADS);
		glVertex3f(distance,-sqedge/2, -(ledge*sqedge)/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,-sqedge/2, (ledge*sqedge)/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, (ledge*sqedge)/2 + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(distance,sqedge/2, -(ledge*sqedge)/2 + Stimulus.arenaheight*heightfactor/2);
	glEnd();

	glPopMatrix();

};


void TrackStim::drawWedge_horizontal(GLfloat angle, GLfloat width, GLfloat height, GLfloat added)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	if (added > 0) // get into range of [-height,height]
		added -= int(added/2/height)*2*height;
	else
		added += int(-added/2/height + 1)*2*height;
	added -= height;

	if (added > 12)  // don't let it get out of bounds, biatch
		added = 12;
	if (added < -12)
		added = -12; 
	Sfull.x = added;
	added *= PI/180;
	
	angle *= PI/180;
	width *= PI/180;
	height *= PI/180; // get it all in radians...

	float lowang = -75*PI/180; // pm 60 degrees should work
	float highang = 75*PI/180;
	float MIDDLE = Stimulus.arenaheight*heightfactor/2;


	float curroff;
	float currlow;
	float currhigh;
	glBegin(GL_QUADS);
		for (curroff = lowang; curroff < highang; curroff += 2*height) 
		{
			currlow = radius * tan(curroff + added);
			currhigh = radius * tan(curroff + height + added);
			glVertex3f(radius*cos(angle),radius*sin(angle), MIDDLE + currlow); // low points
			glVertex3f(radius*cos(angle+width),radius*sin(angle+width), MIDDLE + currlow);
			glVertex3f(radius*cos(angle+width),radius*sin(angle+width), MIDDLE + currhigh); // high points
			glVertex3f(radius*cos(angle),radius*sin(angle), MIDDLE + currhigh);
		};
	glEnd();

};

void TrackStim::drawWedge_above0(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	if (angle < 0)
	{
		if (width + angle < 0)
		{
			return;
		};
		width = width + angle;
		angle = 0; // just away from 0
	};
	if ((width + angle < 0) && (angle > 0))
	{
		width = -angle; // just back to 0
	};


	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

};

void TrackStim::drawWedge_below0(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	if (-angle < 0)
	{
		if (-(width + angle) < 0)
		{
			return;
		};
		width = width + angle;
		angle = 0; // just away from 0
	};
	if ((-(width + angle) < 0) && (-angle > 0))
	{
		width = -angle; // just back to 0
	};


	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

};

void TrackStim::drawWedge_inside_top(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=25.0f;
	GLfloat heightfactor = 10.0f;

	GLfloat OFFSET = 0;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), OFFSET + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), OFFSET + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), OFFSET + Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), OFFSET + Stimulus.arenaheight*heightfactor);
	glEnd();

};

void TrackStim::drawWedge_inside_bottom(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=25.0f;
	GLfloat heightfactor = 10.0f;

	GLfloat OFFSET = 0;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), OFFSET + 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), OFFSET + 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), OFFSET + Stimulus.arenaheight*heightfactor/2);
		glVertex3f(radius*cos(angle),radius*sin(angle), OFFSET + Stimulus.arenaheight*heightfactor/2);
	glEnd();

};

void TrackStim::drawWedge_inside(GLfloat angle, GLfloat width)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=25.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

};

void TrackStim::drawSawtooth(GLfloat angle, GLfloat width, GLfloat STphase, GLfloat NDphase) // ** draws a saw tooth in place of wedge...
{
	static bool FIRST = TRUE;
	static float imagedata[128*3];
	if (FIRST)
	{
		for (int ii=0; ii<128; ii++)
		{
			imagedata[ii*3]=0.0f;
			imagedata[ii*3+1]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*(ii/128.0f-.5)*2),'G');
			imagedata[ii*3+2]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*(ii/128.0f-.5)*2),'B');
		}
		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);
	};
	FIRST = FALSE;

	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	

	glBegin(GL_QUADS);
		glTexCoord1f(STphase);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f(NDphase);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glTexCoord1f(NDphase);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(STphase);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

	glDisable(GL_TEXTURE_1D);

};



void TrackStim::drawSawtooth_setlum(GLfloat angle, GLfloat width, GLfloat STphase, GLfloat NDphase, GLfloat lumoffset) // ** draws a saw tooth in place of wedge...
{
	static float imagedata[128*3];


	for (int ii=0; ii<128; ii++)
	{
		imagedata[ii*3]=0.0f;
		imagedata[ii*3+1]=getDLPColor(Stimulus.lum + lumoffset + Stimulus.lum * Stimulus.contrast*(ii/128.0f-.5)*2,'G');
		imagedata[ii*3+2]=getDLPColor(Stimulus.lum + lumoffset + Stimulus.lum * Stimulus.contrast*(ii/128.0f-.5)*2,'B');
	};
	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);


	GLfloat radius=30.0f;
	GLfloat heightfactor = 10.0f;

	angle *= PI/180;
	width *= PI/180;

//	glEnable(GL_TEXTURE_1D);
//	glBindTexture(GL_TEXTURE_1D, 1);
	

	glBegin(GL_QUADS);
		glTexCoord1f(STphase);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f(NDphase);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), 0.0f);
		glTexCoord1f(NDphase);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(STphase);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	glEnd();

	glDisable(GL_TEXTURE_1D);

};


void TrackStim::drawWedgeSquare(GLfloat angle, GLfloat width, GLfloat height)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;

	angle *= PI/180;
	width *= PI/180;

	glBegin(GL_QUADS);
		glVertex3f(radius*cos(angle),radius*sin(angle), height);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), height);
		glVertex3f(radius*cos(angle+width),radius*sin(angle+width), height + radius*width/2);
		glVertex3f(radius*cos(angle),radius*sin(angle), height + radius*width/2);
	glEnd();

};

void TrackStim::drawWedgeSquare_nobegin(GLfloat angle, GLfloat width, GLfloat height)  // ** draws wedge: everything in degrees, till this point
{
	GLfloat radius=30.0f;

	angle *= PI/180;
	width *= PI/180;

	glVertex3f(radius*cos(angle),radius*sin(angle), height);
	glVertex3f(radius*cos(angle+width),radius*sin(angle+width), height);
	glVertex3f(radius*cos(angle+width),radius*sin(angle+width), height + radius*width);
	glVertex3f(radius*cos(angle),radius*sin(angle), height + radius*width);

};

void TrackStim::drawBar()
{
	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	// glTranslatef(Sfull.x,Sfull.y,0.0f);
	// glRotatef(added,0.0f,0.0f,1.0f);
	TrackStim::drawWedge(0.0f + added,Stimulus.spacing);
};
void TrackStim::drawCylinder(float tc)   // pinwheel 
{
	int ii;
	//float tc = TrackStim::queryCurrTime() - epochtimezero;
	float angle = Stimulus.spacing;  // in degrees
	//float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	float added = 0; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//added += Saccumulate.theta;
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency
	for (ii=0; ii<360/angle/2; ii++)
	{
		TrackStim::drawWedge(ii*angle*2 + added,angle);
	}

};

void TrackStim::drawLoomAzimuthal(float tc)   // pinwheel 
{
	float sqedge = Stimulus.spacing;  // in crazy units
	float added = tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	float fardist = 30;
	float currdist = fardist - added;
	if (currdist < 2) // collision -- move square out of view
	{
		currdist = 1099;
	};

	//Scurr.x = currdist;
	//Scurr.y = Stimulus.spacing2;

	TrackStim::drawSquare(sqedge, currdist, Stimulus.stimrot.phase);
	//TrackStim::drawWedge(Stimulus.stimrot.phase,10);

};

void TrackStim::drawLoomAzimuthal_je(float tc)   // this function differs from the original azimuth loom in that a square appears, is still, and then moves
{ 
	float sqedge = Stimulus.spacing;  // in crazy units
	float radius = 30;
	float added = 0; //used to animate dot
	
	if (tc>Stimulus.tau)
	{
		added = (tc-Stimulus.tau)*Stimulus.stimtrans.mean;  // Stimulu.stimtrans.mean is the way to move the square closer or farther
	};
	if (tc>Stimulus.tau2)
	{
		added = (Stimulus.tau2-Stimulus.tau)*Stimulus.stimtrans.mean;
	};

	float currdist = radius - added;
	if (currdist < 2) // collision -- move square out of view
	{
		currdist = 1099;
	};

	TrackStim::drawSquare(sqedge, currdist, Stimulus.stimrot.phase);
	

};

void TrackStim::drawMovingSquare(float tc)    
{
	float radius = 30;
	float sqedge = 2*radius*tan(Stimulus.spacing*PI/180/2);  // spacing in degrees
	float added = 0;

	if (tc>Stimulus.tau)
	{
		added = (tc-Stimulus.tau)*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	};
	if (tc>Stimulus.tau2)
	{
		added = (Stimulus.tau2-Stimulus.tau)*Stimulus.stimtrans.mean;
	};


	//Scurr.x = currdist;
	//Scurr.y = Stimulus.spacing2;

	TrackStim::drawSquare(sqedge, radius, Stimulus.stimrot.phase + added);

};

void TrackStim::drawCenteringSquare(float tc)    
{
	float radius = 30;
	float sqedge = 2*radius*tan(Stimulus.spacing*PI/180/2);  // spacing in degrees
	float left = -90;
	float right = 90;
	float center = 0;
	//Scurr.x = currdist;
	//Scurr.y = Stimulus.spacing2;

	TrackStim::drawSquare(sqedge, radius, left);
	TrackStim::drawSquare(sqedge, radius, right);
	TrackStim::drawSquare(sqedge, radius, center);

};

void TrackStim::drawMovingBar_je(float tc)    
{
	float radius = 30;
	float sqedge = 2*radius*tan(Stimulus.spacing*PI/180/2);  // spacing in degrees
	float ledge = Stimulus.spacing2;
	float added = 0;

	if (tc>Stimulus.tau)
	{
		added = (tc-Stimulus.tau)*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	};
	if (tc>Stimulus.tau2)
	{
		added = (Stimulus.tau2-Stimulus.tau)*Stimulus.stimtrans.mean;
	};


	//Scurr.x = currdist;
	//Scurr.y = Stimulus.spacing2;

	TrackStim::drawBar_je(sqedge, ledge, radius, Stimulus.stimrot.phase + added);

};

void TrackStim::drawMovingSquare_elevation_je(float tc)    
{
	float radius = 30;
	float sqedge = 2*radius*tan(Stimulus.spacing*PI/180/2);  // spacing in degrees
	float elevation = Stimulus.spacing2;
	float added = 0;

	if (tc>Stimulus.tau)
	{
		added = (tc-Stimulus.tau)*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	};
	if (tc>Stimulus.tau2)
	{
		added = (Stimulus.tau2-Stimulus.tau)*Stimulus.stimtrans.mean;
	};

	//Scurr.x = currdist;
	//Scurr.y = Stimulus.spacing2;

	TrackStim::drawSquare_elevation(sqedge, elevation, radius, Stimulus.stimrot.phase + added);

};

void TrackStim::drawCylinder_new(float time)   // pinwheel 
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added = 0; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	added = time*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	for (ii=0; ii<360/angle/2; ii++)
	{
		TrackStim::drawWedge(ii*angle*2 + added,angle);
	}

};

void TrackStim::drawCylinder_apertured(float time)   // pinwheel 
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added = time*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	for (ii=0; ii<360/angle/2; ii++)
	{
		TrackStim::drawWedge(ii*angle*2 + added,angle);
	};

	
	glColor3f(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2]);

	float start1 = mod(Stimulus.aperture.location + Stimulus.aperture.width/2,360);
	float start2 = mod(Stimulus.aperture.location - Stimulus.aperture.width/2,360);
	if (start1>start2)
	{
		for (float wedge = start1; wedge < 360 - 0.01; wedge = wedge + (360 - start1)/10) 
			TrackStim::drawWedge_inside(wedge,(360-start1)/10);

		for (float wedge = start2; wedge > 0.01; wedge = wedge + (0 - start2)/10) 
			TrackStim::drawWedge_inside(wedge,-start2/10);
	}
	else
	{
		for (float wedge = start1; wedge < start2 - 0.01; wedge = wedge + (360 - Stimulus.aperture.width)/10) 
			TrackStim::drawWedge_inside(wedge,(360 - Stimulus.aperture.width)/10);

	};
	


};

void TrackStim::drawExpandingBars(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 2; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added); // added < 0
		};
	};
	
};

void TrackStim::drawExpandingBars_half(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 2; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added); // added < 0
		};
	};

	// draw gray for top or bottom half...
	glColor3f(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2]);
	if (Stimulus.stimrot.phase < 0.5)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge_inside_top(ii*angle,angle); 
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge_inside_bottom(ii*angle,angle); // added < 0
		};
	};

};

void TrackStim::drawExpandingBars_latinh_test(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing2; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	static bool RESHIFT=TRUE;
	static float offset = 0;

	if ((tc<0.001) && (RESHIFT))
	{
		offset = rand() % int(angle);
		RESHIFT = FALSE;
	};
	
	if (tc>0.1)
	{
		RESHIFT = TRUE; // wait till next round
	};



	if (tc<=Stimulus.tau2)
	{
		added=0;
	}
	else
	{
		if (tc<=Stimulus.tau)
		{
			added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
		}
		else
		{
			added = (Stimulus.tau-Stimulus.tau2)*Stimulus.stimrot.mean; // stop it after Stimulus.tau seconds
		};
	};

	if (Stimulus.stimrot.mean>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(offset+ii*angle-startwidth/2,startwidth+added);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(offset+ii*angle-startwidth/2+added,startwidth-added); // added < 0
		};
	};
	
};

void TrackStim::drawExpandingBars_staircase(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 2; // starts at 2 degrees
	float step = Stimulus.spacing2;
	float angle = Stimulus.spacing;  // in degrees
	float added,addedfull; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>=0)
	{
		addedfull = startwidth + int(added/step)*step;
	}
	else
	{
		addedfull = startwidth + int((-added)/step)*step;
	};

	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,addedfull);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle+startwidth/2-addedfull,addedfull); // added < 0
		};
	};
	
};

void TrackStim::drawExpandingBars_staircase0(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 0; // starts at 2 degrees
	float step = Stimulus.spacing2;
	float angle = Stimulus.spacing;  // in degrees
	float added,addedfull; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>=0)
	{
		addedfull = startwidth + int(added/step)*step;
	}
	else
	{
		addedfull = startwidth + int((-added)/step)*step;
	};

	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,addedfull);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle+startwidth/2-addedfull,addedfull); // added < 0
		};
	};
	
};

void TrackStim::drawExpandingBars_adapting_competing(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 0; 
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	bool EXP_EDGE;
	bool REVERSE = bool(int(Stimulus.stimrot.phase));
	float w_even, w_odd;

	if (tc<=Stimulus.tau2)
	{
		added=tc*Stimulus.stimrot.mean;
		EXP_EDGE = TRUE;
		if (added > angle) // pull it back down to zero, mo fo
			added -= int(abs(added)/angle)*angle;
		if (added < -angle)
			added += int(abs(added)/angle)*angle;

		Sfull.x = added;
		Sfull.y = 0;
	}
	else
	{
		// NOTE -- DOUBLING HERE TO MATCH THINGS FROM BEFORE -- WISE? HOPE SO.
		//angle = 2*Stimulus.spacing;
		added = 2*(tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
		EXP_EDGE = FALSE;
		if (REVERSE) {
			added *= -1;
		};
		if (added > 2*angle) {
			added -= 2*angle*int(added/2/angle); // bring it back;
		};
		if (added < 0) {
			added += 2*angle*int(-added/2/angle+1);
		};

		w_even = angle - abs(added - angle);
		w_odd  = 0 + abs(added - angle);
		
		Sfull.x = 0;
		Sfull.y = added;

	};

	if (EXP_EDGE)
	{
		if (added>=0){
			for (ii=0; ii<360/angle; ii++){
				TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
			};
		}
		else{
			for (ii=0; ii<360/angle; ii++){
				TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added);
			};
		};
	}
	else // opposing edges...
	{
		if (Stimulus.stimrot.mean>0){
			for (ii=0; ii<int(360.01/angle); ii++){
				if (ii % 2 == 0)
					TrackStim::drawWedge(ii*angle-angle/2,w_even);
				else
					TrackStim::drawWedge(ii*angle-angle/2,w_odd);
			};
		} 
		else{
			for (ii=0; ii<int(360.01/angle); ii++){
				if (ii % 2 == 0)
					TrackStim::drawWedge(ii*angle-w_even-angle*3/2,w_even);
				else
					TrackStim::drawWedge(ii*angle-w_odd-angle*3/2,w_odd);
			};
		};
	};
	
};

void TrackStim::drawExpandingBars_adapting_competing_edgelast(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 0; 
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	bool EXP_EDGE;
	bool REVERSE = bool(int(Stimulus.stimrot.phase));
	float w_even, w_odd;

	if (tc<=Stimulus.tau2)
	{
		// NOTE -- DOUBLING HERE TO MATCH THINGS FROM BEFORE -- WISE? HOPE SO.
		//angle = 2*Stimulus.spacing;
		added = 2*tc*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
		EXP_EDGE = FALSE;
		if (added > 2*angle) {
			added -= 2*angle*int(added/2/angle); // bring it back;
		};
		if (added < 0) {
			added += 2*angle*int(-added/2/angle+1);
		};

		w_even = angle - abs(added - angle);
		w_odd  = 0 + abs(added - angle);
		
		Sfull.x = 0;
		Sfull.y = added;
	}
	else
	{
		added=(tc-Stimulus.tau2)*Stimulus.stimrot.mean;
		if (REVERSE) {
			added *= -1;
		};

		EXP_EDGE = TRUE;
		if (added > angle) // pull it back down to zero, mo fo
			added -= int(abs(added)/angle)*angle;
		if (added < -angle)
			added += int(abs(added)/angle)*angle;

		Sfull.x = added;
		Sfull.y = 0;
	};

	if (EXP_EDGE)
	{
		if (added>=0){
			for (ii=0; ii<360/angle; ii++){
				TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
			};
		}
		else{
			for (ii=0; ii<360/angle; ii++){
				TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added);
			};
		};
	}
	else // opposing edges...
	{
		if (Stimulus.stimrot.mean>0){
			for (ii=0; ii<int(360.01/angle); ii++){
				if (ii % 2 == 0)
					TrackStim::drawWedge(ii*angle-angle/2,w_even);
				else
					TrackStim::drawWedge(ii*angle-angle/2,w_odd);
			};
		} 
		else{
			for (ii=0; ii<int(360.01/angle); ii++){
				if (ii % 2 == 0)
					TrackStim::drawWedge(ii*angle-w_even-angle*3/2,w_even);
				else
					TrackStim::drawWedge(ii*angle-w_odd-angle*3/2,w_odd);
			};
		};
	};
	
};

void TrackStim::drawExpandingBars_adapting_competing_edgelast_gray(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 0; 
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	bool EXP_EDGE;
	bool REVERSE = bool(int(Stimulus.stimrot.phase));
	float w_even, w_odd;
	bool GRAY = FALSE;

	if (tc<=Stimulus.tau2)
	{
		// NOTE -- DOUBLING HERE TO MATCH THINGS FROM BEFORE -- WISE? HOPE SO.
		//angle = 2*Stimulus.spacing;
		added = 2*tc*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
		EXP_EDGE = FALSE;
		if (added > 2*angle) {
			added -= 2*angle*int(added/2/angle); // bring it back;
		};
		if (added < 0) {
			added += 2*angle*int(-added/2/angle+1);
		};

		w_even = angle - abs(added - angle);
		w_odd  = 0 + abs(added - angle);
		
		Sfull.x = 0;
		Sfull.y = added;
	}
	else
	{
		if (tc >= Stimulus.tau)
		{
			added=(tc-Stimulus.tau2)*Stimulus.stimrot.mean;
			if (REVERSE) {
				added *= -1;
			};

			EXP_EDGE = TRUE;
			if (added > angle) // pull it back down to zero, mo fo
				added -= int(abs(added)/angle)*angle;
			if (added < -angle)
				added += int(abs(added)/angle)*angle;

			Sfull.x = added;
			Sfull.y = 0;
		}
		else
		{
			GRAY = TRUE;
		};
	};

	if (GRAY)
	{
		glClearColor(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2],0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
	}
	else
	{
		if (EXP_EDGE) // expanding edges
		{
			// in this case, do things with different contrast. shit.
			glClearColor(getDLPColor(Stimulus.lum*(1-Stimulus.density*Stimulus.contrast),'R'),
				getDLPColor(Stimulus.lum*(1-Stimulus.density*Stimulus.contrast),'G'),
				getDLPColor(Stimulus.lum*(1-Stimulus.density*Stimulus.contrast),'B'),
				0.0f);
			glClear(GL_COLOR_BUFFER_BIT);
			glColor3f(getDLPColor(Stimulus.lum*(1+Stimulus.density*Stimulus.contrast),'R'),
				getDLPColor(Stimulus.lum*(1+Stimulus.density*Stimulus.contrast),'G'),
				getDLPColor(Stimulus.lum*(1+Stimulus.density*Stimulus.contrast),'B'));

			if (added>=0){
				for (ii=0; ii<360/angle; ii++){
					TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
				};
			}
			else{
				for (ii=0; ii<360/angle; ii++){
					TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added);
				};
			};
		}
		else // opposing edges...
		{
			if (Stimulus.stimrot.mean>0){
				for (ii=0; ii<int(360.01/angle); ii++){
					if (ii % 2 == 0)
						TrackStim::drawWedge(ii*angle-angle/2,w_even);
					else
						TrackStim::drawWedge(ii*angle-angle/2,w_odd);
				};
			} 
			else{
				for (ii=0; ii<int(360.01/angle); ii++){
					if (ii % 2 == 0)
						TrackStim::drawWedge(ii*angle-w_even-angle*3/2,w_even);
					else
						TrackStim::drawWedge(ii*angle-w_odd-angle*3/2,w_odd);
				};
			};
		};
	};
	
};




void TrackStim::drawExpandingBars_competing_edges(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!

	// condition added properly -- BUG -- MUST MAKE IT SO NEGATIVE added works?
	// also, has weirdness with phase, where static lines are...
	if (added > 2*angle) {
		added -= 2*angle*int(added/2/angle); // bring it back;
	};
	if (added < 0) {
		added += 2*angle*int(-added/2/angle+1);
	};

	// fix: don't make positive, just shift into (-angle, angle), then no abs value signs?


	float w_even = angle - abs(added - angle);
	float w_odd  = 0 + abs(added - angle);

	if (Stimulus.stimrot.mean>0)
	{
		for (ii=0; ii<int(360.01/angle); ii++)
		{
			if (ii % 2 == 0)
				TrackStim::drawWedge(ii*angle-angle/2,w_even);
			else
				TrackStim::drawWedge(ii*angle-angle/2,w_odd);
		};
	} 
	else
	{
		for (ii=0; ii<int(360.01/angle); ii++)
		{
			if (ii % 2 == 0)
				TrackStim::drawWedge(ii*angle-w_even-angle*3/2,w_even);
			else
				TrackStim::drawWedge(ii*angle-w_odd-angle*3/2,w_odd);
		};
	};
	
};



void TrackStim::drawExpandingBars_competing_edges_translation(float tc)   // bars that expand outwards one direction or other 
{
	// CAREFUL WITH LARGE WAVELENGTH. I PAID A LOT OF ATTENTION TO THE FRONT POLE
	// BUT NOT THE REAR POLE
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!

	// condition added properly -- BUG -- MUST MAKE IT SO NEGATIVE added works?
	// also, has weirdness with phase, where static lines are...
	if (added > 2*angle) {
		added -= 2*angle*int(added/2/angle); // bring it back;
	};
	if (added < 0) {
		added += 2*angle*int(-added/2/angle+1);
	};

	// fix: don't make positive, just shift into (-angle, angle), then no abs value signs?


	float w_even = angle - abs(added - angle); // range = [-angle,angle]
	float w_odd  = 0 + abs(added - angle);

	if (Stimulus.stimrot.mean>0)
	{
		for (ii=0; ii<int(180.01/angle); ii++)
		{
			if (ii % 2 == 0)
			{
				TrackStim::drawWedge_above0(ii*angle-angle/2,w_even);
				TrackStim::drawWedge_below0(-(ii*angle-angle/2),-w_even); // symmetrical one
			}
			else
			{
				TrackStim::drawWedge_above0(ii*angle-angle/2,w_odd);
				TrackStim::drawWedge_below0(-(ii*angle-angle/2),-w_odd);
			}
		};
	} 
	else
	{
		for (ii=1; ii<int(180.01/angle); ii++) // go from 1, b/c subtracting off angle*3/2
		{
			if (ii % 2 == 0)
			{
				TrackStim::drawWedge_above0(ii*angle-w_even-angle*3/2,w_even);
				TrackStim::drawWedge_below0(-(ii*angle-w_even-angle*3/2),-w_even); // symmetrical one
			}
			else
			{
				TrackStim::drawWedge_above0(ii*angle-w_odd-angle*3/2,w_odd);
				TrackStim::drawWedge_below0(-(ii*angle-w_odd-angle*3/2),-w_odd);

			}
		};
	};
	
};




void TrackStim::drawExpandingBars_monocular(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 2; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	bool RIGHTEYE=(abs(Stimulus.rot.phase - 1)<0.01);

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	float contraintensity = (Stimulus.spacing - startwidth - abs(added))/Stimulus.spacing; // fractional intensity from bg
	float bgcol = Stimulus.lum * (1 - Stimulus.contrast);
	float fgcol = Stimulus.lum * (1 + Stimulus.contrast);
	float contracol = fgcol*(1-contraintensity) + bgcol*contraintensity;
	float contraRGB[3];
	contraRGB[0]=0;
	contraRGB[1]=getDLPColor(contracol,'G');
	contraRGB[2]=getDLPColor(contracol,'B');


	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added);
		};
	};

	glColor3f(contraRGB[0],contraRGB[1],contraRGB[2]);
	if (RIGHTEYE)
	{
		for (int ii=0; ii<180.001/angle-1; ii++)
			TrackStim::drawWedge_inside(ii*angle+180,angle);
	}
	else
	{
		for (int ii=0; ii<180.001/angle-1; ii++)
			TrackStim::drawWedge_inside(ii*angle,angle);
	}
	setForeground();

	
};

void TrackStim::drawExpandingBars_equal(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 2; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>=0)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2-added/2,startwidth+added);
		};
	}
	else
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added/2,startwidth-added);
		};
	};
	
};

void TrackStim::drawExpandingBars_fwd(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = 0; // starts at 2 degrees
	float angle = Stimulus.spacing;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>angle)
		added -= angle*int(added/angle);

	if (added<-angle)
		added -= angle*int(added/angle); // subtract negative here, gets it back...

	if (added>=0)
	{
		TrackStim::drawWedge(-added,startwidth+2*added);
		for (ii=1; ii<180/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,startwidth+added);
		};
		for (ii=180/angle; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2-added,startwidth+added);
		};
	}
	else
	{
		// no middle one here...
		for (ii=1; ii<180/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth-added);
		};
		for (ii=180/angle; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2,startwidth-added);
		};
	};
	
};

void TrackStim::drawCylinderBars_wait(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing; // starts at 2 degrees
	float angle = Stimulus.spacing*2;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	for (ii=0; ii<360.01/angle; ii++)
	{
		TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth);
	};
	
};

void TrackStim::drawCylinderBars_wait_eye_orientation(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing; // starts at 2 degrees
	float angle = Stimulus.spacing*2;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	int chooseorientation = int(Stimulus.rot.phase+0.01);

	if (chooseorientation == 0)
	{
		for (ii=0; ii<360.01/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth);
		};
	}
	else
	{
		for (ii=0; ii<36; ii++)
		{
			TrackStim::drawWedge_horizontal(ii*10,10,startwidth, added);
		};
	};


	// draw gray over one side of screen, as required
	int chooseeye = int(Stimulus.stimrot.phase + 0.01);
	switch (chooseeye)
	{
		case 0:
			break;
		case 1:
			glColor3f(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2]);
			for (ii=0; ii<18; ii++) // 0 to 180
			{
				TrackStim::drawWedge_inside(ii*10,10);
			};
			break;
		case 2:
			glColor3f(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2]);
			for (ii=18; ii<36; ii++) // 0 to 180
			{
				TrackStim::drawWedge_inside(ii*10,10);
			};
			break;
	};

	
};

void TrackStim::drawCylinderBars_wait_bw(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing; // starts at 2 degrees
	float angle = Stimulus.spacing*2;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (tc<Stimulus.tau)
	{
		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth);
		};
	}
	else
	{	// just draw the appropriate color
		if (Stimulus.density > 0.5) // white
			setForeground();
		else
			setBackgroundColor(); // black

		for (ii=0; ii<360/angle; ii++)
		{
			TrackStim::drawWedge(ii*angle,angle);
		};
	};
	
};

void TrackStim::drawBarSteps(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float spatperiod = Stimulus.stimrot.per;
	float initialwidth = Stimulus.spacing; 
	float finalwidth = initialwidth + Stimulus.spacing2; 
	int direction = int(Stimulus.stimrot.phase+0.01); // 0 or 1

	if (direction)
	{
		initialwidth *= -1;
		finalwidth *= -1;
	};

	if (tc<=Stimulus.tau)
	{ // just draw things statically
		for (ii=0; ii<360/spatperiod; ii++)
		{
			TrackStim::drawWedge(ii*spatperiod-initialwidth/2,initialwidth);
		};
	}	
	else
	{
		for (ii=0; ii<360/spatperiod; ii++)
		{
			TrackStim::drawWedge(ii*spatperiod-initialwidth/2,finalwidth);
		};
	}	

};

void TrackStim::drawCylinderBars_wait_monocular(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing; // starts at 2 degrees
	float angle = Stimulus.spacing*2;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	
	bool RIGHTEYE=(abs(Stimulus.rot.phase - 1)<0.01);

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	float contracol = Stimulus.lum;
	float contraRGB[3];
	contraRGB[0]=0;
	contraRGB[1]=getDLPColor(contracol,'G');
	contraRGB[2]=getDLPColor(contracol,'B');

	for (ii=0; ii<360/angle; ii++)
	{
		TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth);
	};

	glColor3f(contraRGB[0],contraRGB[1],contraRGB[2]);
	if (RIGHTEYE)
	{
		for (int ii=0; ii<180.001/angle-1; ii++)
			TrackStim::drawWedge_inside(ii*angle+180,angle);
	}
	else
	{
		for (int ii=0; ii<180.001/angle-1; ii++)
			TrackStim::drawWedge_inside(ii*angle,angle);
	};
	setForeground();

	
};

void TrackStim::drawCylinderBars_wait_fwd(float tc)   // bars that expand outwards one direction or other 
{
	int ii;
	float startwidth = Stimulus.spacing; // starts at 2 degrees
	float angle = Stimulus.spacing*2;  // in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc<=Stimulus.tau2)
		added=0;
	else
		added = (tc-Stimulus.tau2)*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	
	if (added>0)
		added=mod(added,angle);
	if (added<0)
		added=mod(added,angle);

	for (ii=0; ii<180/angle; ii++)  // this gets the correct number
	{
		if (ii*angle + added < 180 - startwidth)
		{
			TrackStim::drawWedge(ii*angle+added, startwidth);
		}
		else
		{
			TrackStim::drawWedge(ii*angle+added, startwidth - added);  // last wedge
			TrackStim::drawWedge(0, added - startwidth);
		};
	}
	// second half of it
	for (ii=0; ii<180/angle; ii++)  // this gets the correct number
	{
		if (ii*angle + added < 180 - startwidth)
		{
			TrackStim::drawWedge(360 - (ii*angle+added), -startwidth);
		}
		else
		{
			TrackStim::drawWedge(360 - (ii*angle+added), -(startwidth - added));  // last wedge
			TrackStim::drawWedge(360 - 0, -(added - startwidth));
		};
	}
};

void TrackStim::drawCylinder_Sawtooth(float time)   // pinwheel 
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added = 0; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	added = time*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	float ST,ND;

	if (abs(Stimulus.stimtrans.phase - 1) < 0.01)
	{
		ST = 0; ND = 1;
	}
	else
	{
		ST = 1; ND = 0;
	};
	
	for (ii=0; ii<360/angle/2; ii++)
	{
		TrackStim::drawSawtooth(ii*angle*2 + added,angle*2,ST,ND);
	}

};

void TrackStim::drawCylinder_Sawtooth_randlum(float time)   // pinwheel 
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added = 0; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...
	added = time*Stimulus.stimrot.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	float ST,ND;

	float currlum = r4.readLastValue() * Stimulus.lumsigma;
	Sfull.x = currlum;

	if (abs(Stimulus.stimtrans.phase - 1) < 0.01)
	{
		ST = 0; ND = 1;
	}
	else
	{
		ST = 1; ND = 0;
	};
	
	for (ii=0; ii<360/angle/2; ii++)
	{
		TrackStim::drawSawtooth_setlum(ii*angle*2 + added,angle*2,ST,ND,currlum);
	}

};

void TrackStim::drawSineCylinder(float tc)
{
	static bool FIRST = TRUE;
	static float imagedata[128*3];
	if (FIRST)
	{
		for (int ii=0; ii<128; ii++)
		{
			imagedata[ii*3]=0.0f;
			imagedata[ii*3+1]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*sin(ii*2*PI/128.0f)),'G');
			imagedata[ii*3+2]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*sin(ii*2*PI/128.0f)),'B');
		}
		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);
	};
	FIRST = FALSE;

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = Stimulus.spacing;
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

};

void TrackStim::drawSineCylinderWait_gray(float tc)
{
	static bool FIRST = TRUE;
	static float imagedata[128*3];
	
	if ((FIRST) & (tc<0.1))
	{	
		for (int ii=0; ii<128; ii++)
		{
			imagedata[ii*3]=0.0f;
			imagedata[ii*3+1]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)),'G');
			imagedata[ii*3+2]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)),'B');
		}
		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;

	};
	if (tc>0.15)
	{
		FIRST=TRUE;
	};

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = Stimulus.spacing;
	float radius = 30;
	float heightfactor = 10;
	float angle;

	static float added = 0;

	if (tc<Stimulus.tau2) // first frame wait
	{
		added = 0;
	}
	else
	{
		added = Stimulus.stimrot.mean*(tc - Stimulus.tau2)*PI/180; // update to correct orientation  
	}; // corrected to rads on 101222, after data collection -- data is much faster than noted... discard.
	//Sfull.theta = added;

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);
	
	if (tc<Stimulus.tau)
	{
		glColor3f(meangroundtriple[0],meangroundtriple[1],meangroundtriple[2]);
		for (int ii=0;ii<18;ii++)
		{
			drawWedge_inside(ii*20,20);
		};
	};

};



void TrackStim::drawSineCylinderWait(float tc)
{
	static bool FIRST = TRUE;
	static float imagedata[128*3];
	
	if ((FIRST) & (tc<0.1))
	{	
		for (int ii=0; ii<128; ii++)
		{
			imagedata[ii*3]=0.0f;
			imagedata[ii*3+1]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)),'G');
			imagedata[ii*3+2]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)),'B');
		}
		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;

	};
	if (tc>0.15)
	{
		FIRST=TRUE;
	};

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = Stimulus.spacing;
	float radius = 30;
	float heightfactor = 10;
	float angle;

	static float added = 0;

	if (tc<Stimulus.tau2) // first frame wait
	{
		added = 0;
	}
	else
	{
		added = Stimulus.stimrot.mean*(tc - Stimulus.tau2)*PI/180; // update to correct orientation  
	}; // corrected to rads on 101222, after data collection -- data is much faster than noted... discard.
	//Sfull.theta = added;

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

};

void TrackStim::drawSineCylinderWait_2D(float tc)
{
	static bool FIRST = TRUE;
	static float imagedata[128][128*3];

	static float lum_last = 0; // save old values here
	static float cont_last = 0;

	if (  ((FIRST) & (tc<0.1)) 
		& ((abs(lum_last - Stimulus.lum)>0.001) || (abs(cont_last - Stimulus.contrast)>0.001))  )
	{	///////// ONLY UPDATE THIS AS REQUIRED -- we're first and lum or cont is different
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<128; jj++)
			{
				imagedata[jj][ii*3]=0.0f;
				imagedata[jj][ii*3+1]=getDLPColor(Stimulus.lum *
					(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)*cos(jj*2*PI/128.0f)),'G');
				imagedata[jj][ii*3+2]=getDLPColor(Stimulus.lum *
					(1 + Stimulus.contrast*cos(ii*2*PI/128.0f)*cos(jj*2*PI/128.0f)),'B');
			}
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 128, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
		lum_last = Stimulus.lum; // set these guys up, so they're saved...
		cont_last = Stimulus.contrast;

	};
	if (tc>0.15)
	{
		FIRST=TRUE;
	};

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = Stimulus.spacing;
	float vertper = Stimulus.spacing2;
	float radius = 30;
	float heightfactor = 10;
	float angle;

	static float added = 0;

	if (tc<Stimulus.tau2) // first frame wait
	{
		added = 0;
	}
	else
	{
		added = Stimulus.stimrot.mean*(tc - Stimulus.tau2)*PI/180; // update to correct orientation  
	}; // corrected to rads on 101222, after data collection -- data is much faster than noted... discard.
	//Sfull.theta = added;
	
	float vertical_periods, vertical_phase;
	if (vertper>1)
	{
		vertical_periods = (Stimulus.arenaheight*heightfactor)/(2*radius*tan(vertper*PI/180/2));
		// with phase, center it at equator, arenaheight/2
		vertical_phase = vertical_periods/2; //?
	}
	else
	{
		vertical_periods = 1;
		vertical_phase = 0;
	};

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord2f(vertical_phase , ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f(vertical_phase , (ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f(vertical_phase + vertical_periods , (ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(vertical_phase + vertical_periods , ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};


void TrackStim::drawSineCylinderWait_tanh(float tc)
{
	static bool FIRST = TRUE;
	static float imagedata[128*3];
	float gain = Stimulus.stimrot.phase;
	if ((FIRST) & (tc<0.1))
	{
		for (int ii=0; ii<128; ii++)
		{
			imagedata[ii*3]=0.0f;
			imagedata[ii*3+1]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*tanh(gain*cos(ii*2*PI/128.0f))/tanh(gain)),'G');
			imagedata[ii*3+2]=getDLPColor(Stimulus.lum *(1 + Stimulus.contrast*tanh(gain*cos(ii*2*PI/128.0f))/tanh(gain)),'B');
		}
		glEnable(GL_TEXTURE_1D);
		glBindTexture(GL_TEXTURE_1D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		//glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 128, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;

	};
	if (tc>0.15)
	{
		FIRST=TRUE;
	};

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = Stimulus.spacing;
	float radius = 30;
	float heightfactor = 10;
	float angle;

	static float added = 0;

	if (tc<Stimulus.tau2) // first frame wait
	{
		added = 0;
	}
	else
	{
		added = Stimulus.stimrot.mean*(tc - Stimulus.tau2)*PI/180; // update to correct orientation  
	}; // corrected to rads on 101222, after data collection -- data is much faster than noted... discard.
	//Sfull.theta = added;

	glEnable(GL_TEXTURE_1D);
	glBindTexture(GL_TEXTURE_1D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord1f((ii+1)*20.0f/per);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord1f(ii*20.0f/per);
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_1D);

};

void TrackStim::drawScintillate(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if ((FIRST) || (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=backgroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=backgroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=backgroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			/*
			//a
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			*/
			//b
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];

			
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean
		for (int ii=0; ii<numsquares; ii++)
		{
			/*
			//a
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			*/
			//b
			imagex[ii] += int(Stimulus.stimtrans.mean); // just the shift...
			imagex[ii] %= xlim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];

			
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};
	
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_random_posneg(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],dotsign[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine

	// DECIDE WHAT TO DO
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	// SET UP THE MATRIX WITH GRAY
	if ((FIRST) || (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			dotsign[ii]=rand() % 2;

			if (dotsign[ii])
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	//Sfull.x=dotsign[0];
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_gradient_posneg(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],dotsign[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	// divide by 2 here, because draw two squares each time...
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH) / 2; //

	// DECIDE WHAT TO DO
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	// SET UP THE MATRIX WITH GRAY
	if ((FIRST) || (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			dotsign[ii]=1;

			int coord1 = 3*(imagex[ii]+imagey[ii]*IMW);
			imagedata[coord1+0]=foregroundtriple[0];
			imagedata[coord1+1]=foregroundtriple[1];
			imagedata[coord1+2]=foregroundtriple[2];

			int coord2 = 3*(imagex[ii]+int(Stimulus.spacing2+0.01)+imagey[ii]*IMW) % (3*IMW*IMH); // should work here.
			if (coord2<0)
			{
				coord2 += (3*IMW*IMH)-1; // correct.
			};
			imagedata[coord2+0]=backgroundtriple[0];
			imagedata[coord2+1]=backgroundtriple[1];
			imagedata[coord2+2]=backgroundtriple[2];


		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	//Sfull.x=dotsign[0];
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_nograd_posneg(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],dotsign[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	// divide by 2 here, because draw two squares each time...
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH) / 2; //

	// DECIDE WHAT TO DO
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	// SET UP THE MATRIX WITH GRAY
	if ((FIRST) || (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			dotsign[ii]=(rand() % 2)*2-1;

			int coord1 = 3*(imagex[ii]+imagey[ii]*IMW);
			if (dotsign[ii]==1)
			{
				imagedata[coord1+0]=foregroundtriple[0];
				imagedata[coord1+1]=foregroundtriple[1];
				imagedata[coord1+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[coord1+0]=backgroundtriple[0];
				imagedata[coord1+1]=backgroundtriple[1];
				imagedata[coord1+2]=backgroundtriple[2];
			}


			int coord2 = 3*(imagex[ii]+int(Stimulus.spacing2+0.01)+imagey[ii]*IMW) % (3*IMW*IMH); // should work here.
			if (coord2<0)
			{
				coord2 += (3*IMW*IMH)-1; // correct.
			};
			if (Stimulus.stimrot.phase*dotsign[ii]>0)
			{
				imagedata[coord2+0]=foregroundtriple[0];
				imagedata[coord2+1]=foregroundtriple[1];
				imagedata[coord2+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[coord2+0]=backgroundtriple[0];
				imagedata[coord2+1]=backgroundtriple[1];
				imagedata[coord2+2]=backgroundtriple[2];
			};


		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	//Sfull.x=dotsign[0];
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_N2parity_posneg(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],dotsign[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine

	// DECIDE WHAT TO DO
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	// SET UP THE MATRIX WITH GRAY
	if ((FIRST) || (UPDATE && (loopnum % (int)(Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			dotsign[ii]=(rand() % 2)*2-1; // plus or minus 1

			if (dotsign[ii]==1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == (int)(Stimulus.tau+0.5)))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=(imagex[ii] + int(Stimulus.spacing2)) % xlim;
			if (imagex[ii]<0)
				imagex[ii] += xlim; // wrap appropriately

			dotsign[ii]*=Stimulus.stimrot.phase; // plus or minus 1

			if (dotsign[ii]==1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	//Sfull.x = float(dotsign[0])/100;

	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_N2parity_posneg_waittime(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],dotsign[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine

	int flashdur = int(Stimulus.tau+0.5);
	int flashint = int(Stimulus.tau2+0.5);
	int totperiod = flashdur*2 + flashint; // all in frames here...

	// DECIDE WHAT TO DO
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	int locloop = loopnum % totperiod;
	bool GOGRAY = ( (locloop == 0) || (locloop == flashdur) );
	bool DRAW1 = ( locloop == 0 );
	bool DRAW2 = ( locloop == flashdur + flashint );

	// SET UP THE MATRIX WITH GRAY
	if ((FIRST) || (UPDATE && GOGRAY))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && DRAW1)   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			dotsign[ii]=(rand() % 2)*2-1; // plus or minus 1

			if (dotsign[ii]==1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	if (UPDATE && DRAW2)   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=(imagex[ii] + int(Stimulus.spacing2)) % xlim;
			if (imagex[ii]<0)
				imagex[ii] += xlim; // wrap appropriately

			dotsign[ii]*=Stimulus.stimrot.phase; // plus or minus 1

			if (dotsign[ii]==1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};

	//Sfull.x = float(dotsign[0])/100;

	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};



void TrackStim::drawScintillate_complex_continuous(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	static int stimchoose = 0; // coding for the stimulus in x!
	static int xshift = 0; // moving left or right
	static float firsttriple[3], secondtriple[3]; // holds the frame contrast

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{ // INITIALIZE
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		/*
		0-3 to teh right; 4-7 to the left
		0 = 1-1; 1 = 0-0; 2 = 1-0; 3 = 0-1;
		*/
		stimchoose = rand() % 8; // 8 possible ones here
		if (stimchoose<4)
		{
			xshift = 1;
		}
		else
		{
			xshift = -1;
		}
		switch (stimchoose % 4) // 0 is 0 or 2
		{
			case 0:
			case 2:
				firsttriple[0]=foregroundtriple[0];
				firsttriple[1]=foregroundtriple[1];
				firsttriple[2]=foregroundtriple[2];
				break;
			case 1:
			case 3:
				firsttriple[0]=backgroundtriple[0];
				firsttriple[1]=backgroundtriple[1];
				firsttriple[2]=backgroundtriple[2];
				break;
		};
		switch (stimchoose % 4) // 0 is 0 or 1
		{
			case 0:
			case 3:
				secondtriple[0]=foregroundtriple[0];
				secondtriple[1]=foregroundtriple[1];
				secondtriple[2]=foregroundtriple[2];
				break;
			case 1:
			case 2:
				secondtriple[0]=backgroundtriple[0];
				secondtriple[1]=backgroundtriple[1];
				secondtriple[2]=backgroundtriple[2];
				break;
		};

		
		for (int ii=0; ii<numsquares; ii++)
		{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2]; 
		};
		
		for (int ii=0; ii<numsquares; ii++)
		{
			//b
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=firsttriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=firsttriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=firsttriple[2];

			
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean

		for (int ii=0; ii<numsquares; ii++)
		{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2]; 
		};
		

		for (int ii=0; ii<numsquares; ii++)
		{
			//b
			imagex[ii] += xshift; // just the shift...
			imagex[ii] %= xlim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=secondtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=secondtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=secondtriple[2];

			
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};
	
	Sfull.x = stimchoose;
	Sfull.y = secondtriple[1];

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};



void TrackStim::drawScintillate_dtheta(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	static float dtheta = 0; // amount moved during update...
	bool AFTERWAIT = (tc > Stimulus.tau2);


	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if ((FIRST) || (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=backgroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=backgroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=backgroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			/*
			//a
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			*/
			//b
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];

			
		};

		dtheta = 0; 

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean
		if (AFTERWAIT)
			dtheta = Stimulus.spacing2; // the shift in degrees
		else
			dtheta = 0;
			

		/*glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);*/
	};
	
	
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += dtheta; // move it around, or not
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added*PI/180;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};


void TrackStim::drawScintillate_complex_dtheta(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	static float dtheta = 0; // amount moved during update...

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++) // initialize array
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			
			//a -- reset old pixels at once
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2];
		}

		for (int ii=0; ii<numsquares; ii++)
		{
			
			//b -- find new pixels
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c -- change new pixels
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];

			
		};

		dtheta = 0; 

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean
		dtheta = Stimulus.spacing2; // the shift in degrees @ this point in time regardless...
			
		if (Stimulus.trans.phase == -1) // in this case, switch to opposite shading for second dot...
		{
			for (int ii=0; ii<numsquares; ii++)
			{
				//c
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};
			
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, 1);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);
		};
	};
		
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += dtheta; // move it around, or not
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added*PI/180;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};




void TrackStim::drawScintillate_complex_dtheta_wait(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	static float dtheta = 0; // amount moved during update...
	bool AFTERWAIT = (tc > Stimulus.tau2);

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if ((FIRST) || (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			/*
			//a
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			*/
			//b
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];

			
		};

		dtheta = 0; 

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean
		if (AFTERWAIT)  // only shift it if we're after the initial period, defined by Stimulus.tau2, as above
			dtheta = Stimulus.spacing2; // the shift in degrees
		else
			dtheta = 0;
			
		if (Stimulus.trans.phase == -1) // in this case, switch to opposite shading for second dot...
		{
			for (int ii=0; ii<numsquares; ii++)
			{
				//c
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};
			
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, 1);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);
		};
	};
	
	
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += dtheta; // move it around, or not
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added*PI/180;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};




void TrackStim::drawScintillate_complex_rand_dtheta(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	//static float imagedata[64][128][3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	static float dtheta = 0; // amount moved during update...

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if ((FIRST) || (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0)))   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  // (a) zero old positions (b) reassign those squares to random positions (c) make new positions white
		for (int ii=0; ii<numsquares; ii++)
		{
			/*
			//a
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			*/
			//b
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			//c
			if (rand() % 2 == 1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};

			
		};

		dtheta = 0; 

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  // in this case, shift everything by Stimulus.stimtrans.mean
		dtheta = Stimulus.spacing2; // the shift in degrees
			
		/*
		if (Stimulus.trans.phase == -1) // in this case, switch to opposite shading for second dot...
		{
			for (int ii=0; ii<numsquares; ii++)
			{
				//c
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};
			
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, 1);
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);
		};
		*/
	};
	
	
	
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += dtheta; // move it around, or not
	//added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	//glTranslatef(Sfull.x,Sfull.y,0.0f);
	//glRotatef(added + Sfull.theta*180/PI,0.0f,0.0f,1.0f); // should do this stuff outside this function for clarity and consistency

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added*PI/180;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScintillate_3pt_allchoice(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	int stimchoose = int(Stimulus.stimrot.phase+0.1); // coding for the stimulus in x!
	// 0,1,2,3 give 11,00,10,01 on diagonals, as usual; add 10 to add a third point 1->2, and 20 to add a third point 2->1
	int xshift = (Stimulus.stimrot.per>0)*2-1; // moving left or right:+/-1
	int pt3choose = int(Stimulus.stimtrans.phase+0.1); // third point: 0, 1, or 2 for (-1,0,1)
	static float firsttriple[3], secondtriple[3], thirdtriple[3]; // holds the frame contrast

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{ // INITIALIZE the whole thing
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};

	// set up all the intensities properly...
	switch (stimchoose % 10) // 0,1,2,3
	{
		case 0:
		case 2:
			firsttriple[0]=foregroundtriple[0];
			firsttriple[1]=foregroundtriple[1];
			firsttriple[2]=foregroundtriple[2];
			break;
		case 1:
		case 3:
			firsttriple[0]=backgroundtriple[0];
			firsttriple[1]=backgroundtriple[1];
			firsttriple[2]=backgroundtriple[2];
			break;
	};
	switch (stimchoose % 10) // 0,1,2,3
	{
		case 0:
		case 3:
			secondtriple[0]=foregroundtriple[0];
			secondtriple[1]=foregroundtriple[1];
			secondtriple[2]=foregroundtriple[2];
			break;
		case 1:
		case 2:
			secondtriple[0]=backgroundtriple[0];
			secondtriple[1]=backgroundtriple[1];
			secondtriple[2]=backgroundtriple[2];
			break;
	};
	switch (pt3choose)
	{
		case 0:
			thirdtriple[0]=backgroundtriple[0];
			thirdtriple[1]=backgroundtriple[1];
			thirdtriple[2]=backgroundtriple[2];
			break;
		case 1:
			thirdtriple[0]=meangroundtriple[0];
			thirdtriple[1]=meangroundtriple[1];
			thirdtriple[2]=meangroundtriple[2];
			break;
		case 2:
			thirdtriple[0]=foregroundtriple[0];
			thirdtriple[1]=foregroundtriple[1];
			thirdtriple[2]=foregroundtriple[2];
			break;
	};

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  
		// GRAY all, don't know previous xshift!
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};

		/*for (int ii=0; ii<numsquares; ii++)
		{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2]; 
		};*/
		
		// NEW POSITIONS, ASSIGN FIRST POINT
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;

			if ((imagex[ii]+xshift < 0) || (imagex[ii]+xshift > xlim)) // CHECK IN RANGE
			{
				imagex[ii] -= xshift;
			};

			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=firsttriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=firsttriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=firsttriple[2];
		};

		// THIRD POINT ADDED NOW, AS SPECIFIED
		if (stimchoose/10 == 2)
		{ 
			for (int ii=0; ii<numsquares; ii++)
			{
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+0]=thirdtriple[0];
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+1]=thirdtriple[1];
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+2]=thirdtriple[2];
			};
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  
		// GRAY EVERYWHERE
		for (int ii=0; ii<numsquares; ii++)
		{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2]; 
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+0]=meangroundtriple[0];
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+1]=meangroundtriple[1];
				imagedata[3*(imagex[ii]+xshift+imagey[ii]*IMW)+2]=meangroundtriple[2];
		};

		// ADD THIRD POINT APPROPRIATELY IN THIS CASE
		if (stimchoose/10 == 1)
		{ 
			for (int ii=0; ii<numsquares; ii++)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=thirdtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=thirdtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=thirdtriple[2];
			};
		};

		// PEFORM THE SHIFT and DRAW PT2
		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii] += xshift; // just the shift...
			imagex[ii] %= xlim;

			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=secondtriple[0];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=secondtriple[1];
			imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=secondtriple[2];
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};
	
	Sfull.x = stimchoose;
	Sfull.y = xshift*(pt3choose+1);

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};

void TrackStim::drawScint_equilum(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64],imagec[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // NOT CORRECT -- GET FROM ANIMATION routine
	int xshift = (Stimulus.stimrot.per>0)*2-1; // moving left or right:+/-1
	int pt3choose = int(Stimulus.stimtrans.phase+0.1); // value of 1 leave old point
	int flipc = floor(Stimulus.stimtrans.per+0.1);

	int n1,n2;
	n1 = int(numsquares/2); // two populations of dots
	n2 = n1;
	
	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		//if (Stimulus.stimtrans.mean > 0.5) // TAKE THIS OUT LATER
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{ // INITIALIZE the whole thing
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
				
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};


	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == 0))   // update only on the taus -- if even, just move, otherwise, reassign
	{  
		// GRAY all, don't know previous xshift!
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata[3*(ii+jj*IMW)+0]=meangroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=meangroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=meangroundtriple[2]; 
			};
		
		// NEW POSITIONS, ASSIGN FIRST batch -- equal number of bright and dark
		for (int ii=0; ii<n1; ii++) // random first
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			imagec[ii]=-1;
		};
		for (int ii=n1; ii<n1+n2; ii++) // then ones with enforced correlations
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
			imagec[ii]=1;
		};

		for (int ii=0; ii<n1+n2; ii++)
		{
			if (imagec[ii] == 1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	if (UPDATE && (loopnum % (int)(2*Stimulus.tau+0.5) == int(Stimulus.tau+0.5)))
	{  
		// GRAY replacement, if necessary -- otherwise leave it
		if (!pt3choose)
		{
			for (int ii=0; ii<n1+n2; ii++)
			{
					imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=meangroundtriple[0];
					imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=meangroundtriple[1];
					imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=meangroundtriple[2]; 
			};
		};

		// PEFORM THE SHIFT and DRAW PT2
		for (int ii=n1; ii<n1+n2; ii++)
		{
			imagex[ii] += xshift; // just the shift...
			imagex[ii] %= xlim;
			imagec[ii] *= flipc;
		};
		for (int ii=0; ii<n1; ii++)
		{
			imagex[ii] = rand() % xlim; // random
			imagey[ii] = rand() % ylim;
			imagec[ii] *= flipc;
		};
		for (int ii=0; ii<n1+n2; ii++)
		{
			if (imagec[ii] == 1)
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			}
			else
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];
			};
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	};
	
	Sfull.x = pt3choose;
	Sfull.y = flipc;

	// DRAW ALL THAT ONTO WEDGES
	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};









void TrackStim::drawRasterRDotRot(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=128;
	static int IMH=64;
	static float imagedata[128*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex[128*64],imagey[128*64];
	int xlim = min(360/Stimulus.spacing,128);
	int ylim = min(10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*xlim*ylim),IMW*IMH); // 
	// bool AFTERWAIT = (tc > Stimulus.tau2);

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<128; ii++)
			for (int jj=0; jj<64; jj++)
			{
				
				imagedata[3*(ii+jj*IMW)+0]=backgroundtriple[0];
				imagedata[3*(ii+jj*IMW)+1]=backgroundtriple[1];
				imagedata[3*(ii+jj*IMW)+2]=backgroundtriple[2]; 
				
			};

		for (int ii=0; ii<numsquares; ii++)
		{
			imagex[ii]=rand() % xlim;
			imagey[ii]=rand() % ylim;
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

		FIRST = FALSE;
	};
	

	if (UPDATE)   // update only on the taus -- if even, just move, otherwise, reassign
	{  
		for (int ii=0; ii<numsquares; ii++)
		{
			if (float(rand())/float(RAND_MAX) < (1/ScreenUpdateRate/2)/Stimulus.tau) // P{change} = dt/tau
			{
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=backgroundtriple[0]; // zero old position
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=backgroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=backgroundtriple[2];

				imagex[ii]=rand() % xlim; // find new position
				imagey[ii]=rand() % ylim;

				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+0]=foregroundtriple[0]; // make new position good
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+1]=foregroundtriple[1];
				imagedata[3*(imagex[ii]+imagey[ii]*IMW)+2]=foregroundtriple[2];
			};
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 128, 64, 0, GL_RGB, GL_FLOAT, imagedata);

	}

	float wedgewidth = 20*PI/180; // width of each wedge
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;

	float added = 0; //Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth + added*PI/180;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);

};


void TrackStim::drawRasterRDotHemispheres(float tc)
{
	static bool FIRST = TRUE; // initializes everything
	static bool UPDATE = TRUE;
	static int IMW=64;
	static int IMH=64;
	static float imagedata_left[64*64*3]={0.0f};  // width come before height!
	static float imagedata_right[64*64*3]={0.0f};
	static float lasttc = 1; // used to tell if need to update...
	static int loopnum = -1;
	static int imagex_left[64*64],imagey_left[64*64],imagex_right[64*64],imagey_right[64*64];
	int xlim = min(360/Stimulus.spacing,64); // just use this for the alignment, but fill in all squares...
	int ylim = min(2*10*Stimulus.arenaheight/(30*Stimulus.spacing*PI/180),64); // radius 30 distance equivalent
	int numsquares = min(int(-1*log(1-Stimulus.density)*64*64),64*64); // same for left and right
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0;
	static float  currold_right=0;

	if (abs(tc - lasttc) < 0.001) // i.e., they are during the same frame... this routine gets called multiple times per frame (like 6)...
		UPDATE = FALSE;
	else
	{
		UPDATE = TRUE;     // go into update routine
		loopnum++;         // update the loop number
		lasttc = tc;       // set this so that it won't update again until tc changes...
		old_left = currold_left;
		old_right = currold_right;
	};

	if (FIRST)   // NEED to know format of pixels in image data array! figure that out and we're golden...
	{
		for (int ii=0; ii<64; ii++)
			for (int jj=0; jj<64; jj++)
			{
				imagedata_left[3*(ii+jj*IMW)+0]=backgroundtriple[0];
				imagedata_left[3*(ii+jj*IMW)+1]=backgroundtriple[1];
				imagedata_left[3*(ii+jj*IMW)+2]=backgroundtriple[2];

				imagedata_right[3*(ii+jj*IMW)+0]=backgroundtriple[0];
				imagedata_right[3*(ii+jj*IMW)+1]=backgroundtriple[1];
				imagedata_right[3*(ii+jj*IMW)+2]=backgroundtriple[2];

			};

		for (int ii=0; ii<numsquares; ii++)
		{
			imagex_left[ii]=rand() % 64;
			imagey_left[ii]=rand() % 64;

			imagex_right[ii]=rand() % 64;
			imagey_right[ii]=rand() % 64;
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // LEFT is 1
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 64, 64, 0, GL_RGB, GL_FLOAT, imagedata_left);

		glBindTexture(GL_TEXTURE_2D, 2); // RIGHT is 2
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 64, 64, 0, GL_RGB, GL_FLOAT, imagedata_right);

		FIRST = FALSE;
	};
	

	if (UPDATE)   // update only on the taus -- if even, just move, otherwise, reassign
	{  
		for (int ii=0; ii<numsquares; ii++)
		{
			if (float(rand())/float(RAND_MAX) < (1/ScreenUpdateRate/2)/Stimulus.tau) // P{change} = dt/tau
			{
				// first do for LEFT
				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+0]=backgroundtriple[0]; // zero old position
				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+1]=backgroundtriple[1];
				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+2]=backgroundtriple[2];

				imagex_left[ii]=rand() % 64; // find new position
				imagey_left[ii]=rand() % 64;

				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+0]=foregroundtriple[0]; // make new position good
				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+1]=foregroundtriple[1];
				imagedata_left[3*(imagex_left[ii]+imagey_left[ii]*IMW)+2]=foregroundtriple[2];

				// then for RIGHT
				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+0]=backgroundtriple[0]; // zero old position
				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+1]=backgroundtriple[1];
				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+2]=backgroundtriple[2];

				imagex_right[ii]=rand() % 64; // find new position
				imagey_right[ii]=rand() % 64;

				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+0]=foregroundtriple[0]; // make new position good
				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+1]=foregroundtriple[1];
				imagedata_right[3*(imagex_right[ii]+imagey_right[ii]*IMW)+2]=foregroundtriple[2];


			};
		};

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1); // rebind the new textures -- LEFT
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 64, 64, 0, GL_RGB, GL_FLOAT, imagedata_left);

		glBindTexture(GL_TEXTURE_2D, 2); // rebind the new textures -- RIGHT
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 64, 64, 0, GL_RGB, GL_FLOAT, imagedata_right);


	};

	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = -(FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*tc*PI/180) +
		Scurr.theta + Stimulus.stimrot.mean*tc*PI/180 +
		old_left; // in radians...; 
	float add_right = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*tc*PI/180) + 
		Scurr.theta + Stimulus.stimrot.mean*tc*PI/180 +
		old_right;
	currold_left = -Stimulus.stimtrans.mean*tc*PI/180 + Stimulus.stimrot.mean*tc*PI/180;
	currold_right = Stimulus.stimtrans.mean*tc*PI/180 + Stimulus.stimrot.mean*tc*PI/180;

	
	float wedgewidth = 20*PI/180; // width of each wedge = 20 degrees
	float per = 180; // period of this sucker, for each half!
	float radius = 30;
	float heightfactor = 10;
	float angle;

	
	// LEFT && RIGHT
	float fracdif;
	for (int pattern = 1; pattern < 3; pattern++) // patterns 1 and 2
	{
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, pattern);
		glBegin(GL_QUADS);

		for (int ii=0; ii<9; ii++) // through each wedge...
		{
			angle = ii*wedgewidth;
			if (pattern == 1) // LEFT case
			{
				fracdif = (ii*20.0f/per) + (add_left/PI);
			}
			else
			{
				fracdif = (ii*20.0f/per) + (add_right/PI);
				angle += PI;
			};

			// go through the slices, but offset things properly
			glTexCoord2f(fracdif*float(xlim)/float(IMW),0.0f);
			glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);
			glTexCoord2f((fracdif+20.0f/per)*float(xlim)/float(IMW),0.0f);
			glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);
			glTexCoord2f((fracdif+20.0f/per)*float(xlim)/float(IMW),float(ylim)/float(IMH)*2);
			glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
			glTexCoord2f(fracdif*float(xlim)/float(IMW),float(ylim)/float(IMH)*2);
			glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);

		}
		glEnd(); // end quads for each one...
		glDisable(GL_TEXTURE_2D);
	}

	

};





void TrackStim::drawBrightCylinder(float tc)
{
	int ii;
	float angle = Stimulus.spacing;  // in degrees
	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	for (ii=0; ii<360/angle; ii++)
	{
		TrackStim::drawWedge(ii*angle + added,angle);
	}
};

void TrackStim::drawCheckerboard(float tc)   // checkerboard with varying phase 
{
	int ii, jj;
	//float tc = TrackStim::queryCurrTime() - epochtimezero;
	float angle = Stimulus.spacing;  // in degrees
	float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	added += tc*Stimulus.stimtrans.mean;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!
	float radius = 30;
	float squareheight=angle*PI/180*radius; // units, not angle -- angle/2 for arbitrary initial reasons -- change this...


	int numSquares = int(360/angle/2)*int(Stimulus.arenaheight*10.0f/squareheight);
	float *pVertex;
	pVertex = vertexArray;

	// cylinder radius is 30, height is 10x arena height; vertical extent is 1/2 of horizontal to get more squares in
	// this is VERY SLOW. use an array, and then glDrawArrays(), as Alex did...
	float locwidth = angle *PI/180;
	float locang;
	
	//glBegin(GL_QUADS);
	for (ii=0; ii<360/angle/2; ii++)
	{
		for (jj=0; jj<(Stimulus.arenaheight*10.0f/squareheight); jj++)
		{
			// AS FAR AS I CAN TELL, DOING IT IN A MATRIX DOESN'T HELP SPEED ANY, OR MUCH <4000 SQUARES IS MAX TO GET 240 Hz DRAW RATE: 3d, float, less than Alex. weird.
			// seems to work down to about 10 degrees = 1026 squares... not much more.
			if ((jj % 2)==1)
				locang = ii*angle*2 + added + Stimulus.stimrot.phase;
			else
				locang = ii*angle*2 + added;
			
			
			locang *= PI/180;
	        *pVertex++ = radius*cos(locang);  // fields of pointer
	        *pVertex++ = radius*sin(locang);
			*pVertex++ = jj*squareheight;  // 

		    *pVertex++ = radius*cos(locang+locwidth);  // fields of pointer
	        *pVertex++ = radius*sin(locang+locwidth);
			*pVertex++ = jj*squareheight;  // z value of 0

	        *pVertex++ = radius*cos(locang+locwidth);  // fields of pointer
	        *pVertex++ = radius*sin(locang+locwidth);
			*pVertex++ = (jj+1)*squareheight;  // 

	        *pVertex++ = radius*cos(locang);  // fields of pointer
	        *pVertex++ = radius*sin(locang);
			*pVertex++ = (jj+1)*squareheight;  // 
			

			//drawWedgeSquare_nobegin(locang,angle,jj*squareheight);
		}
	}
	//glEnd();
    glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertexArray);  // do my stuff this way -- may be faster!
    glDrawArrays(GL_QUADS,0, numSquares*4);
    glDisableClientState(GL_VERTEX_ARRAY);
	glFlush();

};

void TrackStim::drawColorTest()
{
	float dist = 10;
	//static bool coloron = TRUE;

	if (framenumber % 2)
		setForeground();  // wouldn't normally want to change these colors here, but in this routine we do
	else
		setBackgroundColor();

	GLUquadricObj *cyl = gluNewQuadric();
	glPushMatrix();
	gluSphere(cyl,5,4,4);
	glPopMatrix();

	//coloron = !coloron; // switch each frame...

};

void TrackStim::drawSingleBar(float timein)
{
	
	drawWedge(Stimulus.stimrot.mean - Stimulus.spacing/2 + Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per),Stimulus.spacing);

	Sfull.x = Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per);


};



void TrackStim::drawDoubleBars(float timein)
{
	
	drawWedge(Stimulus.stimrot.mean - Stimulus.spacing/2 + Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per),Stimulus.spacing);

	Sfull.x = Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per);

	drawWedge(-Stimulus.stimrot.mean - Stimulus.spacing/2 + Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per),Stimulus.spacing);


};

void TrackStim::drawDoubleBars_int(float timein)
{
	
	float spacing=Stimulus.spacing2;
	if (spacing < 0.00001)
		spacing = 1;

	float ang=spacing*(int)((Stimulus.stimrot.mean - Stimulus.spacing/2 + Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per))/spacing);
	drawWedge(ang,Stimulus.spacing);

	//Sfull.x = Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per);
	Sfull.x = ang;
	
	ang=spacing*(int)((-Stimulus.stimrot.mean - Stimulus.spacing/2 + Stimulus.stimrot.amp*sin(2*PI*timein/Stimulus.stimrot.per))/spacing);
	drawWedge(ang,Stimulus.spacing);


};

void TrackStim::drawDoubleBars_rand(float timein)
{
	
	drawWedge(Stimulus.stimrot.mean - Stimulus.spacing/2 + r1.readLastValue()*Stimulus.stimrot.amp,Stimulus.spacing);

	Sfull.x = r1.readLastValue()*Stimulus.stimrot.amp;
	Sfull.y = r1.readLastValue()*Stimulus.stimrot.amp;

	drawWedge(-Stimulus.stimrot.mean - Stimulus.spacing/2 + r1.readLastValue()*Stimulus.stimrot.amp,Stimulus.spacing);


};

void TrackStim::drawDoubleBars_rand_indep(float timein)
{
	
	drawWedge(Stimulus.stimrot.mean - Stimulus.spacing/2 + r1.readLastValue()*Stimulus.stimrot.amp,Stimulus.spacing);

	Sfull.x = r1.readLastValue()*Stimulus.stimrot.amp;
	Sfull.y = r2.readLastValue()*Stimulus.stimrot.amp;

	drawWedge(-Stimulus.stimrot.mean - Stimulus.spacing/2 + r2.readLastValue()*Stimulus.stimrot.amp,Stimulus.spacing);


};

void TrackStim::drawSymmetricBars(float timein)
{

	float wid = Stimulus.spacing;
	float init_off = Stimulus.spacing2;
	float speed = Stimulus.stimrot.mean;
	float pos = Stimulus.stimrot.amp;
	float offset;
	
	if (timein<Stimulus.tau)
	{
		offset = -init_off*(speed/abs(speed)); // sign corrected
	}
	else
	{
		offset = -init_off*(speed/abs(speed)) + speed * timein;
	};
	
	// draw the two bars...
	drawWedge(pos + offset - wid/2,wid);
	drawWedge((360 - pos) + offset - wid/2,wid);


};

void TrackStim::updateReichardtIntensity()
{
	static int counter = 0;
	if (counter==0)
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
		r1.setdt(1/240.0f); // this must get called twice per frame now...
		r1.setTau(Stimulus.tau);
		r2.setdt(1/240.0f);
		r2.setTau(Stimulus.tau);
	};
	r1.Rrandnt();
	r2.Rrandnt();

	counter++;

};

void TrackStim::updateLumIntensity()
{
	static int counter = 0;
	if (counter==0)
	{
		r4.resetSeed(0);
		r4.setdt(1/240.0f); // this must get called twice per frame now...
		r4.setTau(Stimulus.tau);
	};
	r4.Rrandnt();

	counter++;

};

void TrackStim::updateReichardtIntensity_binary()
{
	static int counter = 0;
	if (counter == 0) // tau is in frames here, not ms...
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
	};
	if ((counter % (int)(Stimulus.tau+0.1)) == 0)
	{
		r1.Rrand();
		r2.Rrand();
	};
	r1.setCount(counter);
	r2.setCount(counter);

	counter++;  // updates every BG, at 240 Hz

};

void TrackStim::updateReichardtIntensity_binary_isolating()
{
	static int counter = 0;
	if (counter == 0) // tau is in frames here, not ms...
	{
		r1.resetSeed(0);
	};
	if ((counter % (int)(Stimulus.tau+0.1)) == 0)
	{
		r1.Rrand();
	};
	r1.setCount(counter);

	counter++;  // updates every BG, at 240 Hz

};

void TrackStim::updateReichardtIntensity_countonly()
{
	static int counter = 0;
	if (counter == 0) // tau is in frames here, not ms...
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
	};
	r1.setCount(counter);
	r2.setCount(counter);

	counter++;  // updates every BG, at 240 Hz

};

void TrackStim::updateReichardtIntensity_discrete_gaussian()
{
	static int counter = 0;
	if (counter == 0) // tau is in frames here, not ms...
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
	};
	if ((counter % (int)(Stimulus.tau+0.1)) == 0)
	{
		r1.Rrandn();
		r2.Rrandn();
	};
	r1.setCount(counter);
	r2.setCount(counter);

	counter++;  // updates every BG, at 240 Hz

};

void TrackStim::updateReichardtIntensity_binary_3bars()
{
	static int counter = 0;
	if (counter == 0) // tau is in frames here, not ms...
	{
		r1.resetSeed(0);
		r2.resetSeed(1);
		r3.resetSeed(2);
	};
	if ((counter % (int)(Stimulus.tau+0.1)) == 0)
	{
		r1.Rrand();
		r2.Rrand();
		r3.Rrand();
	};
	r1.setCount(counter);
	r2.setCount(counter);
	r3.setCount(counter);

	counter++;  // updates every BG, at 240 Hz

};

void TrackStim::updateReichCorr()
{
	static int counter = 0;
	int rem = counter % int(Stimulus.tau*2+0.01); // tau in frames!

	if (rem == 0)
		RCorr.randomizepixels();
	if (rem == Stimulus.tau)
		RCorr.updatepixels();

	counter++; // updates at 240 Hz
};

void TrackStim::updateReichCorr_everytime()
{

	int N = 360/int(Stimulus.spacing+0.01);
	// deal with RCorr
	RCorr.setPwhite(Stimulus.density); // start at 0.5...
	RCorr.setBC(bool(int(Stimulus.stimtrans.phase+0.01))); // 0 or 1 on stimtrans phase
	RCorr.setnumberofpixels(N);
	RCorr.setupdatetype(int(Stimulus.stimrot.phase+0.01));
	RCorr.setDir(int(Stimulus.rot.phase+0.01));
	RCorr.setparity(int(Stimulus.rot.per+1.01)-1);
	RCorr.setparityNUM(int(Stimulus.trans.per+0.01));
	RCorr.setBackFrames(int(Stimulus.stimrot.per+0.01));

	static int counter = 0;
	int rem = counter % int(Stimulus.tau+0.01); // tau in frames!

	if (rem == 0)
		RCorr.updatepixels();

	counter++; // updates at 240 Hz
};

void TrackStim::updateReichCorrGrid_everytime()
{
	static int counter = 0;
	int rem = counter % int(Stimulus.tau+0.01); // tau in frames!

	if (rem == 0)
		for (int ii=0; ii<8; ii++)
			RCorr2D[ii].updatepixels();

	counter++; // updates at 240 Hz
};

void TrackStim::drawReichCorrFullBars()
{
	float pixv[360];
	float mlum=Stimulus.lum;
	int N = 360/int(Stimulus.spacing+0.01);
	// deal with RCorr in update section, earlier

	for (int ii=0;ii<N;ii++)
		pixv[ii]=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(ii))-0.5f)*2);

	//Sfull.x=pixv[1]*100; Sfull.y=pixv[2]*100;  // store values here? perhaps not wisest, but should work.
	Sfull.x = 0;
	for (int ii=0;ii<12; ii++)
	{
		Sfull.x *= 2; // convert to base 2 to get it out...
		Sfull.x += RCorr.getpixelvalue(ii);
	};

	/*if (int(Stimulus.rot.phase+0.01) == 1)
	{
		for (int ii=0;ii<N;ii++)
			pixv[ii]=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(ii))-0.5f)*2);
	}
	else // going the opposite direction
	{
		for (int ii=0;ii<N;ii++)
			pixv[ii]=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(N-1-ii))-0.5f)*2);
	};*/

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

		// second, draw one set of panels
	for (int ii=0; ii<N; ii++)
	{
		glColor3f(0.0f,getDLPColor(pixv[ii],'G'),getDLPColor(pixv[ii],'B'));
		drawWedge(ii*Stimulus.spacing, Stimulus.spacing);
	}
};

void TrackStim::drawReichCorrFullGrid()
{
	static int IMW=256;
	static int IMH=8;
	static float imagedata[256*8*3]={1.0f};

	float mlum=Stimulus.lum;
	int N = 360/int(Stimulus.spacing+0.01);
	// deal with RCorr2D -- update it every time...
	for (int ii=0; ii<8; ii++)
	{
		RCorr2D[ii].setPwhite(Stimulus.density); // start at 0.5...
		RCorr2D[ii].setBC(bool(int(Stimulus.stimtrans.phase+0.01))); // 0 or 1 on stimtrans phase
		RCorr2D[ii].setnumberofpixels(N);
		RCorr2D[ii].setupdatetype(int(Stimulus.stimrot.phase+0.01));
		RCorr2D[ii].setDir(int(Stimulus.rot.phase+0.01));
		RCorr2D[ii].setparity(int(Stimulus.rot.per+1.01)-1);
		RCorr2D[ii].setparityNUM(int(Stimulus.trans.per+0.01));
	};

	// output a few pixels for checking purposes -- base 2 in x column
	Sfull.x = 0;
	for (int ii=0;ii<12; ii++)
	{
		Sfull.x *= 2; // convert to base 2 to get it out...
		Sfull.x += RCorr2D[0].getpixelvalue(ii);
	};

	if (RCorr2D[0].PIXUPDATED)
	{
		// update all pixels...
		for (int ii=0; ii<N; ii++)
			for (int jj=0; jj<8; jj++)
			{
				if (RCorr2D[jj].getpixelvalue(ii)==0)
				{
					imagedata[3*(ii+jj*IMW)+0]=backgroundtriple[0];
					imagedata[3*(ii+jj*IMW)+1]=backgroundtriple[1];
					imagedata[3*(ii+jj*IMW)+2]=backgroundtriple[2]; 
				}
				else
				{
					imagedata[3*(ii+jj*IMW)+0]=foregroundtriple[0];
					imagedata[3*(ii+jj*IMW)+1]=foregroundtriple[1];
					imagedata[3*(ii+jj*IMW)+2]=foregroundtriple[2]; 
				};
			};
		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, IMW, IMH, 0, GL_RGB, GL_FLOAT, imagedata);

		RCorr2D[0].PIXUPDATED=FALSE; // only do this ONCE!
	}


	// DO THE DRAWING WITH THE PANELS
	float wedgewidth = 20*PI/180; // width of each wedge = 20 degrees
	float per = 360.0f; // period of this sucker
	float radius = 30;
	float heightfactor = 10;
	float angle;
	int xlim = N; // keeps track of how many squares to show in the texture, x axis
	float ylim = Stimulus.arenaheight*heightfactor/radius; // arena height in closest-approach-radians
	ylim = ylim/(Stimulus.spacing*PI/180); // number of squares in the y direction

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
	for (int ii=0; ii<18; ii++) // through each wedge...
	{
		angle = ii*wedgewidth;
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle),radius*sin(angle), 0.0f);

		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),0.0f);
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), 0.0f);

		glTexCoord2f((ii+1)*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle+wedgewidth),radius*sin(angle+wedgewidth), Stimulus.arenaheight*heightfactor);
		
		glTexCoord2f(ii*20.0f/per*float(xlim)/float(IMW),float(ylim)/float(IMH));
		glVertex3f(radius*cos(angle),radius*sin(angle), Stimulus.arenaheight*heightfactor);
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
};

void TrackStim::drawReichCorr2Bars_nocenter()
{
	// deal with RCorr
	RCorr.setPwhite(Stimulus.density); // start at 0.5...
	RCorr.setBC(bool(int(Stimulus.stimtrans.phase+0.01))); // 0 or 1 on stimtrans phase
	RCorr.setnumberofpixels(2);
	RCorr.setupdatetype(int(Stimulus.stimrot.phase+0.01));
	RCorr.setDir(int(Stimulus.rot.phase+0.01));


	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1, lum2;

	lum1=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(0))-0.5f)*2);
	lum2=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(1))-0.5f)*2);


	/*if (int(Stimulus.rot.phase+0.01) == 1)
	{
		lum1=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(0))-0.5f)*2);
		lum2=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(1))-0.5f)*2);
	}
	else // going the opposite direction
	{
		lum1=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(1))-0.5f)*2);
		lum2=Stimulus.lum*(1+Stimulus.contrast*(float(RCorr.getpixelvalue(0))-0.5f)*2);
	};*/


	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};


void TrackStim::drawReichardtSlits()
{

	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	drawWedge(90.0f-Stimulus.spacing, Stimulus.spacing);
	drawWedge(270.0f-Stimulus.spacing, Stimulus.spacing);

	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	drawWedge(90.0f, Stimulus.spacing);
	drawWedge(270.0f, Stimulus.spacing);
};

void TrackStim::drawReichardtSlits_multiple()
{
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};

void TrackStim::drawReichardtSlits_multiple_nocenter()
{
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};


void TrackStim::drawReichardtSlits_multiple_nocenter_1static(float tc)
{
	static bool LEFT = TRUE;
	static bool UP = TRUE;
	static bool RESET = TRUE;
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());

	if ((mod(tc,Stimulus.tau2) < 0.01) & (RESET))
	{
		RESET = FALSE;
		LEFT = !LEFT;
		UP = rand() % 2;
	};
	if (mod(tc,Stimulus.tau2)>0.1)
		RESET=TRUE;

	if (LEFT)
	{
		if (UP)
			lum1=1;
		else
			lum1=0;
	}
	else
	{
		if (UP)
			lum2=1;
		else
			lum2=0;
	};



	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if ((ii*locspace>15) & (ii*locspace<345))
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if ((ii*locspace>15) & (ii*locspace<345))
				drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};

void TrackStim::drawReichardtSlits_deterministic_steps(float tc)
{
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lumlow = Stimulus.lum*(1-Stimulus.contrast);
	float lumhigh = Stimulus.lum*(1+Stimulus.contrast);
	float timeoff = 0.02;

	static float lum1 = lumlow;
	static float lum2 = lumlow;

	static bool CHANGE1 = FALSE;
	static bool CHANGE2 = FALSE;
	
	if (tc < 2*timeoff)
	{
		CHANGE1 = FALSE;
		CHANGE2 = FALSE;
	};

	if ((Stimulus.tau > 0.001))
	{
		if ((tc > Stimulus.duration - Stimulus.tau - timeoff) && (!CHANGE1))
		{
			if (lum1 < mlum)
				lum1 = lumhigh;
			else
				lum1 = lumlow;
			CHANGE1 = TRUE; // these make it run only once
		};

		if ((tc > Stimulus.duration - timeoff) && (!CHANGE2))
		{
			if (lum2 < mlum)
				lum2 = lumhigh;
			else
				lum2 = lumlow;
			CHANGE2 = TRUE;
		};

	};

	if (Stimulus.tau < -0.001)
	{
		if ((tc > Stimulus.duration + Stimulus.tau - timeoff) && (!CHANGE2))
		{
			if (lum2 < mlum)
				lum2 = lumhigh;
			else
				lum2 = lumlow;
			CHANGE2 = TRUE;
		};

		if ((tc > Stimulus.duration - timeoff) && (!CHANGE1))
		{
			if (lum1 < mlum)
				lum1 = lumhigh;
			else
				lum1 = lumlow;
			CHANGE1 = TRUE;
		};
	};

	if (abs(Stimulus.tau) < 0.001)
	{
		if ((tc > Stimulus.duration - timeoff) && (!CHANGE2))
		{
			if (lum2 < mlum)
				lum2 = lumhigh;
			else
				lum2 = lumlow;
			CHANGE2 = TRUE;
		};

		if ((tc > Stimulus.duration - timeoff) && (!CHANGE1))
		{
			if (lum1 < mlum)
				lum1 = lumhigh;
			else
				lum1 = lumlow;
			CHANGE1 = TRUE;
		};
	};



	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	//if (CHANGE1 && CHANGE2)
	//	Sfull.theta = 100*mlum;

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};


void TrackStim::drawReichardtSlits_multiple_nocenter_balanced()
{
	// c1 = 0.5 - c2
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1-Stimulus.contrast*r1.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	lum1=min(1.0,lum1); lum2=min(1.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};




void TrackStim::drawReichardtSlits_multiple_oneeye()
{
	bool USEPERSP = FALSE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	float lum12=max(0.0,(lum1+lum2)/2); // holds the mean luminance for the other sides...
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace>=180) // right eye
					drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
			}
			else
			{
				if (ii*locspace<180) // right eye
					drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
			};
			
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace>=180) // right eye
					drawWedge(ii*locspace, Stimulus.spacing);
			}
			else
			{
				if (ii*locspace<180) // right eye
					drawWedge(ii*locspace, Stimulus.spacing);
			}
		}
		
		// 090103 - WAS MESSED UP! OLDER THAN THIS, MONOCULAR STUFF DIDN'T WORK PROPERLY, PUTTING LUM2 ON THE OPPOSITE EYE. SUCKED!
		glColor3f(0.0f,getDLPColor(lum12,'G'),getDLPColor(lum12,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace<180) // right eye
					drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
			}
			else
			{
				if (ii*locspace>=180) // left eye
					drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
			}
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace>=180)
					drawPerspectiveStrips(ii*locspace-Stimulus.spacing,Stimulus.spacing,0, TRUE);
			}
			else
			{
				if (ii*locspace<180)
					drawPerspectiveStrips(ii*locspace-Stimulus.spacing,Stimulus.spacing,0, TRUE);
			};
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace>=180)
					drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
			}
			else
			{
				if (ii*locspace<180)
					drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
			}

		}

		glColor3f(0.0f,getDLPColor(lum12,'G'),getDLPColor(lum12,'B'));
		for (int ii=1; ii<360/locspace; ii++)
		{
			if (abs(Stimulus.rot.phase - 1)<0.01)
			{
				if (ii*locspace<180) // right eye
					drawPerspectiveStrips(ii*locspace-Stimulus.spacing, Stimulus.spacing*2,0,FALSE);
			}
			else
			{
				if (ii*locspace>=180) // left eye
					drawPerspectiveStrips(ii*locspace-Stimulus.spacing, Stimulus.spacing*2,0,FALSE);
			}
		}
	};
};

void TrackStim::drawReichardtSlits_multiple_dots(float tc)
{
	static bool FIRSTUSE = TRUE;
	bool RESETDOTS = FALSE;
	bool USEDOTS = TRUE;
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (tc < 0.002)  // reset dots on first call each time...
		if (FIRSTUSE)
			RESETDOTS = TRUE;
		else
			FIRSTUSE = FALSE;
	else
		FIRSTUSE = TRUE; 


	if (!USEDOTS)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge(ii*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawReichardtDots_nopersp(ii*locspace,Stimulus.spacing,0, TRUE, RESETDOTS);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawReichardtDots_nopersp(ii*locspace,Stimulus.spacing,0, FALSE, RESETDOTS);
		}
	};
};


void TrackStim::drawReichardtSlits_multiple_widespace()
{
	bool USEPERSP = FALSE;
	static float locspace = 45.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (!USEPERSP)
	{
		// second, draw one set of panels
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge((0.5 + ii)*locspace-Stimulus.spacing, Stimulus.spacing);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawWedge((ii+0.5)*locspace, Stimulus.spacing);
		}
	}
	else // draw them with perspective, always 5 deg height, same width with height increase...
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, TRUE);
		}


		// third, draw second set of panels
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		for (int ii=0; ii<360/locspace; ii++)
		{
			drawPerspectiveStrips(ii*locspace,Stimulus.spacing,0, FALSE);
		}
	};
};

void TrackStim::drawPerspectiveStrips(float angle, float width, float gap, bool LEFT)
{
	float radius = 30.0f;
	static float rads5=5*PI/180;
	float h,dh,leftang,rightang;
	float dh0=radius*tan(rads5);  // basic 5 degrees, largest...
	angle *= PI/180;
	width *= PI/180;
	gap *= PI/180;

	glBegin(GL_QUADS);
	for (int ii=0; ii<12; ii++)
	{
		h=radius*tan(rads5*ii); // top/-bottom of all boxes... 
		dh=radius*tan(rads5*(ii+1))-h;
		if (LEFT)
		{
			leftang = angle - (width + gap/2)/cos(rads5*ii);
			rightang = angle - (gap/2)/cos(rads5*ii);
		}
		else
		{
			leftang = angle + (gap/2)/cos(rads5*ii);
			rightang = angle + (width + gap/2)/cos(rads5*ii);
		};
		
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   + h); // get at right level...
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   + h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   +h+dh);
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   +h+dh);

		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   -h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   -h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   -h-dh);
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   -h-dh);


	};
	glEnd();
};


void TrackStim::drawReichardtDots_nopersp(float angle, float width, float gap, bool LEFT, bool FIRST)
{
//	static bool FIRST1 = TRUE;
	float radius = 30.0f;
	static float rads5=5*PI/180;
	static float offsets[12];
	float h,dh,leftang,rightang;
	float dh0=radius*tan(rads5);  // basic 5 degrees, largest...
	angle *= PI/180;
	width *= PI/180;
	gap *= PI/180;

	if (FIRST)
	{
		// FIRST1 = FALSE;
		for (int ii=0; ii<12; ii++)
			offsets[ii]=(rand() % 360) * PI/180;
	}


	glBegin(GL_QUADS);
	for (int ii=0; ii<12; ii++)
	{
		h=radius*tan(rads5*ii); // top/-bottom of all boxes... 
		dh=radius*tan(rads5*(ii+1))-h;
		if (LEFT)
		{
			leftang = angle - (width + gap/2) + offsets[ii];
			rightang = angle - (gap/2) + offsets[ii];
		}
		else
		{
			leftang = angle + (gap/2) + offsets[ii];
			rightang = angle + (width + gap/2) + offsets[ii];
		};
		
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   + h); // get at right level...
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   + h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   +h+dh);
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   +h+dh);

		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   -h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   -h);
		glVertex3f(radius*cos(rightang),radius*sin(rightang),Stimulus.arenaheight*5   -h-dh);
		glVertex3f(radius*cos(leftang),radius*sin(leftang),Stimulus.arenaheight*5   -h-dh);


	};
	glEnd();
};


void TrackStim::drawReichardtSlits_multiple_spacing()
{
	static float locspace = 30.0f; // offset between bar pairs...
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace-Stimulus.spacing - Stimulus.spacing2/2, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace + Stimulus.spacing2/2, Stimulus.spacing);
	}
};


void TrackStim::drawReichardtSlits_multiple_binary()
{
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter_isolating()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=mlum;
	static float lum2=mlum;
	int count = r1.readCount();
	static int waittime = int(Stimulus.tau + 0.01);
	static bool RIGHT = FALSE;
	
	////// BEGIN PHASE -- choose direction, set bars to mlum
	if ((count % (3*waittime)) == 0)
	{
		RIGHT = (r1.readLastValue()>0.5);
		lum1 = mlum;
		lum2 = mlum;
	};

	////// FIRST BAR CHOSEN
	if ((count % (3*waittime)) == waittime) // first update
	{
		if (RIGHT)
		{
			if (r1.readLastValue()>0.5)
				lum1=Stimulus.lum*(1+Stimulus.contrast);
			else
				lum1=Stimulus.lum*(1-Stimulus.contrast);
			lum2=mlum;
		}
		else
		{
			if (r1.readLastValue()>0.5)
				lum2=Stimulus.lum*(1+Stimulus.contrast);
			else
				lum2=Stimulus.lum*(1-Stimulus.contrast);
			lum1=mlum;
		};
	};

	////// SECOND BAR CHOSEN
	if ((count % (3*waittime)) == 2*waittime) // first update
	{
		if (RIGHT)
		{
			if (r1.readLastValue()>0.5)
				lum2=Stimulus.lum*(1+Stimulus.contrast);
			else
				lum2=Stimulus.lum*(1-Stimulus.contrast);
			lum1=mlum;
		}
		else
		{
			if (r1.readLastValue()>0.5)
				lum1=Stimulus.lum*(1+Stimulus.contrast);
			else
				lum1=Stimulus.lum*(1-Stimulus.contrast);
			lum2=mlum;
		};
	};


	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}

	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};


void TrackStim::drawReichardtSlits_original_gangster(float tc)
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_stagger(float tc)
{
	static int count = 0;
	static bool FIRST = 1;
	static int phase = 0;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	if ((tc<0.08) & (FIRST)) // latch this to change only on the first call each round!
	{
		count ++;
		FIRST = 0;
	};
	if (tc>0.1) // retrigger FIRST
	{
		FIRST = 1;
	};

	phase = 10 * (count % 2); // random phase in degrees


	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};


	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		drawWedge(phase + ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		drawWedge(phase + ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_TRANSLATE_original_gangster_stagger(float tc)
{
	static int count = 0;
	static bool FIRST = 1;
	static int phase = 0;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	if ((tc<0.08) & (FIRST)) // latch this to change only on the first call each round!
	{
		count ++;
		FIRST = 0;
	};
	if (tc>0.1) // retrigger FIRST
	{
		FIRST = 1;
	};

	phase = locspace/2 * (count % 2); // random phase in degrees


	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};


	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<int(180.1/locspace); ii++)
	{
		drawWedge_mirror(phase + ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<int(180.1/locspace); ii++)
	{
		drawWedge_mirror(phase + ii*locspace, Stimulus.spacing);
	}
};




void TrackStim::drawReichardtSlits_original_gangster_randomphase(float tc)
{
	static bool FIRST = 1;
	static int phase = 0;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	if ((tc<0.05) & (FIRST)) // latch this to change only on the first call each round!
	{
		FIRST = 0;
		phase = rand() % 360; // random phase in degrees
	};
	if (tc>0.1) // retrigger FIRST
	{
		FIRST = 1;
	};

	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		drawWedge(phase + ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		drawWedge(phase + ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawGrayWithVerticalNoise(float tc)
{
	static bool FIRST = 1;
	static int phase = 0;
	static int count =0;
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static bool screenvals[100];
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	if ((tc<0.05) & (FIRST)) // latch this to change only on the first call each round!
	{
		FIRST = 0;
		phase = rand() % 360; // random phase in degrees
		count=0;
	};
	if (tc>0.1) // retrigger FIRST
	{
		FIRST = 1;
	};

	if ((tc>2*waittau) & (tc<=3*waittau))
	{
		lum1=Stimulus.lum*(1+Stimulus.contrast);
		lum2=Stimulus.lum*(1-Stimulus.contrast);
	}
	else
	{
		lum1=Stimulus.lum;
		lum2=Stimulus.lum;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	if ((count % int(3*Stimulus.tau+0.01))==0) // kluge-y here -- 3 renders per frame, by my count
	{
		for (int ii=0;ii<100;ii++)
		{
			screenvals[ii]=(rand() % 2); // random choosing
		};
	};
	count++;

	Sfull.x=lum1*100*screenvals[0]; 
	Sfull.y=lum1*100*screenvals[1];  // store values here? perhaps not wisest, but should work.

	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	for (int ii=0; ii<int(360.1/Stimulus.spacing); ii++)
	{
		if (screenvals[ii])
		{
			glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		}
		else
		{
			glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		};
		drawWedge(phase + ii*Stimulus.spacing, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_offcenter(float tc)
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;
	float direction = 1;

	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
		direction = 1;
	}
	else
	{
		lum2=temp1; lum1=temp2;
		direction = -1;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, direction * Stimulus.spacing);
	}

	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace + direction * Stimulus.spacing, direction * Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_4pt(float tc)
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;
	// note 101119 -- this is wrong -- 2 is 0-1 -- this is corrected in all matlab coding...

	float temp1,temp2,temp0;
	// temp1 is second bar to come on, always
	switch (stimchoose)
	{
		case 0:
			temp0=Stimulus.lum*(1-Stimulus.contrast);
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp0=Stimulus.lum*(1+Stimulus.contrast);
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp0=Stimulus.lum*(1-Stimulus.contrast);
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp0=Stimulus.lum*(1+Stimulus.contrast);
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
	{
	//	temp1 = mlum; // temp1 always holds first value...
		temp1 = temp0; // opposite contrast of step... 4 pts, biatch. should just work.
	};

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_6pt(float tc)
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;
	// note 101119 -- this is wrong -- 2 is 0-1 -- this is corrected in all matlab coding...

	float temp1,temp2,temp0;
	// temp1 is second bar to come on, always
	switch (stimchoose)
	{
		case 0:
			temp0=Stimulus.lum*(1-Stimulus.contrast);
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp0=Stimulus.lum*(1+Stimulus.contrast);
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
	{
	//	temp1 = mlum; // temp1 always holds first value...
		temp1 = temp0; // opposite contrast of step... 4 pts, biatch. should just work.
		temp2 = temp0;
	}
	else {
		if (tc <= 2*waittau)
		{
			temp1 = temp0;
		};
	};

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_6ptlong(float tc)
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;
	// note 101119 -- this is wrong -- 2 is 0-1 -- this is corrected in all matlab coding...

	float temp1,temp2,temp0;
	// temp1 is second bar to come on, always
	switch (stimchoose)
	{
		case 0:
			temp0=Stimulus.lum*(1-Stimulus.contrast);
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp0=Stimulus.lum*(1+Stimulus.contrast);
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
	};
	if (tc <= 5*waittau)
	{
	//	temp1 = mlum; // temp1 always holds first value...
		temp1 = temp0; // opposite contrast of step... 4 pts, biatch. should just work.
		temp2 = temp0;
	}
	else {
		if (tc <= 6*waittau)
		{
			temp1 = temp0;
		};
	};

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_original_gangster_spaced(float tc)
{
	static float locspace = 30.0f;
	float interspacing = Stimulus.spacing2;
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	float waittau = Stimulus.tau2; // how long to wait before second slit appears
	bool ORDER = Stimulus.rot.per; // 1 is clockwise?
	int stimchoose = Stimulus.rot.phase; // 0=1-1; 1=0-0; 2=1-0; 3=0-1;

	float temp1,temp2;
	switch (stimchoose)
	{
		case 0:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
		case 1:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 2:
			temp1=Stimulus.lum*(1+Stimulus.contrast);
			temp2=Stimulus.lum*(1-Stimulus.contrast);
			break;
		case 3:
			temp1=Stimulus.lum*(1-Stimulus.contrast);
			temp2=Stimulus.lum*(1+Stimulus.contrast);
			break;
	};
	if (tc <= waittau)
		temp1 = mlum; // simple as that...

	if (ORDER)
	{
		lum1=temp1; lum2=temp2;
	}
	else
	{
		lum2=temp1; lum1=temp2;
	};

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100; 
		Sfull.y = lum2*100; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	// returned to original, 100615
	for (int ii=0; ii<int(360.1/locspace); ii++)
	// 100614 -- CHANGED TO +/- 15 DEGREES!
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace-Stimulus.spacing-interspacing/2, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<int(360.1/locspace); ii++)
	//for (int ii=int(15.1/locspace); ii<int(345.1/locspace)+1; ii++)
	{
		if ((ii*locspace>=15) && (ii*locspace<=345))
			drawWedge(ii*locspace + interspacing/2, Stimulus.spacing);
	}
};




void TrackStim::drawReichardtSlits_multiple_binary_nocenter_no45()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=int(30.1/locspace); ii<int(330.1/locspace)+1; ii++)
	{
		if ((abs(ii*locspace - 45)>Stimulus.spacing) && (abs(ii*locspace - 315)>Stimulus.spacing))
			drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=int(30.1/locspace); ii<int(330.1/locspace)+1; ii++)
	{
		if ((abs(ii*locspace - 45)>Stimulus.spacing) && (abs(ii*locspace - 315)>Stimulus.spacing))
			drawWedge(ii*locspace, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter_1eye()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0.001)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	float lum12=max(0.0,(lum1+lum2)/2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace>=180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}
		else
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		};
		
	};


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace>=180) // right eye
				drawWedge(ii*locspace, Stimulus.spacing);
		}
		else
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace, Stimulus.spacing);
		};
	};
	
	glColor3f(0.0f,getDLPColor(lum12,'G'),getDLPColor(lum12,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
		}
		else
		{
			if (ii*locspace>=180) // left eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
		};
	};
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter_1eye_quiet()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0.001)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	//float lum12=max(0.0,(lum1+lum2)/2);
	float lum12=mlum;

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace>=180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		}
		else
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		};
		
	};


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace>=180) // right eye
				drawWedge(ii*locspace, Stimulus.spacing);
		}
		else
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace, Stimulus.spacing);
		};
	};
	
	glColor3f(0.0f,getDLPColor(lum12,'G'),getDLPColor(lum12,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01)
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
		}
		else
		{
			if (ii*locspace>=180) // left eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing*2);
		};
	};
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter_1eye_equalspace()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0.001)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	//float lum12=max(0.0,(lum1+lum2)/2);
	//float lum12=mlum; // NO LONGER USE MLUM

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01) // left eye STIM
		{
			if (ii*locspace>=180) // left eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
			else // right eye
				drawWedge(ii*locspace-Stimulus.spacing/2-locspace/4, Stimulus.spacing);
		}
		else // REVERSE
		{
			if (ii*locspace<180) // right eye
				drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
			else // left eye
				drawWedge(ii*locspace-Stimulus.spacing/2-locspace/4, Stimulus.spacing);
		};
		
	};


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		if (abs(Stimulus.rot.phase - 1)<0.01) // left eye STIM
		{
			if (ii*locspace>=180) // left eye
				drawWedge(ii*locspace, Stimulus.spacing);
			else // right eye
				drawWedge(ii*locspace-Stimulus.spacing/2+locspace/4, Stimulus.spacing);
		}
		else // REVERSE
		{
			if (ii*locspace<180) // left eye
				drawWedge(ii*locspace, Stimulus.spacing);
			else
				drawWedge(ii*locspace-Stimulus.spacing/2+locspace/4, Stimulus.spacing);
		};
	};
	
};

void TrackStim::drawReichardtSlits_multiple_binary_nocenter_spacing()
{
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace+Stimulus.spacing2/2, Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_multiple_binary_4corners()
{
	static float locspace = 30.0f;
	if (Stimulus.spacing2 > 0)
		locspace = Stimulus.spacing2; // 091105
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	drawWedge(60-Stimulus.spacing, Stimulus.spacing);
	drawWedge(240-Stimulus.spacing, Stimulus.spacing);

	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	drawWedge(120, Stimulus.spacing);
	drawWedge(300, Stimulus.spacing);

};

void TrackStim::drawReichardtSlits_multiple_binary_3bars()
{
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lum3=0;
	static float lumold1=0;
	static float lumold2=0;
	static float lumold3=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
		lumold3 = lum3;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);

	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	if (r3.readLastValue()>0.5)
		lum3=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum3=Stimulus.lum*(1-Stimulus.contrast);


	lum1=max(0.0,lum1); lum2=max(0.0,lum2); lum3=max(0.0,lum3);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; 
		Sfull.theta = lum3*100 + lumold3*1000;
		// therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (black, then white), 1000 (white, then black), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
		Sfull.theta=lum3*100;
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace - Stimulus.spacing*1.5, Stimulus.spacing);
	};

	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace - Stimulus.spacing*0.5, Stimulus.spacing);
	};

	// fourth, draw third set of panels
	glColor3f(0.0f,getDLPColor(lum3,'G'),getDLPColor(lum3,'B'));
	for (int ii=0; ii<360/locspace; ii++)
	{
		drawWedge(ii*locspace + Stimulus.spacing*0.5, Stimulus.spacing);
	};
};


void TrackStim::drawReichardtSlits_multiple_translate_binary()
{
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	static float lum1=0;
	static float lum2=0;
	static float lumold1=0;
	static float lumold2=0;
	static int countold=0;

	// STORE OLD LUM VALUES WHEN UPDATED, always know last value...
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() != countold) || (r1.readCount() == 0)))
	{ // doing the 1 frame thing, frame updated or first frame
		countold = r1.readCount(); // execute this only once per frame...
		lumold1 = lum1; // store old lum values
		lumold2 = lum2;
	};

	// DEAL WITH FIRST INTENSITY
	if (r1.readLastValue()>0.5)
		lum1=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum1=Stimulus.lum*(1-Stimulus.contrast);
	// DEAL WITH SECOND INTENSITY
	if (r2.readLastValue()>0.5)
		lum2=Stimulus.lum*(1+Stimulus.contrast);
	else
		lum2=Stimulus.lum*(1-Stimulus.contrast);

	lum1=max(0.0,lum1); lum2=max(0.0,lum2);

	// STORE THE INTENSITIES IN POSITION OF STIM FILE
	if ((int(Stimulus.tau+0.01) == 1) && ((r1.readCount() + 1) % 2 == 0)) 
	{ // this is case of tau = 1, alternating every frame, so each location must encode two intensities
		Sfull.x = lum1*100 + lumold1*1000; 
		Sfull.y = lum2*100 + lumold2*1000; // therefore, first of two frames is as usual, second is times 1000, not 100
		// possibilities: 0 (both black), 100 (white, then black), 1000 (black, then white), 1100 (white white)
	}
	else
	{
		Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.
	};

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<180/locspace; ii++)  // translation undefined straight ahead -- kill it!
	{
		drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		drawWedge(ii*-locspace,Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<180/locspace; ii++)
	{
		drawWedge(ii*locspace, Stimulus.spacing);
		drawWedge(ii*-locspace - Stimulus.spacing,Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_multiple_translate()
{
	static float locspace = 30.0f;
	float mlum=Stimulus.lum;
	float lum1=Stimulus.lum*(1+Stimulus.contrast*r1.readLastValue());
	float lum2=Stimulus.lum*(1+Stimulus.contrast*r2.readLastValue());
	lum1=max(0.0,lum1); lum2=max(0.0,lum2);
	Sfull.x=lum1*100; Sfull.y=lum2*100;  // store values here? perhaps not wisest, but should work.

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	// second, draw one set of panels
	glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
	for (int ii=1; ii<180/locspace; ii++)  // translation undefined straight ahead -- kill it!
	{
		drawWedge(ii*locspace-Stimulus.spacing, Stimulus.spacing);
		drawWedge(ii*-locspace,Stimulus.spacing);
	}


	// third, draw second set of panels
	glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
	for (int ii=1; ii<180/locspace; ii++)
	{
		drawWedge(ii*locspace, Stimulus.spacing);
		drawWedge(ii*-locspace - Stimulus.spacing,Stimulus.spacing);
	}
};

void TrackStim::drawReichardtSlits_deterministic(float timein)
{
	// timein from 0...
	bool PULSE = FALSE; // false = steps
	bool MULT = TRUE;
	static float locspace = 20.0f;
	//float mlum=Stimulus.lum;
	float mlum=0.0f; // background black...
	float lum1=Stimulus.lum*(1+Stimulus.contrast);
	float lum2=Stimulus.lum*(1+Stimulus.contrast*Stimulus.trans.phase/abs(Stimulus.trans.phase)); // phase == -1 inverts it from other contrast...

	// first, set background properly
	glClearColor(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'),0.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	
	float tau = abs(Stimulus.tau);
	bool AFIRST = (Stimulus.tau > 0); // A comes first if tau > 0, otherwise B first...

	//AAAAAAAAAAAAAAAAAAAAAA
	// second, draw one set of panels
	// if tau<0, then A turns on immediately, with B later
	if (  ( ( ((AFIRST) && (timein < tau)) || ((!AFIRST) && (timein > tau) && (timein < 2*tau)) ) && PULSE ) 
		|| (!PULSE && ((AFIRST) || (!AFIRST && (timein > tau) ) )  )  )
	{
		glColor3f(0.0f,getDLPColor(lum1,'G'),getDLPColor(lum1,'B'));
		Sfull.x = lum1;
	}
	else
	{
		glColor3f(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'));
		Sfull.x = mlum;
	};

	drawWedge(90.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
	drawWedge(270.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);

	if (MULT)
	{
		drawWedge(0.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(30.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(60.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(120.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(150.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(210.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(240.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(300.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(330.0f-Stimulus.spacing-Stimulus.spacing2/2, Stimulus.spacing);
	};

	
	// BBBBBBBBBBBBBBBBBBBBBBBBB
	// third, draw second set of panels
	if (  ( ( ((!AFIRST) && (timein < tau)) || ((AFIRST) && (timein > tau) && (timein < 2*tau)) ) && PULSE )
		|| (!PULSE && ((!AFIRST) || (AFIRST && (timein > tau) ) )  )  )
	{
		glColor3f(0.0f,getDLPColor(lum2,'G'),getDLPColor(lum2,'B'));
		Sfull.y = lum2;
	}
	else
	{
		glColor3f(0.0f,getDLPColor(mlum,'G'),getDLPColor(mlum,'B'));
		Sfull.y = mlum;
	};

	drawWedge(90.0f+Stimulus.spacing2/2, Stimulus.spacing);
	drawWedge(270.0f+Stimulus.spacing2/2, Stimulus.spacing);

	if (MULT)
	{
		drawWedge(0+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(30+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(60+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(120+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(150+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(210+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(240+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(300.0f+Stimulus.spacing2/2, Stimulus.spacing);
		drawWedge(330.0f+Stimulus.spacing2/2, Stimulus.spacing);
	};

};

void TrackStim::drawCylinderHalves(float time)
{
	int ii;
	float angle = Stimulus.spacing; // in degrees
	GLfloat Dscale=1;  // this scales teh distance: weird metric to 
	GLfloat startangle;
	//float time = queryCurrTime() - epochtimezero;
	
	//startangle = (DistanceAccumulated * Dscale) 
	//	+ time*Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per);
	startangle = (DistanceAccumulated * Dscale); 
	startangle -= 2*angle*(int)(startangle / 2/angle);  // stupid -- no mod function for double
	if (startangle < 0)
		startangle += 2*angle;

	glRotatef(Stimulus.stimrot.mean,0,0,1);
	
	// first half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle, angle);
		}
		else
		{
			TrackStim::drawWedge(ii*angle*2+startangle, 2*angle - startangle);  // last wedge
			TrackStim::drawWedge(0, startangle - angle);
		};
	}
	// second half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle)
		{
			TrackStim::drawWedge(360 - (ii*angle*2+startangle), -angle);
		}
		else
		{
			TrackStim::drawWedge(360 - (ii*angle*2+startangle), -(2*angle - startangle));  // last wedge
			TrackStim::drawWedge(360 - 0, -(startangle - angle));
		};
	}
}; // halves that rotate


void TrackStim::drawCylinderHalves_fwd(float time)
{
	int ii;
	float angle = Stimulus.spacing; // in degrees
	GLfloat Dscale=180/9;  // this scales 9 mm = 180 degrees, roughly what ball is... 
	GLfloat startangle;
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static float curroldangle = 0;
	static float oldangle = 0;

	if ((time < 0.001) && (!RESET))
	{
		curroldangle += oldangle;  // happens only once per time<.001; don't update for all 3 calls at the time point...
		RESET = TRUE;
	}
	else
		RESET = FALSE;
	
	//startangle = (DistanceAccumulated * Dscale) 
	//	+ time*Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per);
	startangle = (FwdAccumulated * Dscale); 
	startangle += time*Stimulus.stimtrans.mean + curroldangle;
	oldangle = time*Stimulus.stimtrans.mean; // track how far it's gone in this round, each time...

	startangle -= 2*angle*(int)(startangle / 2/angle);  // stupid -- no mod function for double
	if (startangle < 0)
		startangle += 2*angle;

	glRotatef(Stimulus.stimrot.mean,0,0,1);
	
	// first half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle, angle);
		}
		else
		{
			TrackStim::drawWedge(ii*angle*2+startangle, 2*angle - startangle);  // last wedge
			TrackStim::drawWedge(0, startangle - angle);
		};
	}
	// second half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle)
		{
			TrackStim::drawWedge(360 - (ii*angle*2+startangle), -angle);
		}
		else
		{
			TrackStim::drawWedge(360 - (ii*angle*2+startangle), -(2*angle - startangle));  // last wedge
			TrackStim::drawWedge(360 - 0, -(startangle - angle));
		};
	}
}; // halves that rotate

void TrackStim::drawCylinderHalves_fwd_Sawtooth(float time)
{
	int ii;
	float angle = Stimulus.spacing; // in degrees
	GLfloat Dscale=180/9;  // this scales 9 mm = 180 degrees, roughly what ball is... 
	GLfloat startangle;
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static float curroldangle = 0;
	static float oldangle = 0;

	if ((time < 0.001) && (!RESET))
	{
		curroldangle += oldangle;  // happens only once per time<.001; don't update for all 3 calls at the time point...
		RESET = TRUE;
	}
	else
		RESET = FALSE;
	
	//startangle = (DistanceAccumulated * Dscale) 
	//	+ time*Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per);
	startangle = (FwdAccumulated * Dscale); 
	startangle += time*Stimulus.stimtrans.mean + curroldangle;
	oldangle = time*Stimulus.stimtrans.mean; // track how far it's gone in this round, each time...

	startangle -= 2*angle*(int)(startangle / 2/angle);  // stupid -- no mod function for double
	if (startangle < 0)
		startangle += 2*angle;

	glRotatef(Stimulus.stimrot.mean,0,0,1);

	float ST,ND;
	bool SMFWD;
	if (abs(Stimulus.stimtrans.phase - 1) < 0.01)  // which way the saw tooth waves are facing...
	{
		ST = 0; ND = 1;
		SMFWD = TRUE;
	}
	else
	{
		ST = 1; ND = 0;
		SMFWD = FALSE;
	};
	
	// first half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle*2)
		{
			TrackStim::drawSawtooth(ii*angle*2+startangle, angle*2,ST,ND);
		}
		else
		{
			if (SMFWD) // ******* WORKING HERE -- CORRECT...
			{
				TrackStim::drawSawtooth(ii*angle*2+startangle, 2*angle - startangle,0,(2*angle-startangle)/2/angle);  // last wedge
				TrackStim::drawSawtooth(0, startangle,1-(startangle)/angle/2,1);
			}
			else
			{
				TrackStim::drawSawtooth(ii*angle*2+startangle, 2*angle - startangle,1,1-(2*angle-startangle)/2/angle);  // last wedge
				TrackStim::drawSawtooth(0, startangle,(startangle)/angle/2,0);
			}
		};
	}
	// second half of it
	for (ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle < 180 - angle*2)
		{
			TrackStim::drawSawtooth(360 - (ii*angle*2+startangle), -angle*2,ST,ND);
		}
		else
		{
			if (SMFWD)
			{
				TrackStim::drawSawtooth(360 - (ii*angle*2+startangle), -(2*angle - startangle),0,(2*angle-startangle)/2/angle);  // last wedge
				TrackStim::drawSawtooth(360 - 0, -(startangle),1-(startangle)/angle/2,1);
			}
			else
			{
				TrackStim::drawSawtooth(360 - (ii*angle*2+startangle), -(2*angle - startangle),1,1-(2*angle-startangle)/2/angle);  // last wedge
				TrackStim::drawSawtooth(360 - 0, -(startangle),(startangle)/angle/2,0);
			}
		};
	}
}; // halves that rotate

void TrackStim::drawCylinderHemispheres(float time)
{
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0;
	static float  currold_right=0;

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		old_right += currold_right;
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;
	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) +
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_left; // in radians...; 
	float add_right = -(FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) + 
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_right;
	currold_left = Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180; // save to add to next one...
	currold_right = -Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180;

	float startangle_left; // changed 110722 to correct weird drawings at front left
	startangle_left = add_left*180/PI-int(add_left*180/PI/(2*angle))*2*angle; // in degrees for drawWedge below...
	if (startangle_left<0) // negative case
	{	startangle_left += 2*angle;};

	float startangle_right; 
	startangle_right = add_right*180/PI-int(add_right*180/PI/(2*angle))*2*angle;
	if (startangle_right<0) // negative case
	{	startangle_right += 2*angle;};

	// first half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_left < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
		}
		else
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle - startangle_left);  // last wedge
			TrackStim::drawWedge(0, startangle_left - angle); // draw first one
		};
	}
	// second half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_right < 180 - angle)
		{
			TrackStim::drawWedge(180 + (ii*angle*2+startangle_right), -angle);
		}
		else
		{
			TrackStim::drawWedge(180 + (ii*angle*2+startangle_right), -(angle - startangle_right));  // last wedge
			TrackStim::drawWedge(180 + 0, -(startangle_right - angle));
		};
	}
}; 



void TrackStim::drawCylinderHemispheres_revised(float time)
{
	// revised 12/02/02 to remove small glitch near convergence
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0;
	static float currold_right=0;

	float angle_offset = Stimulus.stimrot.phase;

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		old_right += currold_right;
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;
	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) +
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_left; // in radians...; 
	float add_right = -(FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) + 
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_right;
	currold_left = Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180; // save to add to next one...
	currold_right = -Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180;

	float startangle_left; // changed 110722 to correct weird drawings at front left
	startangle_left = add_left*180/PI-int(add_left*180/PI/(2*angle))*2*angle; // in degrees for drawWedge below...
	if (startangle_left<0) // negative case
	{	startangle_left += 2*angle;};

	float startangle_right; 
	startangle_right = add_right*180/PI-int(add_right*180/PI/(2*angle))*2*angle;
	if (startangle_right<0) // negative case
	{	startangle_right += 2*angle;};

	//////////////////////////////
	// DRAW STUFF
	//////////////////////////////
	glPushMatrix();
	glRotatef(angle_offset,0.0f,0.0f,1.0f);
	
	float boundcase = 0;
	// first half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_left < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_left; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(ii*angle*2+startangle_left, 180 - boundcase); } // last wedge
		};
		if (startangle_left > angle)
		{
			TrackStim::drawWedge(0,  startangle_left - angle); // draw first one
		};
	}
	// second half of it
	
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_right < 180 - angle)
		{
			TrackStim::drawWedge(180+ (ii*angle*2+startangle_right), angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_right; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(180+ (boundcase), 180 - boundcase); } // last wedge
		};
		if (startangle_right > angle)
		{
			TrackStim::drawWedge(180, startangle_right - angle);
		};

	}
	glPopMatrix();
	
}; 

void TrackStim::drawCylinderHemispheres_revised_graycenter(float time)
{
	// revised 12/02/02 to remove small glitch near convergence
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0;
	static float currold_right=0;

	float angle_offset = Stimulus.stimrot.phase;

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		old_right += currold_right;
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;
	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) +
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_left; // in radians...; 
	float add_right = -(FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) + 
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_right;
	currold_left = Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180; // save to add to next one...
	currold_right = -Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180;

	float startangle_left; // changed 110722 to correct weird drawings at front left
	startangle_left = add_left*180/PI-int(add_left*180/PI/(2*angle))*2*angle; // in degrees for drawWedge below...
	if (startangle_left<0) // negative case
	{	startangle_left += 2*angle;};

	float startangle_right; 
	startangle_right = add_right*180/PI-int(add_right*180/PI/(2*angle))*2*angle;
	if (startangle_right<0) // negative case
	{	startangle_right += 2*angle;};

	//////////////////////////////
	// DRAW STUFF
	//////////////////////////////
	glPushMatrix();
	glRotatef(angle_offset,0.0f,0.0f,1.0f);
	
	float boundcase = 0;
	// first half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_left < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_left; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(ii*angle*2+startangle_left, 180 - boundcase); } // last wedge
		};
		if (startangle_left > angle)
		{
			TrackStim::drawWedge(0,  startangle_left - angle); // draw first one
		};
	}
	// second half of it
	
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_right < 180 - angle)
		{
			TrackStim::drawWedge(180+ (ii*angle*2+startangle_right), angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_right; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(180+ (boundcase), 180 - boundcase); } // last wedge
		};
		if (startangle_right > angle)
		{
			TrackStim::drawWedge(180, startangle_right - angle);
		};

	};
	
	float centerlum = Stimulus.density;

	// flicker the center screen every tau2/3
	static int centercount = 0;
	int curriter = 0;
	if (Stimulus.tau2 > 1)
	{
		curriter = centercount % int(Stimulus.tau2+0.01);
		if (curriter <= Stimulus.tau2/2)
		{
			centerlum = Stimulus.lum*(1 + Stimulus.contrast);
		}
		else
		{
			centerlum = Stimulus.lum*(1 - Stimulus.contrast);
		};
		centercount += 1;
	};
	glColor3f(0,getDLPColor(centerlum,'G'),getDLPColor(centerlum,'B'));
	if (Stimulus.spacing2 == 0)
	{
		drawWedge_inside(-10,20);
	}
	else
	{
		drawWedge_inside(-Stimulus.spacing2/2,Stimulus.spacing2);
	};
	glColor3f(foregroundtriple[0],foregroundtriple[1],foregroundtriple[2]);
	glPopMatrix();
	
}; 

void TrackStim::drawCylinderHemispheres_revised_graycenter_wait(float time)
{
	// revised 12/02/02 to remove small glitch near convergence
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0;
	static float currold_right=0;

	float angle_offset = Stimulus.stimrot.phase;

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		old_right += currold_right;
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;

	if (time < Stimulus.tau)
	{
		time = 0;
	}
	else
	{
		time -= Stimulus.tau;
	};

	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) +
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_left; // in radians...; 
	float add_right = -(FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) + 
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_right;
	currold_left = Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180; // save to add to next one...
	currold_right = -Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180;

	float startangle_left; // changed 110722 to correct weird drawings at front left
	startangle_left = add_left*180/PI-int(add_left*180/PI/(2*angle))*2*angle; // in degrees for drawWedge below...
	if (startangle_left<0) // negative case
	{	startangle_left += 2*angle;};

	float startangle_right; 
	startangle_right = add_right*180/PI-int(add_right*180/PI/(2*angle))*2*angle;
	if (startangle_right<0) // negative case
	{	startangle_right += 2*angle;};

	//////////////////////////////
	// DRAW STUFF
	//////////////////////////////
	glPushMatrix();
	glRotatef(angle_offset,0.0f,0.0f,1.0f);
	
	float boundcase = 0;
	// first half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_left < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_left; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(ii*angle*2+startangle_left, 180 - boundcase); } // last wedge
		};
		if (startangle_left > angle)
		{
			TrackStim::drawWedge(0,  startangle_left - angle); // draw first one
		};
	}
	// second half of it
	
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_right < 180 - angle)
		{
			TrackStim::drawWedge(180+ (ii*angle*2+startangle_right), angle);
		}
		else
		{
			boundcase = ii*angle*2+startangle_right; // this is angle to start at...
			if (boundcase < 180)
			{	TrackStim::drawWedge(180+ (boundcase), 180 - boundcase); } // last wedge
		};
		if (startangle_right > angle)
		{
			TrackStim::drawWedge(180, startangle_right - angle);
		};

	};
	
	float centerlum = Stimulus.density;

	// flicker the center screen every tau2/6
	static int centercount = 0;
	int curriter = 0;
	if (Stimulus.tau2 > 1)
	{
		curriter = centercount % int(Stimulus.tau2+0.01);
		if (curriter <= Stimulus.tau2/2)
		{
			centerlum = Stimulus.lum*(1 + Stimulus.contrast);
			//Sfull.x = 0; // for testing period
		}
		else
		{
			centerlum = Stimulus.lum*(1 - Stimulus.contrast);
			//Sfull.x = 1;
		};
		centercount += 1; // 120502 -- this loop executes 6 times per frame, 1 each for the windows, x2 for green and blue
	};

	glColor3f(0,getDLPColor(centerlum,'G'),getDLPColor(centerlum,'G'));
	if (Stimulus.spacing2 == 0)
	{
		drawWedge_inside(-10,20);
	}
	else {
		if ((Stimulus.spacing2 > 0) && (Stimulus.spacing2 <90))
		{
			drawWedge_inside(-Stimulus.spacing2/2,Stimulus.spacing2);
		}
		else
		{
			for (int ii=0; ii<36; ii++)
			{
				drawWedge_inside(10*ii,10);
			};
		}
	};
	glColor3f(foregroundtriple[0],foregroundtriple[1],foregroundtriple[2]);
	glPopMatrix();
	
}; 









void TrackStim::drawCylinder_continuous(float time)
{
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float currold_left=0;

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;
	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = (FwdAccumulated * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) +
		-(Scurr.theta + Stimulus.stimrot.mean*time*PI/180) +
		old_left; // in radians...; 
	currold_left = Stimulus.stimtrans.mean*time*PI/180 + -Stimulus.stimrot.mean*time*PI/180; // save to add to next one...

	float startangle_left; 
	startangle_left = mod(add_left*180/PI,2*angle); // in degrees for drawWedge below...


	// first half of it
	for (int ii=0; ii<360/angle/2; ii++)  // this gets the correct number
	{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
	};
}; 


void TrackStim::drawCylinderHemispheresGain(float time)
{
	float angle = Stimulus.spacing; // in degrees
	static bool RESET = FALSE; // get things to reset properly, keep track of things...
	static int loopnum = 0;
	static float old_left =0;
	static float old_right=0;
	static float currold_left=0; // these keep track of the phases of the total OL stimuli that have come before
	static float  currold_right=0;
	static float fwd0 = 0; // store initial fwd and theta for each bit...
	static float theta0 = 0; 

	if ((time < 0.001) && (!RESET))
	{
		old_left += currold_left;  // happens only once per time<.001
		old_right += currold_right;
		fwd0 = FwdAccumulated;
		theta0 = Scurr.theta;
		RESET = TRUE;
		loopnum = 0;
	}
	else
		if (loopnum > 6) // only do this after it's done all the loops through this... for that first color frame = 3 renderings.
			RESET = FALSE;
	loopnum ++;
	
	float Dfactor = 180/9; // 180 degrees per 9 mm moved...
	float add_left = ( (FwdAccumulated - fwd0) * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) * Stimulus.left.transgain +
		-( (Scurr.theta - theta0) + Stimulus.stimrot.mean*time*PI/180) * Stimulus.left.rotgain +
		old_left; // in radians...; 
	float add_right = -( (FwdAccumulated - fwd0) * Dfactor *PI/180 + Stimulus.stimtrans.mean*time*PI/180) * Stimulus.right.transgain + 
		-( (Scurr.theta - theta0) + Stimulus.stimrot.mean*time*PI/180) * Stimulus.right.rotgain +
		old_right;
	currold_left = add_left - old_left; // save to add to next one...
	currold_right = add_right - old_right;


	float startangle_left = mod(add_left*180/PI,2*angle); // in degrees for drawWedge below...
	float startangle_right = mod(add_right*180/PI,2*angle);

	// first half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_left < 180 - angle)
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, angle);
		}
		else
		{
			TrackStim::drawWedge(ii*angle*2+startangle_left, 2*angle - startangle_left);  // last wedge
			TrackStim::drawWedge(0, startangle_left - angle);
		};
	}
	// second half of it
	for (int ii=0; ii<180/angle/2; ii++)  // this gets the correct number
	{
		if (ii*angle*2 + startangle_right < 180 - angle)
		{
			TrackStim::drawWedge(180 + (ii*angle*2+startangle_right), -angle);
		}
		else
		{
			TrackStim::drawWedge(180 + (ii*angle*2+startangle_right), -(2*angle - startangle_right));  // last wedge
			TrackStim::drawWedge(180 + 0, -(startangle_right - angle));
		};
	}
}; 


void TrackStim::drawInfCorridor(float currtime)
{
	float Dscale = 1; 
	float extent=30*Stimulus.arenasize;
	float depth=30*Stimulus.arenaheight;
	float width=10;
	float spacing = Stimulus.spacing;
	int ii;
	//float time=TrackStim::queryCurrTime();;
	float time = currtime;
	float startpos;
	//float added = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // add the stimulus OL
	
	
	startpos= (DistanceAccumulated * Dscale) 
		+ time*Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per);
	startpos -= spacing*(int)(startpos / spacing); // stupid -- no mod function for double

	glPushMatrix();

	//glRotatef(added,0.0f,0.0f,1.0f);
	glTranslatef(-extent-startpos,0,0);
	for (ii=-(int)(extent/spacing); ii<=(int)(extent/spacing); ii++)
	{
		glBegin(GL_QUADS);
		glVertex3f(0,width,0.0f);
		glVertex3f(0,width,depth);
		glVertex3f((0.5f)*spacing,width,depth);
		glVertex3f((0.5f)*spacing,width,0.0f);
		glEnd();
		glBegin(GL_QUADS);
		glVertex3f(0,-width,0.0f);
		glVertex3f(0,-width,depth);
		glVertex3f((0.5f)*spacing,-width,depth);
		glVertex3f((0.5f)*spacing,-width,0.0f);
		glEnd();
		glTranslatef(spacing,0,0);
	};
	glPopMatrix();
}; // realistic translation

void TrackStim::drawLinearRDots(float dt)
{
	//static float lasttime=TrackStim::queryCurrTime();
	//float currtime = TrackStim::queryCurrTime();
	Rdots.tau = Stimulus.tau;
	Rdots.Zheight = Stimulus.arenaheight; // height factor from cylinder drawing
	Rdots.windowWidth = Stimulus.arenasize;
	Rdots.windowHeight = Stimulus.arenasize;
	Rdots.blobSize = Stimulus.spacing*PI/180*Stimulus.arenaheight;  // in degrees now
	int temp,temp2;
	temp=(int)(Rdots.windowWidth/Rdots.blobSize);
	temp2=(int)(Rdots.windowWidth/Rdots.blobSize);
	Rdots.blobCount = (int)((float)temp*(float)temp2 * -1*log(1-Stimulus.density));  // roughly the number on it
	
	Rdots.speed = (Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per)) * PI/180 * Stimulus.arenaheight; // give it speed to update with -- degrees/s @ nearest point
	float angle = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // gets angle right
	glRotatef(angle,0,0,1.0f);
	Rdots.RenderScene(dt);
}; // random dots

void TrackStim::drawRotRDots(float dt)
{
	//static float lasttime=TrackStim::queryCurrTime();
	//float currtime = TrackStim::queryCurrTime();
	Rdots.tau = Stimulus.tau;
	Rdots.Zheight = Stimulus.arenaheight*10; // height factor from cylinder drawing
	int temp,temp2;
	temp=(int)(Stimulus.arenaheight*10/(Stimulus.spacing*PI/180*30)); // 30 is radius of the circle, st. 30*angle = linear dimension...
	temp2=(int)(360/Stimulus.spacing);
	Rdots.blobCount = (int)((float)temp*(float)temp2 * -1*log(1-Stimulus.density));  // roughly the number on it
	Rdots.blobSize = Stimulus.spacing;  // in degrees
	Rdots.speed = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // give it speed to update with
	Rdots.speed *= PI/180;
	Rdots.RenderScene_rotation(dt);
	//lasttime=currtime;
}; // rotational random dots

void TrackStim::drawRotRDotsHalves(float dt)
{
	//static float lasttime=TrackStim::queryCurrTime();
	//float currtime = TrackStim::queryCurrTime();
	Rdots.speed = Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per); // give it speed to update with
	Rdots.RenderScene_rotation_halves(dt);
	//lasttime=currtime;
};

void TrackStim::drawRotRDotsGrads(float dt)
{
	//static float lasttime=TrackStim::queryCurrTime();
	//float currtime = TrackStim::queryCurrTime();
	Rdots.speed = Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per); // give it speed to update with
	Rdots.RenderScene_rotation_grad(dt);
	//lasttime=currtime;
};

void TrackStim::drawBarCeiling(float time) // moving bar stimulus
{
	float spacing=Stimulus.spacing;
	int ii;
	float extent=3*Stimulus.arenasize;
	float width=3*Stimulus.arenasize;
	//float time;
	float startpos;
	
	//time=TrackStim::queryCurrTime();
	startpos = time*Stimulus.stimtrans.mean + Stimulus.stimtrans.amp*TrackStim::sineWave(Stimulus.stimtrans.per);  // OL part
	startpos -= spacing*(int)(startpos / spacing);  // stupid -- no mod function for double
	float angle = Stimulus.stimrot.mean + Stimulus.stimrot.amp*TrackStim::sineWave(Stimulus.stimrot.per); // gets angle right

	glPushMatrix();
	glRotatef(angle,0,0,1.0f);
	
	glTranslatef(-extent/2-startpos,0,Stimulus.arenaheight);
	for (ii=-(int)(extent/spacing); ii<(int)(extent/spacing); ii++)
	{
		glBegin(GL_QUADS);
		glVertex3f(0,-width,0.0f);
		glVertex3f(spacing/2,-width,0.0f);
		glVertex3f(spacing/2,width,0.0f);
		glVertex3f(0,width,0.0f);
		glEnd();
		glTranslatef(spacing,0,0);

	};
	glPopMatrix();
};  


// 131002 Altered from scintilate function, case number 127 in the stimulus text file
void TrackStim::drawStarField(float tc)
{
	int ii;
	static bool FIRST = true;
#define X_MAX 5
#define Y_MAX 5
	static std::vector<StarClass::Star> stars;
	float imagedata[X_MAX * Y_MAX * 3] = { 0.0f };
	static float imagetime[X_MAX * Y_MAX] = { 0.0f };
	static int x_next = 0;
	static int y_next = 0;
	static int x_curr = 0;
	static int y_curr = 0;
	static float last_tc = 0;
	float dtime = 0.5; // seconds between transition
	float startwidth = Stimulus.spacing; // starts at 2 degrees, 1/2 thickness of the bars?

	
	if (FIRST)
	{
		stars.resize(1);
		last_tc = tc;

		generateNextStar(stars[0]);

		x_curr = stars[0].x;
		y_curr = stars[0].y;

		for (int xi = 0; xi < X_MAX; xi++)
		{
			for (int yi = 0; yi < Y_MAX; yi++)
			{
				if (xi == x_curr && yi == y_curr)
				{
					imagedata[3 * (xi + yi * X_MAX) + 0] = foregroundtriple[0];
					imagedata[3 * (xi + yi * X_MAX) + 1] = foregroundtriple[1];
					imagedata[3 * (xi + yi * X_MAX) + 2] = foregroundtriple[2];
				}
				else
				{
					imagedata[3 * (xi + yi * X_MAX) + 0] = 0.0; // backgroundtriple[0];
					imagedata[3 * (xi + yi * X_MAX) + 1] = 0.5; // backgroundtriple[1];
					imagedata[3 * (xi + yi * X_MAX) + 2] = 0.5; // backgroundtriple[2];
				}
			}
		}

		imagetime[x_curr + y_curr * 2] = dtime;

		FIRST = false;
	}

	if (tc - last_tc > 0.001)
	{
		imagetime[x_curr + y_curr * 2] -= tc - last_tc;

		last_tc = tc;

		if (imagetime[x_curr + y_curr * 2] <= 0)
		{

			generateNextStar(stars[0]);
			x_curr = 1; //stars[0].x;
			y_curr = 1; // stars[0].y;

			imagetime[x_curr + y_curr * 2] = dtime;
		}

		for (int xi = 0; xi < X_MAX; xi++)
		{
			for (int yi = 0; yi < Y_MAX; yi++)
			{
				if (xi == x_curr && yi == y_curr)
				{
					imagedata[3 * (xi + yi * X_MAX) + 0] = foregroundtriple[0];
					imagedata[3 * (xi + yi * X_MAX) + 1] = foregroundtriple[1];
					imagedata[3 * (xi + yi * X_MAX) + 2] = foregroundtriple[2];
				}
				else
				{
					imagedata[3 * (xi + yi * X_MAX) + 0] = 0.0;// backgroundtriple[0];
					imagedata[3 * (xi + yi * X_MAX) + 1] = 0.5;// backgroundtriple[1];
					imagedata[3 * (xi + yi * X_MAX) + 2] = 0.5;// backgroundtriple[2];
				}
			}
		}
	}

	glClearColor(meangroundtriple[0], meangroundtriple[1], meangroundtriple[2], 0);
	glClear(GL_COLOR_BUFFER_BIT);

	
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1); // rebind the new texture...
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, X_MAX, Y_MAX, 0, GL_RGB, GL_FLOAT, imagedata);

	float testwidth = 1.5;
	float testheight = 1.5;
	float testdepth = -2;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 1);
	glBegin(GL_QUADS);
		// tl
		// bl
		// br
		// tr
	glTexCoord2f(0.0f, 1.0f); //glTexCoord2f(0.0f, (float)Y_MAX - 1);
	glVertex3f(0.0f, testheight, testdepth);
	glTexCoord2f(0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, testdepth);
	glTexCoord2f(1.0f, 0.0f); //glTexCoord2f((float)X_MAX - 1, 0.0f);
	glVertex3f(testwidth, 0.0f, testdepth);
	glTexCoord2f(1.0f, 1.0f); //glTexCoord2f((float)X_MAX - 1, (float)Y_MAX - 1);
	glVertex3f(testwidth, testheight, testdepth);

		// glVertex3f(0.0f, testheight, 1.f);
		// glVertex3f(testwidth, testheight, 1.f);
		// glVertex3f(testwidth, 0.0f, 1.f);
		// glVertex3f(0.0f, 0.0f, 1.f);
		
		/*
		glVertex3f(28.0f, 10.3f, 0.0f);
		glVertex3f(25.0f, 15.0f, 0.0f);
		glVertex3f(28.0f, 10.3f, 300.0f);
		glVertex3f(25.0f, 15.0f, 300.0f);
		*/
	glEnd();
	glDisable(GL_TEXTURE_2D);


	/*
	for (ii=0; ii<360.01/angle; ii++) // draw all the bars
	{
		TrackStim::drawWedge(ii*angle-startwidth/2+added,startwidth);
	};*/
};

void TrackStim::updateStars(float tc, StarClass::Star star[])
{
	//for (int ii = 0; ii<)
}

void TrackStim::generateNextStar(StarClass::Star &star)
{
	bool conflict_exists = false;
	float x_next;
	float y_next;
	float x_curr = star.x;
	float y_curr = star.y;

	conflict_exists = true;
	while (conflict_exists)
	{
		x_next = rand() % (X_MAX - 1);
		y_next = rand() % (Y_MAX);
		if (y_next == y_curr && abs(x_next - x_curr) <= 1)
			x_next = rand() % (X_MAX - 1);
		else
			conflict_exists = false;
	}
	star.x = x_next;
	star.y = y_next;
}

/* Case 128
*/
void TrackStim::drawCylinderBarsES(float tc)
{
	int ii;
	//float startpos = Stimulus.spacing; // degrees, where first bar is drawn
	float startpos = 180;
	float startpos_i;
	float width = abs(Stimulus.width); // degrees
	float direction = width / Stimulus.width; // Sign
	float r = 30; // arbitrary radius, doesn't matter for now since everything else is in degrees
	float added; // 090302, added this to get rid of stimrot.mean rotation, not useful for me here...

	if (tc <= Stimulus.tau2) // pause time before things move
		added = 0;
	else
		added = (tc - Stimulus.tau2)*80;  // put in translation mean here to get easy spinning. not ideal... bring amp back to 0 though!

	float height = 5*r; // arbitrary number just to cover the screen
	float x1;
	float z1;
	float x2;
	float z2;

	startpos_i = startpos; //- 80*(tc - 0.55) * direction;
	glPushMatrix();
	for (int i=0; i < 360.1/(abs(width)*2); i++)
	{
		
		float x1 = r*cos(startpos_i*3.14159 / 180);
		float z1 = -r*sin(startpos_i*3.14159 / 180);
		float x2 = r*cos((startpos_i - width)*3.14159 / 180);
		float z2 = -r*sin((startpos_i - width)*3.14159 / 180);

		glEnable(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 1);
		glBegin(GL_QUADS);
		// Draw the vertices in CCW order
		glVertex3f(x1, height / 2, z1); // top left
		glVertex3f(x1, -height / 2, z1); // bottom left
		glVertex3f(x2, -height / 2, z2); // bottom right
		glVertex3f(x2, height / 2, z2); // top right
		glEnd();

		startpos_i -= 2 * width;
		

		//TrackStim::drawWedge(i*(2 * width) - width / 2 + added, width);
	};
	glPopMatrix();
};

// Case 129
void TrackStim::drawSimpleBarCW(float tc)
{
	// Start and end positions of bars in degrees
	float startpos[6];
	float endpos[6];
	int numObjects = 6;
	for (int ii = 0; ii < numObjects; ii++)
	{
		switch (ii)
		{
			case 0:
				startpos[ii] = Stimulus.startpos1;
				break;
			case 1:
				startpos[ii] = Stimulus.startpos2;
				break;
			case 2:
				startpos[ii] = Stimulus.startpos3;
				break;
			case 3:
				startpos[ii] = Stimulus.startpos4;
				break;
			case 4:
				startpos[ii] = Stimulus.startpos5;
				break;
			case 5:
				startpos[ii] = Stimulus.startpos6;
				break;
		};
		endpos[ii] = startpos[ii] - Stimulus.spacing;
	};

	float drawpos;
	float width = Stimulus.width; // degrees

	float r = 2.5; // arbitrary radius, doesn't matter for now since everything else is in degrees
	float height = 5 * r; // arbitrary height just to cover the screen
	// Coordinates for 4 vertices of bar
	float x1;
	float z1;
	float x2;
	float z2;

	// Absolute timings
	float t1 = 0; // Bars appear at initial position (default zero for now, could start with blank screen later)
	float t2 = Stimulus.tau; // Bars disappear
	float t3 = Stimulus.tau + Stimulus.tau2; // Bars appear at new position
	float t4 = 2*Stimulus.tau + Stimulus.tau2; // Bars disappear again

	glClearColor(meangroundtriple[0], meangroundtriple[1], meangroundtriple[2], 0);

	if ( (tc >= t1 && tc < t2) || (tc >= t3 && tc < t4))
	{
		glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);

		//glPushMatrix();
		for (int ii = 0; ii < numObjects; ii++) // Loop through screens
		{
			// Select the correct position to draw the bars
			drawpos = (tc < t2) ? startpos[ii] : endpos[ii];
			
			float x1 = r*cos(drawpos * PI / 180);
			float z1 = -r*sin(drawpos * PI / 180);
			float x2 = r*cos((drawpos - width)*PI / 180);
			float z2 = -r*sin((drawpos - width)*PI / 180);
					

			glBegin(GL_QUADS);
			// Draw the vertices in CCW order
			glVertex3f(x1, height / 2, z1); // top left
			glVertex3f(x1, -height / 2, z1); // bottom left
			glVertex3f(x2, -height / 2, z2); // bottom right
			glVertex3f(x2, height / 2, z2); // top right
			glEnd();
		};
		//glPopMatrix();
		
	};
};

// Case 130
void TrackStim::drawSimpleBarCCW(float tc)
{
	// Start and end positions of bars in degrees
	int numObjects = 6;
	float startpos[6];
	float endpos[6];
	for (int ii = 0; ii < numObjects; ii++)
	{
		switch (ii)
		{
		case 0:
			startpos[ii] = Stimulus.startpos1;
			break;
		case 1:
			startpos[ii] = Stimulus.startpos2;
			break;
		case 2:
			startpos[ii] = Stimulus.startpos3;
			break;
		case 3:
			startpos[ii] = Stimulus.startpos4;
			break;
		case 4:
			startpos[ii] = Stimulus.startpos5;
			break;
		case 5:
			startpos[ii] = Stimulus.startpos6;
			break;
		};
		endpos[ii] = startpos[ii] + Stimulus.spacing;
	};

	float drawpos;
	float width = Stimulus.width; // degrees

	float r = 30; // arbitrary radius, doesn't matter for now since everything else is in degrees
	float height = 5 * r; // arbitrary height just to cover the screen
	// Coordinates for 4 vertices of bar
	float x1;
	float z1;
	float x2;
	float z2;

	// Absolute timings
	float t1 = 0; // Bars appear at initial position (default zero for now, could start with blank screen later)
	float t2 = Stimulus.tau; // Bars disappear
	float t3 = Stimulus.tau + Stimulus.tau2; // Bars appear at new position
	float t4 = 2 * Stimulus.tau + Stimulus.tau2; // Bars disappear again

	glClearColor(meangroundtriple[0], meangroundtriple[1], meangroundtriple[2], 0);
	glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);

	if ((tc >= t1 && tc < t2) || (tc >= t3 && tc < t4))
	{
		//glPushMatrix();
		for (int ii = 0; ii < numObjects; ii++) // Loop through screens
		{
			// Select the correct position to draw the bars
			drawpos = (tc < t2) ? startpos[ii] : endpos[ii];

			float x1 = r*cos(drawpos * PI / 180);
			float z1 = -r*sin(drawpos * PI / 180);
			float x2 = r*cos((drawpos + width)*PI / 180);
			float z2 = -r*sin((drawpos + width)*PI / 180);


			glBegin(GL_QUADS);
			// Draw the vertices in CCW order
			glVertex3f(x1, height / 2, z1); // top left
			glVertex3f(x1, -height / 2, z1); // bottom left
			glVertex3f(x2, -height / 2, z2); // bottom right
			glVertex3f(x2, height / 2, z2); // top right
			glEnd();
		};
		//glPopMatrix();

	};
};


// Case 131
void TrackStim::drawSimpleBarCW_rphi(float tc)
{
	// Start and end positions of bars in degrees
	float startpos[6];
	float endpos[6];
	int numObjects = 6;
	for (int ii = 0; ii < numObjects; ii++)
	{
		switch (ii)
		{
		case 0:
			startpos[ii] = Stimulus.startpos1;
			break;
		case 1:
			startpos[ii] = Stimulus.startpos2;
			break;
		case 2:
			startpos[ii] = Stimulus.startpos3;
			break;
		case 3:
			startpos[ii] = Stimulus.startpos4;
			break;
		case 4:
			startpos[ii] = Stimulus.startpos5;
			break;
		case 5:
			startpos[ii] = Stimulus.startpos6;
			break;
		};
		endpos[ii] = startpos[ii] - Stimulus.spacing;
	};

	float drawpos;
	float width = Stimulus.width; // degrees

	float r = 2.5; // arbitrary radius, doesn't matter for now since everything else is in degrees
	float height = 5 * r; // arbitrary height just to cover the screen
	// Coordinates for 4 vertices of bar
	float x1;
	float z1;
	float x2;
	float z2;

	// Absolute timings
	float t1 = 0; // Bars appear at initial position (default zero for now, could start with blank screen later)
	float t2 = Stimulus.tau; // Bars disappear
	float t3 = Stimulus.tau + Stimulus.tau2; // Bars appear at new position
	float t4 = 2 * Stimulus.tau + Stimulus.tau2; // Bars disappear again

	glClearColor(meangroundtriple[0], meangroundtriple[1], meangroundtriple[2], 0);

	if ((tc >= t1 && tc < t2) || (tc >= t3 && tc < t4))
	{

		//glPushMatrix();
		for (int ii = 0; ii < numObjects; ii++) // Loop through screens
		{
			// Select the correct position and colors to draw the bars
			if (tc < t2)
			{
				drawpos = startpos[ii];
				glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);
			}
			else
			{
				drawpos = endpos[ii];
				glColor3f(backgroundtriple[0], backgroundtriple[1], backgroundtriple[2]);
			}

			float x1 = r*cos(drawpos * PI / 180);
			float z1 = -r*sin(drawpos * PI / 180);
			float x2 = r*cos((drawpos - width)*PI / 180);
			float z2 = -r*sin((drawpos - width)*PI / 180);


			glBegin(GL_QUADS);
			// Draw the vertices in CCW order
			glVertex3f(x1, height / 2, z1); // top left
			glVertex3f(x1, -height / 2, z1); // bottom left
			glVertex3f(x2, -height / 2, z2); // bottom right
			glVertex3f(x2, height / 2, z2); // top right
			glEnd();
		};
		//glPopMatrix();

	};
};

// Case 132
void TrackStim::drawSimpleBarCCW_rphi(float tc)
{
	// Start and end positions of bars in degrees
	int numObjects = 6;
	float startpos[6];
	float endpos[6];
	for (int ii = 0; ii < numObjects; ii++)
	{
		switch (ii)
		{
		case 0:
			startpos[ii] = Stimulus.startpos1;
			break;
		case 1:
			startpos[ii] = Stimulus.startpos2;
			break;
		case 2:
			startpos[ii] = Stimulus.startpos3;
			break;
		case 3:
			startpos[ii] = Stimulus.startpos4;
			break;
		case 4:
			startpos[ii] = Stimulus.startpos5;
			break;
		case 5:
			startpos[ii] = Stimulus.startpos6;
			break;
		};
		endpos[ii] = startpos[ii] + Stimulus.spacing;
	};

	float drawpos;
	float width = Stimulus.width; // degrees

	float r = 30; // arbitrary radius, doesn't matter for now since everything else is in degrees
	float height = 5 * r; // arbitrary height just to cover the screen
	// Coordinates for 4 vertices of bar
	float x1;
	float z1;
	float x2;
	float z2;

	// Absolute timings
	float t1 = 0; // Bars appear at initial position (default zero for now, could start with blank screen later)
	float t2 = Stimulus.tau; // Bars disappear
	float t3 = Stimulus.tau + Stimulus.tau2; // Bars appear at new position
	float t4 = 2 * Stimulus.tau + Stimulus.tau2; // Bars disappear again

	glClearColor(meangroundtriple[0], meangroundtriple[1], meangroundtriple[2], 0);
	glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);

	if ((tc >= t1 && tc < t2) || (tc >= t3 && tc < t4))
	{
		//glPushMatrix();
		for (int ii = 0; ii < numObjects; ii++) // Loop through screens
		{
			// Select the correct position and colors to draw the bars
			if (tc < t2)
			{
				drawpos = startpos[ii];
				glColor3f(foregroundtriple[0], foregroundtriple[1], foregroundtriple[2]);
			}
			else
			{
				drawpos = endpos[ii];
				glColor3f(backgroundtriple[0], backgroundtriple[1], backgroundtriple[2]);
			}

			float x1 = r*cos(drawpos * PI / 180);
			float z1 = -r*sin(drawpos * PI / 180);
			float x2 = r*cos((drawpos + width)*PI / 180);
			float z2 = -r*sin((drawpos + width)*PI / 180);


			glBegin(GL_QUADS);
			// Draw the vertices in CCW order
			glVertex3f(x1, height / 2, z1); // top left
			glVertex3f(x1, -height / 2, z1); // bottom left
			glVertex3f(x2, -height / 2, z2); // bottom right
			glVertex3f(x2, height / 2, z2); // top right
			glEnd();
		};
		//glPopMatrix();

	};
};

/// read in stimulus file, output from excel file
TrackStim::STIMPARAMS * TrackStim::readParams(char * szFile)
{
	FILE *filein;
	int numlines;
	char linename[255];
	char oneline[255];
	float pcurr[40];
	STIMPARAMS Stims[40];
	for (int jj=0; jj<40; jj++) Stims[jj] = TrackStim::initializeStimulus();

	string currlinename=" "; // initialize string objects once
	string complinename=" ";
	string currline=" ";

	//filein = fopen("data/feedbackparams.txt", "rt");				// File To Load Data From
	filein = fopen(szFile, "rt");				// File To Load Data From
	for (int ii=0; ii<260; ii++) TrackStim::paramfilename[ii]=szFile[ii];
	
	// read in initial values about file
	TrackStim::readstring(filein,oneline);
	sscanf(oneline, "PARAMS %d", &numlines);
	TrackStim::readstring(filein,oneline);
	sscanf(oneline, "EPOCHS %d", &numepochs);
	
	for (int loop = 0; loop < numlines; loop++)   // loop on number of parameters
	{
		// go through remaining lines...
		fscanf(filein, "%s", &linename);
		currlinename=linename;

		for (int jj = 0; jj < numepochs; jj++)    // loop on number of parameter sets
		{
			fscanf(filein, "%f", &pcurr[jj]);   // this appears to work
		};
		
		for (int jj = 0; jj < numepochs; jj++)
		{
			complinename="Stimulus.stimtype";
			if ((currlinename == complinename))  // does this work, or do i have to declare a new one here?
				Stims[jj].stimtype = pcurr[jj];
			complinename="Stimulus.lum";
			if ((currlinename == complinename))
				Stims[jj].lum = pcurr[jj];
			complinename="Stimulus.contrast";
			if ((currlinename == complinename))
				Stims[jj].contrast = pcurr[jj];
			complinename="Stimulus.duration";
			if ((currlinename == complinename))
				Stims[jj].duration = pcurr[jj];

			complinename="Stimulus.transgain";
			if ((currlinename == complinename))
				Stims[jj].transgain = pcurr[jj];
			complinename="Stimulus.rotgain";
			if ((currlinename == complinename))
				Stims[jj].rotgain = pcurr[jj];

			complinename="Stimulus.trans.mean";
			if ((currlinename == complinename))
				Stims[jj].trans.mean = pcurr[jj];
			complinename="Stimulus.trans.amp";
			if ((currlinename == complinename))
				Stims[jj].trans.amp = pcurr[jj];
			complinename="Stimulus.trans.per";
			if ((currlinename == complinename))
				Stims[jj].trans.per = pcurr[jj];
			complinename="Stimulus.trans.phase";
			if ((currlinename == complinename))
				Stims[jj].trans.phase = pcurr[jj];

			complinename="Stimulus.rot.mean";
			if ((currlinename == complinename))
				Stims[jj].rot.mean = pcurr[jj];
			complinename="Stimulus.rot.amp";
			if ((currlinename == complinename))
				Stims[jj].rot.amp = pcurr[jj];
			complinename="Stimulus.rot.per";
			if ((currlinename == complinename))
				Stims[jj].rot.per = pcurr[jj];
			complinename="Stimulus.rot.phase";
			if ((currlinename == complinename))
				Stims[jj].rot.phase = pcurr[jj];

			complinename="Stimulus.stimrot.mean";
			if ((currlinename == complinename))
				Stims[jj].stimrot.mean = pcurr[jj];
			complinename="Stimulus.stimrot.amp";
			if ((currlinename == complinename))
				Stims[jj].stimrot.amp = pcurr[jj];
			complinename="Stimulus.stimrot.per";
			if ((currlinename == complinename))
				Stims[jj].stimrot.per = pcurr[jj];
			complinename="Stimulus.stimrot.phase";
			if ((currlinename == complinename))
				Stims[jj].stimrot.phase = pcurr[jj];

			complinename="Stimulus.stimtrans.mean";
			if ((currlinename == complinename))
				Stims[jj].stimtrans.mean = pcurr[jj];
			complinename="Stimulus.stimtrans.amp";
			if ((currlinename == complinename))
				Stims[jj].stimtrans.amp = pcurr[jj];
			complinename="Stimulus.stimtrans.per";
			if ((currlinename == complinename))
				Stims[jj].stimtrans.per = pcurr[jj];
			complinename="Stimulus.stimtrans.phase";
			if ((currlinename == complinename))
				Stims[jj].stimtrans.phase = pcurr[jj];

			complinename="Stimulus.spacing";
			if ((currlinename == complinename))
				Stims[jj].spacing = pcurr[jj];
			complinename="Stimulus.spacing2";
			if ((currlinename == complinename))
				Stims[jj].spacing2 = pcurr[jj];
			complinename="Stimulus.density";
			if ((currlinename == complinename))
				Stims[jj].density = pcurr[jj];
			complinename="Stimulus.tau";
			if ((currlinename == complinename))
				Stims[jj].tau = pcurr[jj];
			complinename="Stimulus.tau2";
			if ((currlinename == complinename))
				Stims[jj].tau2 = pcurr[jj];

			complinename="Stimulus.arenasize";
			if ((currlinename == complinename))
				Stims[jj].arenasize = pcurr[jj];
			complinename="Stimulus.arenaheight";
			if ((currlinename == complinename))
				Stims[jj].arenaheight = pcurr[jj];

			complinename="Stimulus.right.rotgain";
			if ((currlinename == complinename))
				Stims[jj].right.rotgain = pcurr[jj];
			complinename="Stimulus.right.transgain";
			if ((currlinename == complinename))
				Stims[jj].right.transgain = pcurr[jj];
			complinename="Stimulus.left.rotgain";
			if ((currlinename == complinename))
				Stims[jj].left.rotgain = pcurr[jj];
			complinename="Stimulus.left.transgain";
			if ((currlinename == complinename))
				Stims[jj].left.transgain = pcurr[jj];

			complinename="Stimulus.lumsigma";
			if ((currlinename == complinename))
				Stims[jj].lumsigma = pcurr[jj];

			complinename="Stimulus.changemethod";
			if ((currlinename == complinename))
				Stims[jj].changemethod = pcurr[jj];

			complinename="Stimulus.aperture.location";
			if ((currlinename == complinename))
				Stims[jj].aperture.location = pcurr[jj];
			complinename="Stimulus.aperture.width";
			if ((currlinename == complinename))
				Stims[jj].aperture.width = pcurr[jj];

			complinename = "Stimulus.startpos1";
			if ((currlinename == complinename))
				Stims[jj].startpos1 = pcurr[jj];
			complinename = "Stimulus.startpos2";
			if ((currlinename == complinename))
				Stims[jj].startpos2 = pcurr[jj];
			complinename = "Stimulus.startpos3";
			if ((currlinename == complinename))
				Stims[jj].startpos3 = pcurr[jj];
			complinename = "Stimulus.startpos4";
			if ((currlinename == complinename))
				Stims[jj].startpos4 = pcurr[jj];
			complinename = "Stimulus.startpos5";
			if ((currlinename == complinename))
				Stims[jj].startpos5 = pcurr[jj];
			complinename = "Stimulus.startpos6";
			if ((currlinename == complinename))
				Stims[jj].startpos6 = pcurr[jj];
			complinename = "Stimulus.width";
			if ((currlinename == complinename))
				Stims[jj].width = pcurr[jj];
		};

	}

	fclose(filein);

	return Stims;
};

