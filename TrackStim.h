#pragma once
#include <ctime>
#include <cmath>
#include <windows.h>		// Header File For Windows
#include <math.h>			// Math Library Header File
#include <stdio.h>			// Header File For Standard Input/Output
#include <string>
#include <vector>
#include <algorithm>
//#include <iostream>
using namespace std;
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
//#include <gl\glaux.h>		// Header File For The Glaux Library
#include "Animation.h"      // header file for animation class to make random dot patterns.
#include "RandomC.h"
#include "EnforcedCorr1D.h"


class TrackStim
	/*
	Class for creating all the stimuli on the 200 Hz screen (can change to DLP)...
	*/
{

public:

	float ScreenUpdateRate; // throw up at 200 Hz

	struct POSITIONINFO {  // for fly or stimulus
		float x,y,theta; //  current position, THETA IN RADIANS
		float vx,vy,vtheta; //  velocities
		float xold,yold,thetaold; // previous measured position
		float dx,dy,dtheta;
		float tcurr,told; // times of current and prev measurement
		float dt; // tcurr-told
		float dfwd,dlat;
	} ;

	POSITIONINFO Fmeas;  // public "fly measured"

	bool PARAMSFROMFILE; // read parameters from file

	bool VIEWCUBECORNER;
	bool USE_FRUSTUM;

	struct OLPARAMS {
		float mean, amp, per, phase;
	};

	struct GPARAMS {
		float rotgain, transgain;
	};

	struct APPARAMS {
		float location, width;
	};

	struct STIMPARAMS { // everything attached to the stimulus; array of this to get stim sequence by index: "stim epoch"
		int stimtype;
		float lum;
		float lumsigma;
		float contrast;
		float duration;
		float transgain,rotgain;
		// below: changes fly's frame with OL
		OLPARAMS trans;  // OL on fly translation & rotation -- OL on stimulus position and angle
		OLPARAMS rot;
		// below: in fly's frame
		OLPARAMS stimtrans;  // OL on stimulus translation, where applicable: halves, bars, inf corridor
		OLPARAMS stimrot;    // OL on stim rotation, where applicable: bars
		float spacing, spacing2;
		float density;
		float tau, tau2;
		float arenasize, arenaheight;  // arena parameters
		GPARAMS right, left;
		int changemethod;
		APPARAMS aperture;
		float startpos1, startpos2, startpos3, startpos4, startpos5, startpos6;
		float width;
	} Stimulus ;

	int numepochs; // number of epochs in stimulus
	int epochchoose; // control of which after reading them in...
	float epochtimezero;

	struct VIEWPOSITIONS
	{
		int x[4],y[4],h[4],w[4];  // (x,y,h,w) coords for all 4 view windows...
	} ViewPorts;

	struct Cposition
	{
		float x[2],y[2],z[2];
	} ScreenPosition;

	class StarClass
	{
	public:
		// Constructor
		StarClass()
		{
			starArray.resize(2);
		};

		static int starCount; // Number of stars

		struct Star
		{
			float timeToDeath = 0; // seconds, the star is destroyed when it reaches 0
			int generation = 0; // Determines if it's the first or second of star pair

			// Position of star
			float x = 0;
			float y = 0;
			float z = 0;

			// Size of star
			float xdim = 0;
			float ydim = 0;
		};
		
		std::vector<Star> starArray;
	};

private:
	
	// basic
	float PI;

	POSITIONINFO Fcurr; // private "fly current"
	POSITIONINFO Scurr; // private "stim current"
	POSITIONINFO Sfull; // contains open loop info, too
	POSITIONINFO Saccumulate; // for use in pulses -- adds to Sfull with updatestimulus

	// random dots class
	CAnimation Rdots;

	char paramfilename[260]; // stores param file anme to be written out
	char worldfilename[260]; // stores world file name to be written out


	// stimulus data
	long int framenumber;   // counts which frame is going up
	LARGE_INTEGER CounterBeginTime;
	LARGE_INTEGER CounterFreq;
	FILE * monitor_position_file ;   // file to write to for positions
	FILE * monitor_stimulus_file ;   // file to write for stimulus info
	FILE * parameter_file ;          // parameters out, for a check
	FILE * error_file ;

	float GammaCorrect;
	float DistanceAccumulated;
	float FwdAccumulated;

	void readstring(FILE *f,char *str);
	float mod(float x,float y);

	float getDLPColor(float DLPintensity, char channel);

	struct sXYZ    // stores x,y, z coords read in from file...
	{
		float x, y, z;
	} World[1000];    // 1000 objects maximum in world...
	int NumObs;  // number of objects in data file

	struct MAPBUFFER
	{
		float x[500],y[500];
		int currpos;
	} Map;

	float backgroundtriple[3],foregroundtriple[3],blacktriple[3],meangroundtriple[3];
	void setColorTriples();
	float noiseValues[14400];  // get this on a workaround with periodic noise... can't make huge array!
	float retrieveNextNoiseValue();

	float calculateNextNoiseValue();

	RandomC r1, r2, r3, r4;

	EnforcedCorr1D RCorr;

	EnforcedCorr1D RCorr2D[8]; //use this to make grid!


	float vertexArray[120000];

public:
	
	//////////basic
	TrackStim();  // **
	virtual ~TrackStim(); // **

	//////////////set various things in the stimulus
	void setFlyMeasuredPos(float inx, float iny, float intheta, float intime); // ** sets the measured fly position
	void setFlyMeasuredDelta(float dx, float dy, float dtheta, float intime);
	void setFlyCurrentPos(float intime);  // ** sets the current fly position, linearly interpolating to current time...
	void updateStimCoords();
	void checkBoundaries();    // ** makes into torus -- add object avoidance here? what to do about that...
	void setZeroTime();        // ** also done with class first defined...
	void writePos();    // ** writes out measured position, do every measurement
	void writeStim();    // ** writes out stimulus data, do each color?
	void writeError();
	float queryCurrTime(); // ** gets current time
	float sineWave(float freq); // ** makes sine wave
	void incrementFrameNumber();
	STIMPARAMS initializeStimulus();
	VIEWPOSITIONS initializeViewPositions();
	TrackStim::STIMPARAMS * readParams(char szFile[260]);

	void setBackground();  // sets background color
	void setBackgroundColor();
	void setForeground();  // sets foreground color

	/////////////making real stimuli
	void readXYZ(char szFile[260]);   //**
	void readNoise();
	void rotationNoiseUpdate();
	void rotationNoiseCalcUpdate();
	void rotationXYBarCalcUpdate();
	void readViewPositions();  //**
	void writeViewPositions();   //**
	void readScreenPosition();  //**
	void writeScreenPosition();   //**
	void distancePulse();
	void rotationPulse();
	void rotationAccumulate();
	void rotationAccumulateWait();
	void rotationDoublePulsePlus();
	void rotationDoublePulseMinus();
	bool drawScene();  // ************************** executes the drawing set by setStimType
	void drawBackgroundCube();
	void drawWindow(float azim, float zenith, float solid);
	void drawMap();   //**
	void drawPosts();  //**
	void drawPosts2();  //**
	void drawAlignmentGrid();
	void drawCeiling();
	void drawUniform(float time);
	void drawWedge(GLfloat angle, GLfloat width);
	void drawWedge_mirror(GLfloat angle, GLfloat width);
	void drawWedge_horizontal(GLfloat angle, GLfloat width, GLfloat height, GLfloat added);
	void drawWedge_above0(GLfloat angle, GLfloat width);
	void drawWedge_below0(GLfloat angle, GLfloat width);
	void drawWedge_inside_top(GLfloat angle, GLfloat width);
	void drawWedge_inside_bottom(GLfloat angle, GLfloat width);
	void drawWedge_inside(GLfloat angle, GLfloat width);
	void drawSquare(GLfloat sqedge, GLfloat distance, GLfloat angle);
	void drawSawtooth(GLfloat angle, GLfloat width, GLfloat STphase, GLfloat NDphase);
	void drawSawtooth_setlum(GLfloat angle, GLfloat width, GLfloat STphase, GLfloat NDphase, GLfloat lumoffset);
	void drawWedgeSquare(GLfloat angle, GLfloat width, GLfloat height);
	void drawWedgeSquare_nobegin(GLfloat angle, GLfloat width, GLfloat height);
	void drawColorTest();
	void drawBar();
	void drawCylinder(float tc);
	void drawCylinder_new(float time);
	void drawCylinder_apertured(float time);
	void drawCylinder_continuous(float time);
	void drawExpandingBars(float tc);
	void drawExpandingBars_half(float tc);
	void drawExpandingBars_latinh_test(float tc);
	void drawExpandingBars_staircase(float tc);
	void drawExpandingBars_staircase0(float tc);
	void drawExpandingBars_adapting_competing(float tc);
	void drawExpandingBars_adapting_competing_edgelast(float tc);
	void drawExpandingBars_adapting_competing_edgelast_gray(float tc);
	void drawExpandingBars_competing_edges(float tc);
	void drawExpandingBars_competing_edges_translation(float tc);
	void drawExpandingBars_monocular(float tc);
	void drawExpandingBars_equal(float tc);
	void drawExpandingBars_fwd(float tc);
	void drawCylinderBars_wait(float tc);
	void drawCylinderBars_wait_eye_orientation(float tc);
	void drawCylinderBars_wait_bw(float tc);
	void drawCylinderBars_wait_monocular(float tc);
	void drawCylinderBars_wait_fwd(float tc);
	void drawBarSteps(float tc);
	void drawLoomAzimuthal(float tc); // SASKI FUNCTION!
	void drawMovingSquare(float tc); // SASKI FUNCTION!
	void drawMovingBar_je(float tc); //JE FUNCTION
	void drawSineCylinder(float tc);
	void drawSineCylinderWait(float tc);
	void drawSineCylinderWait_gray(float tc);
	void drawSineCylinderWait_2D(float tc);
	void drawSineCylinderWait_tanh(float tc);
	void drawScintillate(float tc);
	void drawScintillate_random_posneg(float tc);
	void drawScintillate_gradient_posneg(float tc);
	void drawScintillate_nograd_posneg(float tc);
	void drawScintillate_N2parity_posneg(float tc);
	void drawScintillate_N2parity_posneg_waittime(float tc);
	void drawScintillate_dtheta(float tc);
	void drawScintillate_complex_continuous(float tc);
	void drawScintillate_complex_dtheta(float tc);
	void drawScintillate_complex_dtheta_wait(float tc);
	void drawScintillate_complex_rand_dtheta(float tc);
	void drawScintillate_3pt_allchoice(float tc);
	void drawScint_equilum(float tc);
	void drawRasterRDotRot(float tc);
	void drawRasterRDotHemispheres(float tc);
	void drawBrightCylinder(float tc);
	void drawCheckerboard(float tc);
	void drawSingleBar(float timein);
	void drawDoubleBars(float timein);
	void drawDoubleBars_int(float timein);
	void drawDoubleBars_rand(float timein);
	void drawDoubleBars_rand_indep(float timein);
	void updateReichardtIntensity();
	void updateReichardtIntensity_binary();
	void updateReichardtIntensity_binary_isolating();
	void updateReichardtIntensity_countonly();
	void updateReichardtIntensity_discrete_gaussian();
	void updateReichCorr();
	void updateReichCorr_everytime();
	void updateReichCorrGrid_everytime();
	void drawReichCorrFullBars();
	void drawReichCorrFullGrid();
	void drawReichCorr2Bars_nocenter();
	void updateLumIntensity();
	void updateReichardtIntensity_binary_3bars();
	void drawReichardtSlits();
	void drawReichardtSlits_multiple();
	void drawReichardtSlits_multiple_nocenter();
	void drawReichardtSlits_multiple_nocenter_1static(float tc);
	void drawReichardtSlits_deterministic_steps(float tc);
	void drawReichardtSlits_multiple_nocenter_balanced();
	void drawReichardtSlits_multiple_oneeye();
	void drawReichardtSlits_multiple_dots(float tc);
	void drawReichardtSlits_multiple_widespace();
	void drawPerspectiveStrips(float angle, float width, float gap, bool LEFT); 
	void drawReichardtDots_nopersp(float angle, float width, float gap, bool LEFT, bool FIRST);
	void drawReichardtSlits_multiple_spacing();
	void drawReichardtSlits_multiple_binary();
	void drawReichardtSlits_multiple_binary_nocenter();
	void drawReichardtSlits_multiple_binary_nocenter_isolating();
	void drawReichardtSlits_original_gangster(float tc);
	void drawReichardtSlits_original_gangster_stagger(float tc);
	void drawReichardtSlits_TRANSLATE_original_gangster_stagger(float tc);
	void drawReichardtSlits_original_gangster_randomphase(float tc);
	void drawReichardtSlits_original_gangster_offcenter(float tc);
	void drawReichardtSlits_original_gangster_spaced(float tc);
	void drawReichardtSlits_original_gangster_4pt(float tc);
	void drawReichardtSlits_original_gangster_6pt(float tc);
	void drawReichardtSlits_original_gangster_6ptlong(float tc);
	void drawGrayWithVerticalNoise(float tc);
	void drawReichardtSlits_multiple_binary_nocenter_no45();
	void drawReichardtSlits_multiple_binary_nocenter_1eye();
	void drawReichardtSlits_multiple_binary_nocenter_1eye_quiet();
	void drawReichardtSlits_multiple_binary_nocenter_1eye_equalspace();
	void drawReichardtSlits_multiple_binary_nocenter_spacing();
	void drawReichardtSlits_multiple_binary_4corners();
	void drawReichardtSlits_multiple_binary_3bars();
	void drawReichardtSlits_multiple_translate_binary();
	void drawReichardtSlits_multiple_translate();
	void drawReichardtSlits_deterministic(float timein);
	void drawCylinderHalves(float time);
	void drawCylinderHalves_fwd(float time);
	void drawCylinderHemispheres(float time);
	void drawCylinderHemispheres_revised(float time);
	void drawCylinderHemispheres_revised_graycenter(float time);
	void drawCylinderHemispheres_revised_graycenter_wait(float time);
	void drawCylinderHalves_fwd_Sawtooth(float time);
	void drawCylinder_Sawtooth(float time);
	void drawCylinder_Sawtooth_randlum(float time);
	void drawCylinderHemispheresGain(float time);
	void drawInfCorridor(float currtime); // realistic translation
	void drawBarCeiling(float time);   // moving bar stimulus
	void drawLinearRDots(float dt); // random dots
	void drawRotRDots(float dt); // rotational random dots
	void drawRotRDotsHalves(float dt);
	void drawRotRDotsGrads(float dt);
	void drawSymmetricBars(float timein);
	void drawBar_je(float sqedge, float ledge, float distance, float angle); // jenn esch function
	void drawCenteringSquare(float tc); //jenn esch function
	void drawSquare_elevation(float sqedge, float elevation, float distance, float angle); //jenn esch
	void drawMovingSquare_elevation_je(float tc); //jenn esch
	void drawLoomAzimuthal_je(float tc); //jenn esch function added 03152013, adapted from Saskia's function
	void drawStarField(float tc); //alterd from drawScintilate -yf 131002
	void drawCylinderBarsES(float tc);
	void drawSimpleBarCW(float tc);
	void drawSimpleBarCCW(float tc);
	void drawSimpleBarCW_rphi(float tc);
	void drawSimpleBarCCW_rphi(float tc);

	// auxiliary functions for stimuli
	void updateStars(float tc, StarClass::Star star[]);
	void generateNextStar(StarClass::Star &star);
} ;