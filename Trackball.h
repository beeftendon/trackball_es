#pragma once

// #include "resource.h"

#define _WIN32_WINNT 0x501
#define WM_INPUT 0x00FF
#define MOUSE_LEFT 1
#define MOUSE_RIGHT 2


#include "stdafx.h"
#include <afxwin.h>
//#include <windows.h>
#include <math.h>			// Math Library Header File
#include <stdio.h>			// Header File For Standard Input/Output
#include <gl\gl.h>			// Header File For The OpenGL32 Library
#include <gl\glu.h>			// Header File For The GLu32 Library
//#include <gl\glaux.h>		// Header File For The Glaux Library
#include "rawmouse.h"
#include "TrackStim.h"
#include <commdlg.h>

HDC			hDC=NULL;		// Private GDI Device Context
HGLRC		hRC=NULL;		// Permanent Rendering Context
HWND		hWnd=NULL;		// Holds Our Window Handle
HINSTANCE	hInstance;		// Holds The Instance Of The Application

bool	keys[256];			// Array Used For The Keyboard Routine
bool	active=TRUE;		// Window Active Flag Set To TRUE By Default
bool	fullscreen=TRUE;	// Fullscreen Flag Set To Fullscreen Mode By Default
bool	useDLP=TRUE;
bool	USEPARAMS= TRUE;
bool	LOADBASICS = FALSE;
bool	VIEWCUBE = FALSE;
float	MAXRUNTIME=1200;   // 20 minutes... 
bool	FRUSTUM_PERSPECTIVE = TRUE;

const float piover180 = 0.0174532925f;
float heading;
float xpos;
float zpos;
float ypos;

float	zrot;

int dmousex[2]; // x coords for first and second mice
int dmousey[2]; // y coords
bool USEMOUSE = FALSE;
bool USETRACKS = FALSE;
bool WRITEOUT = TRUE;
int flytracknum=0;
float flyxmat[10000];
float flyymat[10000];
float flythetamat[10000];

float dx=0;
float dy=0;
float dtheta=0;  // total dx, dy, and dtheta calculated from mousex,mousey, or direct from mouse
float dfwd=0;
float dlat=0;


LRESULT	CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);	// Declaration For WndProc
