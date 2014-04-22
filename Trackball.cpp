// Trackball.cpp : Defines the entry point for the application.
//

#include "Trackball.h"


char * GetFileName()
{
	
	OPENFILENAME ofn;
	char szFile[260]; // Buffer for selected file name.
	szFile[0]=NULL;
	// blah, blah, more funtions...
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = hWnd;
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = "Text:\0*.txt\0";
	ofn.nFilterIndex = 0;
	ofn.lpstrFileTitle = NULL;
	ofn.nMaxFileTitle = 0;
	ofn.lpstrInitialDir = NULL;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

    if (GetOpenFileName(&ofn))
		return ofn.lpstrFile;
	else
		return NULL;

};


void readstr(FILE *f,char *str)
{
	do
	{
		fgets(str, 255, f);
	} while ((str[0] == '/') || (str[0] == '\n'));
	return;
}

void SetupTracks()
{
	float x, y, theta;
	float arenasize=20;
	int maxtracks;
	FILE * filein;
	char oneline[255];
	for (int ii=0; ii<255; ii++) oneline[ii]='\0';
	filein = fopen("data/xytheta.txt", "rt");				// File To Load World Data From

	readstr(filein,oneline);
	sscanf(oneline, "NUMTRACKS %d\n", &maxtracks);

	for (int loop = 0; loop < maxtracks; loop++)
	{
		readstr(filein,oneline);
		sscanf(oneline, "%f %f %f", &x, &y, &theta);
		flyxmat[loop] = x*25;
		flyymat[loop] = y*25;
		flythetamat[loop] = theta/80;
	}
	fclose(filein);
	return;
}

void assignMouseCalcs(float *ldx,float *ldy,float *ldtheta, float currang)
{
	// everything to be output to mm, cause i say so
	// direction is x is forward, right handed coord system
	// only works when dmousex and y and small so that linear approximations work
	//////float DPI = 675; // measured DPI of mice // DAC 080909 -- 675 -- measured 360 rotation...
	//////float DPI = 625; // measured DAC 081120, from 180 deg rotation
	//////float DPI = 700; // measured DAC 110530, from 180 deg rotations, backed it out = pi*d/2
	float DPI = 1146.9;
	//////float mouse1factor = 1.15; // this makes them match in extent...
	float mouse1factor = 1.01; // match very closely with current alignment, 110530 DAC
	float ballR = 0.125;  // inches // DAC 080909 -- 0.125 (correct)
	float dang, fwd, lat;  // all those things in local frame of reference
	float factor=2.0f/sqrt(3.0f);
	dang = -1 * ((float)(long)dmousex[0]*mouse1factor + (float)(long)dmousex[1])/2    / DPI / ballR;   // radians
	// -1 added below on 110701, end of day; all data before that has inverted fwd...
	fwd = -1*((float)(long)dmousey[0]*mouse1factor + (float)(long)dmousey[1])/2 * factor     / DPI * 25.4;  // mm
	lat = ((float)(long)dmousey[0]*mouse1factor - (float)(long)dmousey[1])/2 * factor     / DPI * 25.4;  // mm -- check sign on this one to match theta

	*ldx = fwd * cos(currang) - lat * sin(currang); // currang is integrated angle
	*ldy = fwd * sin(currang) + lat * cos(currang);
	*ldtheta = dang; 
};

void writeMouseCoords()
{
	static FILE * mouse_position_file = NULL; // file pointer
	static int loop = 0;   // loop number
	if (mouse_position_file == NULL)
	{
		mouse_position_file=fopen("_mouse_output.txt","a");
		
		time_t rawtime;
		struct tm * timeinfo;
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		fprintf (mouse_position_file, asctime (timeinfo) );
		
	}
	else
	{	
		fprintf(mouse_position_file,"%8d %8d %8d %8d %8d \n",loop,dmousex[0],dmousey[0],dmousex[1],dmousey[1]);
	};
	loop++;
};

void writeMainSetup()
{
	static FILE * main_booleans_file = NULL;
	// one time printing...
	main_booleans_file = fopen("_main_booleans.txt","a");
	fprintf (main_booleans_file, "useDLP %i\n", useDLP);
	fprintf (main_booleans_file, "USEPARAMS %i\n", USEPARAMS);
	fprintf (main_booleans_file, "USEMOUSE %i\n", USEMOUSE);
	fprintf (main_booleans_file, "USETRACKS %i\n", USETRACKS);
	fprintf (main_booleans_file, "MAXRUNTIME %f\n", MAXRUNTIME);
	fprintf (main_booleans_file, "LOADBASICS %i\n", LOADBASICS);
	fprintf (main_booleans_file, "VIEWCUBE %i\n", VIEWCUBE);
	fprintf (main_booleans_file, "FRUSTUM_PERSECTIVE %i\n", FRUSTUM_PERSPECTIVE);
	
	fclose (main_booleans_file);
};

int ShowHideTaskBar ( BOOL bHide )
{
	// find taskbar window
	CWnd* pWnd = CWnd::FindWindow(_T("Shell_TrayWnd"), _T(""));
	if(!pWnd )
		return 0;
	if( bHide )
	{
		// hide taskbar
		pWnd->ShowWindow(SW_HIDE);
	}
	else
	{
		// show taskbar
		pWnd->ShowWindow(SW_SHOW);
	}
	return 1;
}

GLvoid ReSizeGLScene(GLsizei width, GLsizei height)		// Resize And Initialize The GL Window
{
	// PRETTY MUCH OBSOLETE, AS I DON'T RESIZE WINDOW EVER, REALLY ALLOW FOR THIS EVEN... LEGACY
	if (height==0)										// Prevent A Divide By Zero By
	{
		height=1;										// Making Height Equal One
	}

	glViewport(0,0,width,height);						// Reset The Current Viewport

	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	gluPerspective(90.0f,(GLfloat)width/(GLfloat)height,0.1f,100.0f);

	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
}

int InitGL(GLvoid)										// All Setup For OpenGL Goes Here
{

	//glEnable(GL_TEXTURE_2D);							// Enable Texture Mapping
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE);					// Set The Blending Function For Translucency
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);				// This Will Clear The Background Color To Black
	glClear(GL_COLOR_BUFFER_BIT);
	glClearDepth(1.0);									// Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LESS);								// The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST);							// Enables Depth Testing
	//glShadeModel(GL_SMOOTH);							// Enables Smooth Color Shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);	// Really Nice Perspective Calculations

	// anti-aliasing stuff? DAC 080827
	glEnable (GL_BLEND); // needed!
	//glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);			// Set The Blending Function For Translucency
	//glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA); // not great -- get the seams showing up
	glBlendFunc(GL_ONE,GL_ZERO);                         // aliased, but no shadows -- one solution to shadows is to make the clear color the background color again
	//glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	//glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);  // can be GL_NICEST, GL_FASTEST, GL_DONT_CARE
	//glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	// SetupPosts();

	// fileout = fopen("data/monitor.txt", "w");

	return TRUE;										// Initialization Went OK
}



GLvoid KillGLWindow(GLvoid)								// Properly Kill The Window
{
	if (fullscreen)										// Are We In Fullscreen Mode?
	{
		ChangeDisplaySettings(NULL,0);					// If So Switch Back To The Desktop
		ShowCursor(TRUE);								// Show Mouse Pointer
	}

	ShowHideTaskBar(FALSE);  // show taskbar when done...

	if (hRC)											// Do We Have A Rendering Context?
	{
		if (!wglMakeCurrent(NULL,NULL))					// Are We Able To Release The DC And RC Contexts?
		{
			MessageBox(NULL,"Release Of DC And RC Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}

		if (!wglDeleteContext(hRC))						// Are We Able To Delete The RC?
		{
			MessageBox(NULL,"Release Rendering Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		}
		hRC=NULL;										// Set RC To NULL
	}

	if (hDC && !ReleaseDC(hWnd,hDC))					// Are We Able To Release The DC
	{
		MessageBox(NULL,"Release Device Context Failed.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hDC=NULL;										// Set DC To NULL
	}

	if (hWnd && !DestroyWindow(hWnd))					// Are We Able To Destroy The Window?
	{
		MessageBox(NULL,"Could Not Release hWnd.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hWnd=NULL;										// Set hWnd To NULL
	}

	if (!UnregisterClass("OpenGL",hInstance))			// Are We Able To Unregister Class
	{
		MessageBox(NULL,"Could Not Unregister Class.","SHUTDOWN ERROR",MB_OK | MB_ICONINFORMATION);
		hInstance=NULL;									// Set hInstance To NULL
	}
}

/*	This Code Creates Our OpenGL Window.  Parameters Are:					*
 *	title			- Title To Appear At The Top Of The Window				*
 *	width			- Width Of The GL Window Or Fullscreen Mode				*
 *	height			- Height Of The GL Window Or Fullscreen Mode			*
 *	bits			- Number Of Bits To Use For Color (8/16/24/32)			*
 *	fullscreenflag	- Use Fullscreen Mode (TRUE) Or Windowed Mode (FALSE)	*/
 
BOOL CreateGLWindow(char* title, int width, int height, int bits, bool fullscreenflag)
{
	GLuint		PixelFormat;			// Holds The Results After Searching For A Match
	WNDCLASS	wc;						// Windows Class Structure
	DWORD		dwExStyle;				// Window Extended Style
	DWORD		dwStyle;				// Window Style
	RECT		WindowRect;				// Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left=(long)0;			// Set Left Value To 0
	WindowRect.right=(long)width;		// Set Right Value To Requested Width
	WindowRect.top=(long)0;				// Set Top Value To 0
	WindowRect.bottom=(long)height;		// Set Bottom Value To Requested Height

	fullscreen=fullscreenflag;			// Set The Global Fullscreen Flag

	hInstance			= GetModuleHandle(NULL);				// Grab An Instance For Our Window
	wc.style			= CS_HREDRAW | CS_VREDRAW | CS_OWNDC;	// Redraw On Size, And Own DC For Window.
	wc.lpfnWndProc		= (WNDPROC) WndProc;					// WndProc Handles Messages
	wc.cbClsExtra		= 0;									// No Extra Window Data
	wc.cbWndExtra		= 0;									// No Extra Window Data
	wc.hInstance		= hInstance;							// Set The Instance
	wc.hIcon			= LoadIcon(NULL, IDI_WINLOGO);			// Load The Default Icon
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);			// Load The Arrow Pointer
	wc.hbrBackground	= NULL;									// No Background Required For GL
	wc.lpszMenuName		= NULL;									// We Don't Want A Menu
	wc.lpszClassName	= "OpenGL";								// Set The Class Name

	if (!RegisterClass(&wc))									// Attempt To Register The Window Class
	{
		MessageBox(NULL,"Failed To Register The Window Class.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;											// Return FALSE
	}
	
	if (fullscreen)												// Attempt Fullscreen Mode?
	{
		DEVMODE dmScreenSettings;								// Device Mode
		memset(&dmScreenSettings,0,sizeof(dmScreenSettings));	// Makes Sure Memory's Cleared
		dmScreenSettings.dmSize=sizeof(dmScreenSettings);		// Size Of The Devmode Structure
		dmScreenSettings.dmPelsWidth	= width;				// Selected Screen Width
		dmScreenSettings.dmPelsHeight	= height;				// Selected Screen Height
		dmScreenSettings.dmBitsPerPel	= bits;					// Selected Bits Per Pixel
		dmScreenSettings.dmDisplayFrequency = 120;              // DAC 080819 -- see if this works?
		dmScreenSettings.dmFields=DM_BITSPERPEL|DM_PELSWIDTH|DM_PELSHEIGHT;

		// Try To Set Selected Mode And Get Results.  NOTE: CDS_FULLSCREEN Gets Rid Of Start Bar.
		if (ChangeDisplaySettings(&dmScreenSettings,CDS_FULLSCREEN)!=DISP_CHANGE_SUCCESSFUL)
		{
			// If The Mode Fails, Offer Two Options.  Quit Or Use Windowed Mode.
			if (MessageBox(NULL,"The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?","NeHe GL",MB_YESNO|MB_ICONEXCLAMATION)==IDYES)
			{
				fullscreen=FALSE;		// Windowed Mode Selected.  Fullscreen = FALSE
			}
			else
			{
				// Pop Up A Message Box Letting User Know The Program Is Closing.
				MessageBox(NULL,"Program Will Now Close.","ERROR",MB_OK|MB_ICONSTOP);
				return FALSE;									// Return FALSE
			}
		}
	}

	if (useDLP)
	{
		if (!ShowHideTaskBar(TRUE)) // to hide the task bar!
			MessageBox(NULL,"Could not hide the taskbar","Taskbar error",MB_OK | MB_ICONINFORMATION);
	};

	if (fullscreen)												// Are We Still In Fullscreen Mode?
	{
		dwExStyle=WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle=WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	}
	else
	{
		dwExStyle=WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;			// Window Extended Style
		dwStyle=WS_OVERLAPPEDWINDOW;							// Windows Style
		ShowCursor(FALSE);                                      // DAC 080819
	}

	if (useDLP)
	{
		dwExStyle=WS_EX_APPWINDOW;								// Window Extended Style
		dwStyle=WS_POPUP;										// Windows Style
		ShowCursor(FALSE);										// Hide Mouse Pointer
	};

	AdjustWindowRectEx(&WindowRect, dwStyle, FALSE, dwExStyle);		// Adjust Window To True Requested Size

	// Create The Window
	if (!useDLP)
	{
		if (!(hWnd=CreateWindowEx(	dwExStyle,							// Extended Style For The Window
									"OpenGL",							// Class Name
									title,								// Window Title
									dwStyle |							// Defined Window Style
									WS_CLIPSIBLINGS |					// Required Window Style
									WS_CLIPCHILDREN,					// Required Window Style
									0, 0,								// Window Position -- can modify this! DAC
									WindowRect.right-WindowRect.left,	// Calculate Window Width
									WindowRect.bottom-WindowRect.top,	// Calculate Window Height
									NULL,								// No Parent Window
									NULL,								// No Menu
									hInstance,							// Instance
									NULL)))								// Dont Pass Anything To WM_CREATE
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL,"Window Creation Error.","ERROR",MB_OK|MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		};
	}
	else
	{
		if (!(hWnd=CreateWindowEx(	dwExStyle,							// Extended Style For The Window
									"OpenGL",							// Class Name
									title,								// Window Title
									//dwStyle |							// Defined Window Style
									//WS_CLIPSIBLINGS |					// Required Window Style
									//WS_CLIPCHILDREN | 
									WS_BORDER | WS_POPUP,					// Required Window Style
									1680, 0,								// Window Position -- can modify this! DAC
									WindowRect.right-WindowRect.left,	// Calculate Window Width
									WindowRect.bottom-WindowRect.top,	// Calculate Window Height
									NULL,								// No Parent Window
									NULL,								// No Menu
									hInstance,							// Instance
									NULL)))								// Dont Pass Anything To WM_CREATE
		{
			KillGLWindow();								// Reset The Display
			MessageBox(NULL,"Window Creation Error.","ERROR",MB_OK|MB_ICONEXCLAMATION);
			return FALSE;								// Return FALSE
		};
	}
	

	static	PIXELFORMATDESCRIPTOR pfd=				// pfd Tells Windows How We Want Things To Be
	{
		sizeof(PIXELFORMATDESCRIPTOR),				// Size Of This Pixel Format Descriptor
		1,											// Version Number
		PFD_DRAW_TO_WINDOW |						// Format Must Support Window
		PFD_SUPPORT_OPENGL |						// Format Must Support OpenGL
		PFD_DOUBLEBUFFER,							// Must Support Double Buffering
		PFD_TYPE_RGBA,								// Request An RGBA Format
		bits,										// Select Our Color Depth
		0, 0, 0, 0, 0, 0,							// Color Bits Ignored
		0,											// No Alpha Buffer
		0,											// Shift Bit Ignored
		0,											// No Accumulation Buffer
		0, 0, 0, 0,									// Accumulation Bits Ignored
		16,											// 16Bit Z-Buffer (Depth Buffer)  
		0,											// No Stencil Buffer
		0,											// No Auxiliary Buffer
		PFD_MAIN_PLANE,								// Main Drawing Layer
		0,											// Reserved
		0, 0, 0										// Layer Masks Ignored
	};
	
	if (!(hDC=GetDC(hWnd)))							// Did We Get A Device Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Device Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(PixelFormat=ChoosePixelFormat(hDC,&pfd)))	// Did Windows Find A Matching Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Find A Suitable PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!SetPixelFormat(hDC,PixelFormat,&pfd))		// Are We Able To Set The Pixel Format?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Set The PixelFormat.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if (!(hRC=wglCreateContext(hDC)))				// Are We Able To Get A Rendering Context?
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Create A GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	if(!wglMakeCurrent(hDC,hRC))					// Try To Activate The Rendering Context
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Can't Activate The GL Rendering Context.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	ShowWindow(hWnd,SW_SHOW);						// Show The Window
	SetForegroundWindow(hWnd);						// Slightly Higher Priority
	SetFocus(hWnd);									// Sets Keyboard Focus To The Window
	ReSizeGLScene(width, height);					// Set Up Our Perspective GL Screen

	if (!InitGL())									// Initialize Our Newly Created GL Window
	{
		KillGLWindow();								// Reset The Display
		MessageBox(NULL,"Initialization Failed.","ERROR",MB_OK|MB_ICONEXCLAMATION);
		return FALSE;								// Return FALSE
	}

	return TRUE;									// Success
}

LRESULT CALLBACK WndProc(	HWND	hWnd,			// Handle For This Window
							UINT	uMsg,			// Message For This Window
							WPARAM	wParam,			// Additional Message Information
							LPARAM	lParam)			// Additional Message Information
{
	switch (uMsg)									// Check For Windows Messages
	{
		case WM_ACTIVATE:							// Watch For Window Activate Message
		{
			if (!HIWORD(wParam))					// Check Minimization State
			{
				active=TRUE;						// Program Is Active
			}
			else
			{
				active=FALSE;						// Program Is No Longer Active
			}

			return 0;								// Return To The Message Loop
		}

		case WM_SYSCOMMAND:							// Intercept System Commands
		{
			switch (wParam)							// Check System Calls
			{
				case SC_SCREENSAVE:					// Screensaver Trying To Start?
				case SC_MONITORPOWER:				// Monitor Trying To Enter Powersave?
				return 0;							// Prevent From Happening
			}
			break;									// Exit
		}

		case WM_CLOSE:								// Did We Receive A Close Message?
		{
			PostQuitMessage(0);						// Send A Quit Message
			return 0;								// Jump Back
		}

		case WM_KEYDOWN:							// Is A Key Being Held Down?
		{
			keys[wParam] = TRUE;					// If So, Mark It As TRUE
			return 0;								// Jump Back
		}

		case WM_KEYUP:								// Has A Key Been Released?
		{
			keys[wParam] = FALSE;					// If So, Mark It As FALSE
			return 0;								// Jump Back
		}

		case WM_SIZE:								// Resize The OpenGL Window
		{
			ReSizeGLScene(LOWORD(lParam),HIWORD(lParam));  // LoWord=Width, HiWord=Height
			//if (!fileout == NULL)
			//	fprintf (fileout,"WM_SIZE called\n");
			return 0;								// Jump Back
		}

		case WM_INPUT:    // accumulates stuff properly...
		{
			add_to_raw_mouse_x_and_y((HRAWINPUT)lParam);
			//get_xypos(&mousex,&mousey); // obsolete now
			//fprintf (fileout, "WM_INPUT called case\n");
			return 0;
		}
	}

	// Pass All Unhandled Messages To DefWindowProc
	return DefWindowProc(hWnd,uMsg,wParam,lParam);
}

int WINAPI WinMain(	HINSTANCE	hInstance,			// Instance
					HINSTANCE	hPrevInstance,		// Previous Instance
					LPSTR		lpCmdLine,			// Command Line Parameters
					int			nCmdShow)			// Window Show State
{
	MSG		msg;									// Windows Message Structure
	BOOL	done=FALSE;								// Bool Variable To Exit Loop
	char * paramfilename, paramfilekeep[260], * worldfilename, worldfilekeep[260];    // both will be passed to trackstim and used to write things out...
	TrackStim::STIMPARAMS AllStims[40];
	TrackStim::STIMPARAMS * AllStimsTemp;
	float	currtime;
	
	srand(time(0));

	// Ask The User Which Screen Mode They Prefer
	/*if (MessageBox(NULL,"Would You Like To Run In Fullscreen Mode?", "Start FullScreen?",MB_YESNO|MB_ICONQUESTION)==IDNO)
	{
		fullscreen=FALSE;							// Windowed Mode
	};*/

	fullscreen=FALSE;
	if (MessageBox(NULL,"Put it onto the DLP", "Start FullScreen?",MB_YESNO|MB_ICONQUESTION)==IDNO)
	{
		useDLP=FALSE;							// Windowed Mode
	};

	// set up stimulus here...
	TrackStim Onscreen;   // instantiate class

	if (USETRACKS) SetupTracks();

	if (!useDLP) Onscreen.ScreenUpdateRate=60;  // we're doing it on the monitor, so adjust; otherwise, it initializes to 120 for DLP
	 
	if (LOADBASICS)
	{
		Onscreen.readNoise();
		
		// read in world data...
		//worldfilename = GetFileName() ;
		//for (int ii=0; ii<260; ii++) worldfilekeep[ii]=worldfilename[ii]; // klugey
		//Onscreen.readXYZ(worldfilekeep);
		Onscreen.readXYZ("E:/Code/trackball_final/Trackball/Debug/Data/posts.txt");
	};

	// read in parameter data... 
	paramfilename = GetFileName() ;
	for (int ii=0; ii<260; ii++) paramfilekeep[ii]=paramfilename[ii];
	AllStimsTemp = Onscreen.readParams(paramfilekeep);
	if (USEPARAMS)
		for (int jj=0; jj<Onscreen.numepochs; jj++) AllStims[jj]=AllStimsTemp[jj];
	else
		for (int jj=0; jj<Onscreen.numepochs; jj++) AllStims[jj]=Onscreen.initializeStimulus();  // set them all to the initialized one
	Onscreen.epochchoose=0;
	Onscreen.Stimulus = AllStims[Onscreen.epochchoose]; 

	// read in viewport data & screen positioning data
	Onscreen.readViewPositions(); // if no file exists, just initializes it...
	Onscreen.readScreenPosition();

	VIEWCUBE = Onscreen.VIEWCUBECORNER;
	FRUSTUM_PERSPECTIVE = Onscreen.USE_FRUSTUM;

	writeMainSetup();  // writes out file with set up of main...

	// Create Our OpenGL Window
	if (!CreateGLWindow("TrackBall",800,600,16,fullscreen)) // check this, might want 24 bits/pixel...
	{
		return 0;									// Quit If Window Was Not Created
	}
	
	// initialize mice...
	if (!init_raw_mouse(1, 0, 1)) { // registers for (sysmouse=yes,  terminal services mouse=no, HID_mice=yes)
	  MessageBox(NULL, "RawInput not supported by Operating System.  Exiting." , "Error!" ,MB_OK);
	  return 0; 
	}

	Onscreen.setZeroTime(); // set time to be zero
	Onscreen.epochtimezero=Onscreen.queryCurrTime();    // zero it...


	while(!done)									// Loop That Runs While done=FALSE
	{
		if (PeekMessage(&msg,NULL,0,0,PM_REMOVE))	// Is There A Message Waiting?
		{
			if (msg.message==WM_QUIT)				// Have We Received A Quit Message?
			{
				done=TRUE;							// If So done=TRUE
			}
			else									// If Not, Deal With Window Messages
			{
				TranslateMessage(&msg);				// Translate The Message
				DispatchMessage(&msg);				// Dispatch The Message
			}
		}
		else										// If There Are No Messages
		{
			// Draw The Scene.  Watch For ESC Key And Quit Messages From DrawGLScene()
			//if ((active && !Onscreen.drawScene()) || keys[VK_ESCAPE])	// not most recent mouse commands in draw?
			if (keys[VK_ESCAPE]) 
			{
				done=TRUE;							// ESC or DrawGLScene Signalled A Quit
			}
			else									// Not Time To Quit, Update Screen
			{
				// do all the good stuff here, updating things, looking for things, etc., including mouse update here...
				currtime = Onscreen.queryCurrTime();

				// control of program here...
				//if (FALSE)   
				if ((currtime - Onscreen.epochtimezero > Onscreen.Stimulus.duration) && (Onscreen.numepochs > 1))
				{
					if (Onscreen.Stimulus.changemethod == 0) // randomize it, interleaved with zero
						if (Onscreen.epochchoose == 0)
							Onscreen.epochchoose = ( rand() % (Onscreen.numepochs - 1) ) + 1;  // from 1 .. numepochs - 1
						else
							Onscreen.epochchoose = 0; // set back to background...

					if (Onscreen.Stimulus.changemethod == 1) // go in order
						Onscreen.epochchoose = (Onscreen.epochchoose + 1) % Onscreen.numepochs;

					if (Onscreen.Stimulus.changemethod == 2) // randomize all, no interleave
						Onscreen.epochchoose = (rand() % Onscreen.numepochs);

					Onscreen.Stimulus = AllStims[Onscreen.epochchoose];
					Onscreen.epochtimezero=currtime;   // switch them up
				}

				if (USEMOUSE)
				{
					dx = (float)(long)get_raw_mouse_x_delta(0); // something to do with 64-bit conversions. weird. but if it works, use it, even if an indication of my shitty coding
					dy = (float)(long)get_raw_mouse_y_delta(0);
				} 
				else
				{
					dmousex[0]=get_raw_mouse_x_delta(MOUSE_LEFT);
					dmousey[0]=get_raw_mouse_y_delta(MOUSE_LEFT);
					dmousex[1]=get_raw_mouse_x_delta(MOUSE_RIGHT);
					dmousey[1]=get_raw_mouse_y_delta(MOUSE_RIGHT);

					assignMouseCalcs(&dx,&dy,&dtheta,Onscreen.Fmeas.theta); // procedure to assign mouse calcs from dmousex/y
					// feed these to function to get dx, dy, and dtheta
				};
				
				if (!USETRACKS)
				{
					Onscreen.setFlyMeasuredDelta(xpos + dx, ypos + dy, zrot + dtheta, currtime);
					xpos=0; ypos=0; zrot=0;  // so that they accumulate only the delta x
				}
				else
				{
					Onscreen.setFlyMeasuredPos(flyxmat[(flytracknum ) % sizeof(flyxmat)/sizeof(float)],
							flyymat[(flytracknum ) % sizeof(flyxmat)/sizeof(float)],
							flythetamat[(flytracknum ) % sizeof(flyxmat)/sizeof(float)],
							currtime);
					flytracknum++; 
				}

				// update stuff now so that writes are correct new values... not needed -- done in drawscene...
				//Onscreen.setFlyCurrentPos(currtime);
				//Onscreen.updateStimCoords();

				
				Onscreen.drawScene();  // has own Fcurr and Sfull update lines...
				SwapBuffers(hDC);					// Swap Buffers (Double Buffering)
				
				if (WRITEOUT) { // do after drawscene to incorporate functions that modify the Sfull...
					Onscreen.writePos();
					Onscreen.writeStim();
					if (useDLP)
						Onscreen.writeError();
					writeMouseCoords();
				};
				Onscreen.incrementFrameNumber(); // important to keep track of some stuff...

				// peace out...
				if (currtime > MAXRUNTIME)
					done = TRUE;

				if (keys[' ']) WRITEOUT = TRUE; // starts as false, space bar to make write out...

				if (keys['J'])
				{
					zrot -= 0.05f;
				}

				if (keys['K'])
				{
					zrot += 0.05f;
				}

				if (keys[VK_PRIOR])
				{
				}

				if (keys[VK_NEXT])
				{
				}

				if (keys[VK_UP] && (!(keys['Q'] || keys['W'] || keys['E'])))
				{
					ypos += 0.1f;
				}

				if (keys[VK_DOWN] && (!(keys['Q'] || keys['W'] || keys['E'])))
				{
					ypos -= 0.1f;
				}

				if (keys[VK_RIGHT] && (!(keys['Q'] || keys['W'] || keys['E'])))
				{
					xpos += 0.1f;
				}

				if (keys[VK_LEFT] && (!(keys['Q'] || keys['W'] || keys['E'])))
				{
					xpos -= 0.1f;
				}

				if (keys['X']) {Onscreen.writeScreenPosition();}; // write them out

				if (keys['Z']) {Onscreen.readScreenPosition();};  // read them in -- file must exist in order not to crash!

				if (keys['E'])
				{
					if (keys['2'])   // right window!
					{
						if (keys[VK_UP]) {Onscreen.ScreenPosition.y[0] += 0.01;};
						if (keys[VK_DOWN]) {Onscreen.ScreenPosition.y[0] -= 0.01;};
						if (keys[VK_RIGHT]) {Onscreen.ScreenPosition.x[0] += 0.01;};
						if (keys[VK_LEFT]) {Onscreen.ScreenPosition.x[0] -= 0.01;};
						if (keys[VK_PRIOR]) {Onscreen.ScreenPosition.z[0] += 0.01;};
						if (keys[VK_NEXT]) {Onscreen.ScreenPosition.z[0] -= 0.01;};						
					};
					if (keys['4'])     // left window!
					{
						if (keys[VK_UP]) {Onscreen.ScreenPosition.y[1] += 0.01;};
						if (keys[VK_DOWN]) {Onscreen.ScreenPosition.y[1] -= 0.01;};
						if (keys[VK_RIGHT]) {Onscreen.ScreenPosition.x[1] += 0.01;};
						if (keys[VK_LEFT]) {Onscreen.ScreenPosition.x[1] -= 0.01;};
						if (keys[VK_PRIOR]) {Onscreen.ScreenPosition.z[1] += 0.01;};
						if (keys[VK_NEXT]) {Onscreen.ScreenPosition.z[1] -= 0.01;};						
					};
				};

				if (keys['S']) {Onscreen.writeViewPositions();}; // write them out

				if (keys['A']) {Onscreen.readViewPositions();};  // read them in -- file must exist in order not to crash!
				
				if (keys['Q'])  // window adjust -- make it difficult!
				{
					if (keys['2'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.y[0]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.y[0]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.x[0]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.x[0]--;};
					}
					if (keys['3'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.y[1]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.y[1]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.x[1]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.x[1]--;};
					}
					if (keys['4'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.y[2]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.y[2]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.x[2]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.x[2]--;};
					}
					if (keys['5'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.y[3]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.y[3]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.x[3]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.x[3]--;};
					}
				}

				if (keys['W'])  // window adjust -- make it difficult!
				{
					if (keys['2'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.h[0]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.h[0]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.w[0]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.w[0]--;};
					}
					if (keys['3'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.h[1]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.h[1]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.w[1]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.w[1]--;};
					}
					if (keys['4'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.h[2]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.h[2]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.w[2]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.w[2]--;};
					}
					if (keys['5'])
					{
						if (keys[VK_UP]) {Onscreen.ViewPorts.h[3]++;};
						if (keys[VK_DOWN]) {Onscreen.ViewPorts.h[3]--;};
						if (keys[VK_RIGHT]){Onscreen.ViewPorts.w[3]++;};
						if (keys[VK_LEFT]) {Onscreen.ViewPorts.w[3]--;};
					}
				}

				if (keys[VK_F1])						// Is F1 Being Pressed?
				{
					keys[VK_F1]=FALSE;					// If So Make Key FALSE
					KillGLWindow();						// Kill Our Current Window
					fullscreen=!fullscreen;				// Toggle Fullscreen / Windowed Mode
					// Recreate Our OpenGL Window
					if (!CreateGLWindow("Checkerboard",800,600,16,fullscreen))
					{
						return 0;						// Quit If Window Was Not Created
					}
				}
			}
		}
	}

	// Shutdown
	KillGLWindow();										// Kill The Window
	return (msg.wParam);								// Exit The Program
}