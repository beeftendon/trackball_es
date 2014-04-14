#include "stdafx.h"
#include "EnforcedCorr1D.h"

/////////////// functions to set and write stimulus stuff from fly measurements //////////////////////////////

EnforcedCorr1D::EnforcedCorr1D()
{

	Pwhite = 0.5f;
	numberofpixels=360;
	for (int ii=0; ii<360; ii++)
	{
		pixvals[ii]=0;
		pixvalsnew[ii]=0;
		parityvals[0][ii]=-1;
		parityvals[1][ii]=2*(rand() % 2)-1;
		parityvals[2][ii]=2*(rand() % 2)-1;
		parityvals[3][ii]=2*(rand() % 2)-1;
		parityvals[4][ii]=2*(rand() % 2)-1;
		parityvals[5][ii]=2*(rand() % 2)-1;
		parityvals[6][ii]=2*(rand() % 2)-1;
		parityvals[7][ii]=2*(rand() % 2)-1;
		parityvals[8][ii]=2*(rand() % 2)-1;
		parityvals[9][ii]=2*(rand() % 2)-1;
	}
	updatetype = 0;
	circularBC = FALSE;

	direction = 1;
	parity = 1;
	parityNUM = 2;

	select_frames = 1;

} ;

EnforcedCorr1D::~EnforcedCorr1D()
{
};

void EnforcedCorr1D::updateparityvals()
{
	for (int ii=0; ii<numberofpixels; ii++)
	{
		parityvals[9][ii]=parityvals[8][ii];
		parityvals[8][ii]=parityvals[7][ii];
		parityvals[7][ii]=parityvals[6][ii];
		parityvals[6][ii]=parityvals[5][ii];
		parityvals[5][ii]=parityvals[4][ii];
		parityvals[4][ii]=parityvals[3][ii];
		parityvals[3][ii]=parityvals[2][ii];
		parityvals[2][ii]=parityvals[1][ii];
		parityvals[1][ii]=parityvals[0][ii];
		if (pixvals[ii]==1)
			parityvals[0][ii]=1;
		else
			parityvals[0][ii]=-1;
	};
	PIXUPDATED=TRUE;
};

void EnforcedCorr1D::setnumberofpixels(int numin)
{
	numberofpixels=numin;
};

void EnforcedCorr1D::setPwhite(float Pin)
{
	Pwhite = Pin;
};

void EnforcedCorr1D::setupdatetype(int updatein)
{
	updatetype = updatein;
};

void EnforcedCorr1D::setBC(bool circBC)
{
	circularBC = circBC;
};

void EnforcedCorr1D::setparity(int parityin)
{
	parity = parityin;
};

void EnforcedCorr1D::setparityNUM(int parityNUMin)
{
	parityNUM = parityNUMin;
};

void EnforcedCorr1D::updatepixels()
{
	int xoff[10];
	int toff[10];
	switch (updatetype)
	{
	case 0:
		randomizepixels();
		break;
	case 1:
		update1(1,1);
		break;
	case 2:
		update1(0,0);
		break;
	case 3:
		update1(1,0);
		break;
	case 4:
		update1(0,1);
		break;
	case 5:
		update1corr(1,1,0.5);
		break;
	case 6:
		update1corr(0,0,0.5);
		break;
	case 7:
		update1corr(1,0,0.5);
		break;
	case 8:
		update1corr(0,1,0.5);
		break;
	case 9:
		updateparitypixvals(parity,parityNUM);
		break;
	case 10: // 3 pt, 1 to 2
		xoff[0]=-1; xoff[1]=-1; // -1, -1
		toff[0]=0;  toff[1]=1;  // 0, 1
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 11: // 3 pt, 2 to 1
		xoff[0]=-1; xoff[1]=0;
		toff[0]=1;  toff[1]=1;
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 12: // 3 pt, zig zag
		xoff[0]=-1; xoff[1]=0; 
		toff[0]=1;  toff[1]=2;  
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 13: // 3 pt, fwd, diag
		xoff[0]=-1; xoff[1]=-1; // -1, -1
		toff[0]=1;  toff[1]=2;  // 0, 1
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 14: // 3 pt, diag, fwd
		xoff[0]=0; xoff[1]=-1; // -1, -1
		toff[0]=1;  toff[1]=2;  // 0, 1
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 15: // 4 pt, vertical tetris piece
		xoff[0]=0; xoff[1]=-1; xoff[2]=-1;// 0, -1, -1
		toff[0]=1;  toff[1]=1; toff[2]=2;// 1, 1, 2
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 16: // horizontal knight 2-1
		xoff[0]=-1; xoff[1]=-2; // -1, -1
		toff[0]=1;  toff[1]=1;  // 0, 1
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 17: // horizontal knight 1-2
		xoff[0]=-1; xoff[1]=-2; // -1, -1
		toff[0]=0;  toff[1]=1;  // 0, 1
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 18: // 4 pt, horizontal tetris piece
		xoff[0]=-1; xoff[1]=-1; xoff[2]=-2;// 0, -1, -1
		toff[0]=0;  toff[1]=1; toff[2]=1;// 1, 1, 2
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;
	case 19: // 2 pt, variable toff
		xoff[0]=-1; // nearest neighbor
		toff[0]=select_frames;  // some number back...
		updateparitypixvals_arbitrarypattern(parity,parityNUM,toff,xoff);
		break;


};
};

void EnforcedCorr1D::randomizepixels()
{
	bool getswhite;
	for (int ii=0;ii<numberofpixels;ii++)
	{
		getswhite=(((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
		if (getswhite)
			pixvals[ii]=1;
		else
			pixvals[ii]=0;
	};
	updateparityvals();
};

void EnforcedCorr1D::updateparitypixvals(int P,int NUM)
{
	bool getswhite;
	int currpar;
	int startpos, endpos;
	
	// deal with edges
	if (circularBC)
	{
		startpos=0;
		endpos=numberofpixels;
	}
	else
	{
		if (direction == 1)
		{
			startpos=NUM-1;
			endpos=numberofpixels;

			for (int ii=0;ii<startpos;ii++)
			{
				getswhite=(((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
				if (getswhite)
					pixvals[ii]=1;
				else
					pixvals[ii]=0;
			};
		}
		else
		{
			startpos=0;
			endpos=numberofpixels-NUM+1;

			for (int ii=endpos;ii<numberofpixels;ii++)
			{
				getswhite=(((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
				if (getswhite)
					pixvals[ii]=1;
				else
					pixvals[ii]=0;
			};
		};
	};

	// deal with core values...
	for (int ii=startpos;ii<endpos;ii++)
	{
		currpar=computeparity(NUM,ii);

		if (currpar == P)
			pixvals[ii]=1;
		else
			pixvals[ii]=0;

	};

	updateparityvals();
};

void EnforcedCorr1D::updateparitypixvals_arbitrarypattern(int P, int NUM, int (&toff)[10], int (&xoff)[10])
{
	// NOTE ON THIS FUNCTION -- FOR CONTEMPORARY TIME, XOFF MUST ALWAYS BE NEGATIVE!
	// also problem: must initialize first pixvals in loop randomly if toff has a 0... to be changed!
	bool getswhite;
	int currpar;
	int startpos, endpos;
	
	// INITIALIZE TO RANDOM VALUES; when loop starts if toff has a 0, first is random, then uses neighboring values. all good?
	for (int ii=0;ii<numberofpixels;ii++)
	{
		getswhite=(((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
		if (getswhite)
			pixvals[ii]=1;
		else
			pixvals[ii]=0;
	};

	// deal with edges
	if (circularBC)
	{
		startpos=0;
		endpos=numberofpixels;
		if (direction == 1)
		{
			for (int ii=startpos;ii<endpos;ii++)
			{
				currpar=computeparity_arbitrarypattern(NUM,ii,toff,xoff);

				if (currpar == P)
					pixvals[ii]=1;
				else
					pixvals[ii]=0;

			}
		}
		else // go other direction
		{
			for (int ii=endpos-1;ii>startpos-1;ii--)
			{
				currpar=computeparity_arbitrarypattern(NUM,ii,toff,xoff);

				if (currpar == P)
					pixvals[ii]=1;
				else
					pixvals[ii]=0;

			};
		};

	}
	else // this code so far only works with periodic BC... if not, make all white to be clear to check it out!
	{
		for (int ii=0;ii<numberofpixels;ii++)
		{
			pixvals[ii]=1;
		};
	};

	// deal with core values...

	if (updatetype==10)
	{
		;;
	};
	
	updateparityvals();
};


int EnforcedCorr1D::computeparity(int NUM, int loc)
// assumes circularBC; update must take care of this...
{
	int parityout=1;

	for (int jj=0; jj<NUM-1; jj++)
	{
		int currloc;
		if (direction==1) // something wrong here?
			currloc = (loc-1-jj);
		else
			currloc = (loc+1+jj) % numberofpixels;

		if (currloc<0)
			currloc = numberofpixels + currloc;

		parityout *= parityvals [jj] [currloc];
	};
	return parityout;	
};


int EnforcedCorr1D::computeparity_arbitrarypattern(int NUM, int loc,int (&toff)[10], int (&xoff)[10])
// assumes circularBC; update must take care of this...
{
	int parityout=1;
	int currloc;

	//static int looksee[10];

	for (int jj=0; jj<NUM-1; jj++)
	{
		if (direction==1) // compute current location...
		{
			currloc = (loc+xoff[jj]);
		}
		else
		{
			currloc = (loc-xoff[jj]);
		};

		// special cases
		if (currloc<0)
		{
			currloc = numberofpixels + currloc;
		};
		if (currloc>=numberofpixels)
		{
			currloc = currloc % numberofpixels; // doing mod carefully in here...
		};
		
		// use pix vals rather than parity vals...
		if (toff[jj] == 0) // current time frame; must be moving in correct direction... appears correct.
		{
			parityout *= ( (pixvals [currloc] * 2) - 1);
		}
		else
		{
			parityout *= parityvals [toff[jj]-1] [currloc];
		};

		/*
		if (toff[jj] == 0) // current time frame; must be moving in correct direction... appears correct.
			looksee[jj] = (pixvals [currloc] * 2 - 1);
		else
			looksee[jj] = parityvals [toff[jj]-1] [currloc];
		*/

	};
	return parityout;	
};


void EnforcedCorr1D::update1(int x1t1,int x2t2)
{
	bool getswhite;

	if (direction == 1)
	for (int ii=0;ii<numberofpixels;ii++)
	{
		if ((ii==0) && (!circularBC))
		{
			getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
			if (getswhite)
				pixvalsnew[ii]=1;
			else
				pixvalsnew[ii]=0;
		}
		else
		{
			if (ii>0)
				if (pixvals[ii-1]==x1t1)
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
			else
				if (pixvals[numberofpixels-1]==x1t1)
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
		};
		
	}
	else
	for (int ii=0;ii<numberofpixels;ii++)
	{
		if ((ii==numberofpixels-1) && (!circularBC))
		{
			getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
			if (getswhite)
				pixvalsnew[ii]=1;
			else
				pixvalsnew[ii]=0;
		}
		else
		{
			if (ii<numberofpixels-1)
				if (pixvals[ii+1]==x1t1)
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
			else
				if (pixvals[0]==x1t1)
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
		};
		
	}

	copypixvalsnewtopixvals();
	updateparityvals();

};

void EnforcedCorr1D::update1corr(int x1t1,int x2t2, float c)
{
	bool getswhite;
	bool updateit;

	if (direction == 1)
	for (int ii=0;ii<numberofpixels;ii++)
	{
		updateit = (((float)rand())/((float)RAND_MAX+1.0f) < c);
		getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
		if ((ii==0) && (!circularBC))
		{
			if (getswhite)
				pixvalsnew[ii]=1;
			else
				pixvalsnew[ii]=0;
		}
		else
		{
			if (ii>0)
				if ((pixvals[ii-1]==x1t1) && (updateit))
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
			else
				if ((pixvals[numberofpixels-1]==x1t1) && (updateit))
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
		};
		
	}
	else
	for (int ii=0;ii<numberofpixels;ii++)
	{
		updateit = (((float)rand())/((float)RAND_MAX+1.0f) < c);
		getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
		if ((ii==numberofpixels-1) && (!circularBC))
		{
			if (getswhite)
				pixvalsnew[ii]=1;
			else
				pixvalsnew[ii]=0;
		}
		else
		{
			if (ii<numberofpixels-1)
				if ((pixvals[ii+1]==x1t1) && (updateit))
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
			else
				if ((pixvals[0]==x1t1) && (updateit))
				{
					pixvalsnew[ii]=x2t2;
				}
				else
				{
					if (getswhite)
						pixvalsnew[ii]=1;
					else
						pixvalsnew[ii]=0;
				}
		};
		
	}

	copypixvalsnewtopixvals();
	updateparityvals();

};






void EnforcedCorr1D::update2(int parity)
{
	static int xtrip1[4][3]={{1,1,1},{1,0,0},{0,1,0},{0,0,1}}; // 2 different parities
	static int xtrip2[4][3]={{0,0,0},{0,1,1},{1,0,1},{1,1,0}};
	int (*xuse)[4][3];
	bool getswhite;
	bool changed;

	if (parity == 1)
		xuse=&xtrip1; // POINT to it
	else
		xuse=&xtrip2; // point to it

	//meat of loop
	for (int ii=0; ii<numberofpixels;ii++)
	{
		changed=FALSE;
		for (int jj=0;jj<4;jj++)
			if (test2((*xuse)[jj][0],(*xuse)[jj][1],ii))
			{
				pixvalsnew[ii]=(*xuse)[jj][2];
				changed=TRUE;
			};

		if (!changed)
		{
			getswhite = (((float)rand())/((float)RAND_MAX+1.0f) < Pwhite);
			if (getswhite)
				pixvalsnew[ii]=1;
			else
				pixvalsnew[ii]=0;
		}

	};
	copypixvalsnewtopixvals();
	updateparityvals();

};

bool EnforcedCorr1D::test2(int v1,int v2,int position)
{
	if (position>0)
		return ((pixvals[position-1]==v1) && (pixvals[position]==v2));
	else
	{ // first position...
		if (!circularBC)
			return (FALSE);
		else
			return ((pixvals[numberofpixels-1]==v1) && (pixvals[position]==v2));
	};
	
};

void EnforcedCorr1D::copypixvalstopixvalsnew()
{
	for (int ii=0; ii<numberofpixels;ii++)
		pixvalsnew[ii]=pixvals[ii];
};

void EnforcedCorr1D::copypixvalsnewtopixvals()
{
	for (int ii=0; ii<numberofpixels;ii++)
		pixvals[ii]=pixvalsnew[ii];
};

int EnforcedCorr1D::getpixelvalue(int pixnum)
{
	return pixvals[pixnum];
};

void EnforcedCorr1D::setDir(int dir)
{
	direction = dir;
};

void EnforcedCorr1D::setBackFrames(int inframes)
{
	select_frames = inframes;
};














