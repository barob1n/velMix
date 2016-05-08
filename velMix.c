#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "defns.h"
#include "segy.h"
#include "segyIO_class.h"
#include <time.h>
#include <omp.h>

#define FLUSH while (getchar() !='\n');
void velMix(float ***velTr, float ***velOutTr,int ns, int numTr,int iNdx ,int window,int maxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, float thresh, int rate);
int slopeChk( float slopeOld, float slopeIXl, float dipMax);
int nearest (float input);
int boundCheck (int sampPlusWin, int slope, int ns);

/***********************************************************************
 * Function: main
 * Input/output:  input seismic file name, output dip file name, 
 * 		output dip file name, number of traces used in semblance calc,
 * 		dip as max ms/trace allowable, semblance window size
 * 	      
 * Descrip: Mixes velocity voulume along geologic dips.
 * ********************************************************************/
int main(int argc, char *argv[])
{
omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(7);

FILE *finVel=NULL, *finX=NULL, *finI=NULL, *fsemb=NULL, *foutVel=NULL;
int   fmt, i, j , itr, ierr, k;
int   ns;	/* Number of samples in a trace	*/
float dts;	/* Sample rate in seconds	*/
segy  thdr;	/* Segy trace header		*/
segy thdrI, thdrX, thdrS;
bhed bhI,bhX,bhS;
bhed  bh;	/* Segy file binary header	*/
char chdr[3202];	/* segy file character header */
int   endian=1; /* Which side of the egg	*/
int   status=0; /* Return status for segy calls, I should check*/
char  txthdr[3201];	/* The character file header*/

int choice;
char buff[80];

segy ***trGthHdr;
float *tempTrI; /*holding place for traces as they are read in*/
float *tempTrX;
float *tempTrS;
float *tempTrV;
float ***dipTrI;
float ***dipTrX;
float ***semTr;
float ***velTr;
float ***velOutTr;

int *xNdxMx; /*Keeps track of number of xline on each inline*/
int *xNdxMn; /*first xline on each inline*/
int statusR; /*Status of last segyRead*/
int statusW; /*Status of last segyWrite*/
int numTr; /*Radius of traces, about the trace being examined, used in
             the semblance-dip calculation*/
int maxDip; /*maximum allowed dip, IN SAMPLES PER TRACE*/
int window; /*Size of sembalnce window in number of samps*/
int shift; /*shifts starting sample point down for dip calc.  if you
				started at the 0th sample, the dip calc would run into
				* negative sample points*/
int iOrigin = 0; /*origin of inline grid*/
int xOrigin = 0; /*origin of xline grid*/
int iMax = 0;
int xMax = 0;
int line=0;
int iline = 0;
int xline = 0;
int ilineOld;			
int ilineCnt=0;
float thresh = 0;
int skip = 0;
int rate;
int ii;
int num10Percent;
int nsX,nsI,nsS;

time_t now;

float temp;
endian = checkEndian();

thdr.iline = thdr.ep;
thdr.xline=thdr.cdpt;

	if(argc ==14){
/*Open fileS as readable/writable binary*/
		finVel = fopen(argv[1], "rb");
		if (finVel == NULL) {
			printf("Unable to open the Input file.  Please check the name.\n");
		return -1;
		}
		finI = fopen(argv[2], "rb");
		if (finI == NULL) {
			printf("Unable to open the iline dip file.  Please check the name.\n");
		return -1;
		}
		finX = fopen(argv[3], "rb");
		if (finX == NULL) {
			printf("Unable to open the xline dip file.  Please check the name.\n");
		return -1;
		}
		fsemb = fopen(argv[4], "rb");
		if (fsemb == NULL) {
			printf("Unable to open the semblance input file.  Please check the name.\n");
		return -1;
		}
			
		foutVel=fopen(argv[5],"r");
		if(foutVel!=NULL){
			printf("\noutput velocity file Seems to exist. ");
			printf("\nWould you like to overwrite it (-1 = yes)?");
			fgets(buff,80,stdin);
			sscanf(buff,"%d",&choice);
			if (choice == -1){
				fclose(foutVel);
				foutVel=fopen(argv[5],"wb");
			}else{
				fclose(foutVel);
				return 0;
			}
		choice=0;
		}else {
			foutVel=fopen(argv[5],"wb");
		}

		if(sscanf(argv[6],"%d", &numTr)!=1){
			printf("\nnumTr entered not valid. Exiting... ");
			return -1;
		}
		if(sscanf(argv[7],"%d", &maxDip)!=1){
			printf("\nmaxDip entered not valid. Exiting... ");
			return -1;
		}
		if(sscanf(argv[8],"%d", &window) !=1){
			printf("\nwindow entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[9],"%f", &thresh) !=1){
			printf("\nwindow entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[10], "%d", &iOrigin) != 1){
			printf("\niOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[11], "%d", &xOrigin) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[12], "%d", &iMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}
		if(sscanf(argv[13], "%d", &xMax) != 1){
			printf("\nxOrigin entered not valid. Exiting...");
			return -1;
		}		
	} else {
	fprintf(stderr,"\n***************************************************************************\n");
      fprintf(stderr,"Mixes velocities along geologic dip. Currently, the input dip files\n\n");
      fprintf(stderr,"are expected to be in ms/trace. The program rounds dip to nearest sample. \n\n");
      fprintf(stderr,"Program expects the following command line: \n ");
      fprintf(stderr,"sembDip <vel.sgy><dipI.sgy><dipX.sgy><sem.sgy><Out.sgy><Num><maxDip><win>\n");
      fprintf(stderr," <thresh><iOrigin><xOrigin><iMax><xMax>\n\n");
      fprintf(stderr,"vel.sgy: Input velocity file. \n\n");
      fprintf(stderr,"dipI.sgy: Output Inline dip sgy filename.\n\n");
      fprintf(stderr,"dipX.sgy: Output Xline dip sgy filename.\n\n");
      fprintf(stderr,"sem.sgy: Output semblance  sgy filename.\n\n");
      fprintf(stderr,"Out.sgy: Output velocity  sgy filename.\n\n");
      fprintf(stderr,"Num: Number of traces for semblance given as  \n");
      fprintf(stderr,"     distance from trace being analysed. \n");
      fprintf(stderr,"     ex. numTr = 1 uses 9 traces. 3 by 3 block.\n");
      fprintf(stderr,"                 o--o--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--x--o\n");
      fprintf(stderr,"                 |  |  |\n");
      fprintf(stderr,"                 o--o--o\n\n");
      fprintf(stderr,"maxDip: Maximum value for difference for conflicting dips.\n\n");
      fprintf(stderr,"win: Number of samps used for semblance analysis.\n\n");
      fprintf(stderr,"thresh: minimum value for semblance - velocity smearing stops if semblance\n");
      fprintf(stderr,"        is less than this value.\n\n");
      fprintf(stderr,"iOrigin: Starting inline number.\n\n");
      fprintf(stderr,"xOrigin: Starting xnline number.\n\n");
      fprintf(stderr,"iMax: Max iline. \n\n");
      fprintf(stderr,"xMax: Max xline.");
    fprintf(stderr, "\n**************************************************************************\n");
      return 0;
    } 

	/*print parameters. ask user to verify they are correct-option to abort*/
	printf("\nInput paramters are:");
	printf("\n\n  Input Velocity: %s", argv[1]);
	printf("\n  Input dipI: %s", argv[2]);
	printf("\n  Input dipX: %s", argv[3]);
	printf("\n  Input semb: %s", argv[4]);
	printf("\n  Output velocity:  %s",argv[5] );
	printf("\n  Trace radius: %d", numTr);
	printf("\n  Max dip (in samples!): %d", maxDip);
	printf("\n  Window (in samples!): %d", window);
	printf("\n  threshhold: %f", thresh);
	printf("\n  inline origin: %d", iOrigin);
	printf("\n  xline origin: %d", xOrigin);
	printf("\n  iMax : %d", iMax);
	printf("\n  xMax: %d", xMax);
	
	printf("\n\nAre these correct? (-1: exit)");
	fgets(buff,80,stdin);
	sscanf(buff,"%d",&choice);
	if(choice == -1){ return 0;} 
	
	segyReadHeader(finX, chdr, &bh, endian);
	nsX=bh.hns;
	
	segyReadHeader(finI, chdr, &bh, endian);
	nsI=bh.hns;
	
	segyReadHeader(fsemb, chdr, &bh, endian);
	nsS=bh.hns;
	
	segyReadHeader(finVel, chdr, &bh, endian);
	ns=bh.hns;
	rate = bh.hdt;
	rate=rate/1000;
	
	if(nsX !=ns){
		printf("\n\n Xline dip file's number of samps not matching input velocity file's number of samps!");
		return -1;
	}else if(nsI != ns){
		printf("\n\n Inline dip file's number of samps not matching input velocity file's number of samps!");
		return -1;
	}else if(nsS != ns){
		printf("\n\n Semblance file's number of samps not matching input velocity file's number of samps!");
		return -1;
	}
	
	printf("\nNumber of samples: %d", ns);
	
	if(segyWriteHeader(foutVel, &chdr, &bh, endian)!=0){
		printf("\n Unable to write header for output data!! \n");
		return -1;
	}
	
	/*Converting numTr to the width of cube ex. 3x3*/
	numTr=(numTr)*2+1;
	
	/*I prefer the window to be an odd number so as to be symmetric*/
	if(window % 2 ==0){
		window = window + 1; 
	}
	/*allocating memory for the traces used for semblance calc/analysis*/
	velTr = (float***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         velTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 velTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	dipTrI = (float ***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         dipTrI[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 dipTrI[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 dipTrX = (float ***)calloc((numTr) , sizeof(float **));
     for (i=0; i<(numTr); ++i){
         dipTrX[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 dipTrX[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }	
	 semTr = (float ***)calloc((numTr) , sizeof(float **));
     for (i=0; i<(numTr); ++i){
         semTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 semTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	velOutTr = (float***)calloc((numTr) , sizeof(float **));
    for (i=0; i<(numTr); ++i){
         velOutTr[i] = (float **)calloc((xMax-xOrigin+2),sizeof(float*));
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 velOutTr[i][j] = (float *) calloc((ns+1),sizeof(float));
		 }
	 }
	 	
	 trGthHdr = (int ***)calloc((numTr) , sizeof(int **));
	 for (i=0; i<(numTr); ++i){
         trGthHdr[i] = (int **)calloc((xMax-xOrigin+2),sizeof(int*)); 
         for ( j =0; j < (xMax-xOrigin+2); ++j){
			 trGthHdr[i][j] = (int *) calloc((1),sizeof(segy));
		 }
	 }

	tempTrI = (float *)calloc((ns+1),sizeof(float));
	tempTrX = (float *)calloc((ns+1),sizeof(float));
	tempTrS = (float *)calloc((ns+1),sizeof(float)); 
	tempTrV = (float *)calloc((ns+1),sizeof(float));
	 
	xNdxMx = (int *)calloc ((iMax - iOrigin +2 + numTr), sizeof(int)); //to track the starting and ending xlines
	xNdxMn = (int *)calloc ((iMax - iOrigin +2 + numTr), sizeof(int)); //to track the starting and ending xlines

	 /*Reading in frist Trace, I am currently reading in the vel trace last since it is the header I am most concerned with*/
	 /* At somepoint a check should be implemented to be certain the the starting xbin and ybin are the same of all files*/
	
	statusR = segyReadTrace(finI, &bhI, &thdrI,tempTrI, ns, endian);
	statusR = segyReadTrace(finX, &bhX, &thdrX,tempTrX, ns, endian);
	statusR = segyReadTrace(fsemb, &bhS, &thdrS,tempTrS, ns, endian);
	statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
	
	/*im setting the iOrigin to the first inline read in-will be useful later*/
	if(thdr.iline != iOrigin){
		printf("\n\nFirst inline read does not match iOrigin: %d\n",iOrigin);
		printf("Setting iOrigin equal to first iline read in!: %d\n\n",thdr.iline);
		iOrigin=thdr.iline;
	}
	/*are the values sensabale?*/
	if (thdr.iline > iMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}else if (thdr.xline - xOrigin < 0){
		printf("\n\nFirst xline read less than xOrigin.\n\n");
		return 0;
	}else if (thdr.xline > xMax){
		printf("\n\nFirst iline read great than iMax.\n\n");
		return 0;
	}

	/*save the xline values and set iline switch -really only need to set the Mn value*/
	xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;//the first xline will be zero if it is at the xline axis
	ilineOld=thdr.iline;
	
	/*read in the first inline  of traces*/
	while ( thdr.iline - iOrigin < numTr){
		
		memcpy(dipTrI[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrI,(ns+1)*sizeof(float));
		memcpy(dipTrX[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrX,(ns+1)*sizeof(float));		
		memcpy(semTr[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrS,(ns+1)*sizeof(float));
		
		memcpy(velTr[thdr.iline - iOrigin][thdr.xline -xOrigin], tempTrV,(ns+1)*sizeof(float));
		memcpy(trGthHdr[thdr.iline - iOrigin][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
		
		//memcpy(dipTrI[thdr.iline - iOrigin + numTr/2][thdr.xline -xOrigin], tempTrI,(ns+1)*sizeof(float));
		//memcpy(dipTrX[thdr.iline - iOrigin + numTr/2][thdr.xline -xOrigin], tempTrX,(ns+1)*sizeof(float));		
		//memcpy(semTr[thdr.iline - iOrigin + numTr/2][thdr.xline -xOrigin], tempTrS,(ns+1)*sizeof(float));
		
		//memcpy(velTr[thdr.iline - iOrigin + numTr/2][thdr.xline -xOrigin], tempTrV,(ns+1)*sizeof(float));
		//memcpy(trGthHdr[thdr.iline - iOrigin + numTr/2][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
		
		statusR = segyReadTrace(finI, &bhI, &thdrI,tempTrI, ns, endian);
		statusR = segyReadTrace(finX, &bhX, &thdrX,tempTrX, ns, endian);
		statusR = segyReadTrace(fsemb, &bhS, &thdrS,tempTrS, ns, endian);
		statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
		
		if (thdr.iline > iMax){
			printf("\n\nFound iline greater than iMax.\n\n");
			printf("\n\n iline: ", thdr.iline);
			return 0;
		}else if (thdr.iline < iOrigin){
			printf("\n\nFound iline less than iOrigin.\n\n");
			return 0;
		}else if (thdr.xline - xOrigin < 0){
			printf("\n\nFound xline  less than xOrigin.\n\n");
			return 0;
		}else if (thdr.xline > xMax){
			printf("\n\nFound xline greater than xMax.\n\n");
			return 0;
		}	
			
		/*if it gets to end of reading before having read in enough traces abort
		*else end o file and you have enough lines exit while loop */
		if(statusR!=0){
			if(thdr.iline - iOrigin < numTr -1 ){
				printf("\nNot enough ilines for a single execution calc!\n");
				return  0;
			}else{
					
				break; /*minimum number of lines met...exit loop*/
			}
		}
		    /*the xline number just read is large than xNdxMx*/
			if ( (thdr.xline-xOrigin) > xNdxMx[thdr.iline-iOrigin]) xNdxMx[thdr.iline - iOrigin]=thdr.xline - xOrigin;
			if (ilineOld != thdr.iline){
				xNdxMn[thdr.iline-iOrigin]=thdr.xline-xOrigin;
				//printf("\n I found %d %d", thdr.iline);
				ilineOld=thdr.iline;
			}
	}

	/*copy the value read last, which caused while loop to exit, by the way 
	 * I could just compute ilineCnt from the difference between thdr.iline and iOrigin
	 * since this is a bit confusing*/
	ilineCnt=numTr/2;
	
	
	printf("\n\n ****************");
			time(&now);
			printf("\n Start Time: %s Percent Done: %d",  ctime(&now),0 );
	
	/*compute semblance and dips*/
	velMix(velTr, velOutTr, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,rate);
	
	/*outputing status*/
	
	num10Percent = (iMax-iOrigin + 1)/10;
	if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
		printf("\n\n ****************");
		time(&now);
		printf("\n Finished Inline: %d Time: %s Percent Done: ", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline - iOrigin + 1)/num10Percent*10 );
	}

	/*spits out first numTr/2 +1  lines*/
	for (i=0;i<=numTr/2;++i){
		/*output first numTr traces for begining of inline*/
		for(j = xNdxMn[i]; j < xNdxMn[i] + numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMn[i] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[i]+numTr/2; j <= xNdxMx[i] - numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[i] - numTr/2 + 1; j <= xNdxMx[i]; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMx[i] - numTr/2], ns, endian);
		}
	}
	
   while (statusR == 0){ 
   //for (ii=0;ii< iMax-iOrigin + numTr/2;++ii){
	   
		/*Shift all the traces back by one, meaning that
		* trGth[1] -> trGth[0], trGth[2] -> [1], and tempTr -> trGth[numTr-1]*/
		for (k =0 ;k<numTr-1;++k){
			for(j=0;j<xMax-xOrigin+1;++j){

				memcpy(dipTrI[k][j],dipTrI[k+1][j],(ns+1)*sizeof(float));
				memcpy(dipTrX[k][j],dipTrX[k+1][j],(ns+1)*sizeof(float));
				memcpy(semTr[k][j],semTr[k+1][j],(ns+1)*sizeof(float));
				memcpy(velTr[k][j],velTr[k+1][j],(ns+1)*sizeof(float));
				memcpy(trGthHdr[k][j],trGthHdr[k+1][j], sizeof(*trGthHdr[0][0]));
			}
		}

		/*setting last inline to zero before I read a new inline into it*/
		for (j =0; j < xMax-xOrigin +1;++j){
			memset(velTr[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(trGthHdr[numTr-1][j],0,sizeof(*trGthHdr[0][0]));
			memset(dipTrI[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(dipTrX[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(semTr[numTr-1][j],0,(ns+1)*sizeof(float));
			memset(velOutTr[numTr/2][j],0,(ns+1)*sizeof(float));
		}
	
		//while(statusR ==0 && thdr.iline ==ilineOld){
		while (statusR==0 && thdr.iline ==ilineOld){	
			
			memcpy(dipTrI[numTr - 1][thdr.xline - xOrigin], tempTrI, (ns + 1)*sizeof(float));
			memcpy(dipTrX[numTr - 1][thdr.xline - xOrigin], tempTrX, (ns + 1)*sizeof(float));
			memcpy(semTr[numTr - 1][thdr.xline - xOrigin], tempTrS, (ns + 1)*sizeof(float));
			memcpy(velTr[numTr - 1][thdr.xline - xOrigin], tempTrV, (ns + 1)*sizeof(float));
			memcpy(trGthHdr[numTr - 1][thdr.xline -xOrigin], &thdr, sizeof(*trGthHdr[0][0]));
			
			if ((thdr.xline - xOrigin) > xNdxMx[thdr.iline - iOrigin]) xNdxMx[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			
			statusR = segyReadTrace(finI, &bhI, &thdrI,tempTrI, ns, endian);
			statusR = segyReadTrace(finX, &bhX, &thdrX,tempTrX, ns, endian);
			statusR = segyReadTrace(fsemb, &bhS, &thdrS,tempTrS, ns, endian);
			statusR = segyReadTrace(finVel, &bh, &thdr,tempTrV, ns, endian);
			
			if(statusR==0){
				if (thdr.iline > iMax){
					printf("\n\nFound iline greater than iMax.\n\n");
					return 0;
				}else if (thdr.iline < iOrigin){
					printf("\n\nFound iline less than iOrigin.\n\n");
					printf("\n\n iline: %d %d %d %d %d %d", thdr.iline, thdr.xline, status, xNdxMn[ilineCnt], xNdxMx[ilineCnt],ilineCnt);
					return 0;
				}else if (thdr.xline - xOrigin < 0){
					printf("\n\nFound xline  less than xOrigin.\n\n");
					return 0;
				}else if (thdr.xline > xMax){
					printf("\n\nFound xline greater than xMax.\n\n");
					return 0; 
				}	
			}		
		}

		if (statusR == 0){
			xNdxMn[thdr.iline - iOrigin] = thdr.xline - xOrigin;
			
			ilineOld = thdr.iline;
		}
			
		ilineCnt++;
		
		velMix(velTr, velOutTr, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,rate);
		
		if((((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline) - iOrigin + 1)%num10Percent==0){
			printf("\n\n ****************");
			time(&now);
			printf("\n Finished Inline: %d Time: %s Percent Done: %d", (*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline, ctime(&now),((*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline-iOrigin + 1)/num10Percent*10 );
			//printf("\\ Percent Done: 
			//printf("\n Time: %s", ctime(&now));
		}
		
	
		/*output current inline*/
		i=numTr/2;
		for(j = xNdxMn[ilineCnt]; j < xNdxMn[ilineCnt] + numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMn[i] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[ilineCnt]+numTr/2; j <= xNdxMx[ilineCnt] - numTr/2; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[ilineCnt] - numTr/2 + 1; j <= xNdxMx[ilineCnt]; ++j){
			if((*trGthHdr[i][j]).iline == 0){
				(*trGthHdr[i][j]).iline=(*trGthHdr[i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[i][j], velOutTr[numTr/2][xNdxMx[i] - numTr/2], ns, endian);
		}
	}
	
	for (i=1;i<=numTr/2;++i){
		/*output last numTr traces for begining of inline*/
		for(j = xNdxMn[ilineCnt + i]; j < xNdxMn[ilineCnt + i] + numTr/2; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				//(*trGthHdr[ilineCnt + i][j]).iline=(*trGthHdr[numTr/2][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[numTr/2 + i][j]).iline=ilineCnt + i;
				(*trGthHdr[numTr/2 + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][xNdxMn[ilineCnt + i ] + numTr/2], ns, endian);
		}
		
		//for (j = xNdxMn[i]; j <= xNdxMx[i]; ++j){
		for (j= xNdxMn[ilineCnt + i]+numTr/2; j <= xNdxMx[ilineCnt + i] - numTr/2; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				//(*trGthHdr[ilineCnt + i][j]).iline=(*trGthHdr[numTr/2 + i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[ilineCnt + i][j]).iline=ilineCnt + i;
				(*trGthHdr[ilineCnt + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][j], ns, endian);
		}
	
		for(j = xNdxMx[i] - numTr/2 + 1; j <= xNdxMx[i]; ++j){
			if((*trGthHdr[numTr/2 + i][j]).iline == 0){
				//(*trGthHdr[ilineCnt + i][j]).iline=(*trGthHdr[numTr/2 + i][xNdxMn[ilineCnt]]).iline;
				(*trGthHdr[ilineCnt + i][j]).iline=ilineCnt + i;
				(*trGthHdr[ilineCnt + i][j]).xline=xOrigin + j;
			}
			segyWriteTrace(foutVel, trGthHdr[numTr/2 + i][j], velOutTr[numTr/2][xNdxMx[ilineCnt + i] - numTr/2], ns, endian);
		}
	}

printf("\n\nDONE!\n\n");

free(semTr);
free(velTr);
free (velOutTr);
free(trGthHdr);
free(dipTrI);
free(dipTrX);
free(tempTrI);
free(tempTrX);
free(tempTrS);
free(tempTrV);
free(xNdxMx);
free(xNdxMn);
return 0;
}
/*
************************************************************************
*/
void doMessage(char *str)
{
   fprintf(stderr, "%s\n", str);
}
/*
***********************************************************************
*/
/***********************************************************************
 * Function: Velocity Mixer
 * Input: Inline dip volume, Xline dip volume, semblance value volume,
 * 		  velocity volume, numTr, window, maxDip.s
 * 
 * Description:  Mixes velocities witing a radius of numTr along dips.
***********************************************************************/
void velMix(float ***velTr, float ***velOutTr,int ns, int numTr,int iNdx ,int window,int maxDip, float ***dipTrI, float ***dipTrX, float***semTr, int *xNdxMn, int *xNdxMx, int iOrigin, int xOrigin, float thresh, int rate)
{
	int thread_id, nloops;	
#pragma omp parallel private(thread_id, nloops)
{	
	//velMix(velTr, velOutTr, ns, numTr,ilineCnt,window,maxDip,dipTrI,dipTrX,semTr,xNdxMn,xNdxMx,iOrigin,xOrigin,thresh,rate);
	int i,j,k,ii;
	int samp = 0; //sample value
	int xStart=xNdxMn[iNdx];//-xOrigin;//index number of first xline.  if xNdxMn = xOrigin, then xStart =0 
	int xEnd = xNdxMx[iNdx];//-xOrigin; //xEnd indx number of last xline
	int xNdx; //current array value of xline being processed on this particular inline
	int shift = window/2 ; //shift is the starting point for the samples.
	int halfTr = numTr / 2; //should have called this radius
	int halfWin = window / 2; //hafwindow of samples - should hve called winradi
	float velSum=0; /*Current sum of velocities*/
	int slope;
	int slopeDiag=0, slopeXl=0, slopeIl=0;
	float sumVel; //sum of valid velocities
	int numVel; //number of valid velocities 
	float semChk;
	int sampPlusWin; //Current sample including the window.
	float slopeOld;
	int ilDirec,xlDirec, ilMove,xlMove;

	#pragma omp for
	/*walking down the inline.  Remeber, this will start index numTr/2*/
	for (xNdx=xStart+(halfTr);xNdx <= xEnd-(halfTr); ++xNdx){
		/*computing the first sample point, then and window values, then i will loop through the rest*/
		for (samp=shift; samp <= ns-shift; samp++){		
		/*re-setting the number of vels and their sum for each samp of each trace*/
		numVel=0; 
		sumVel=0;	
			/*for the window get velocities - this  should defualt to one*/
			for (k=0; k<window; ++k){		
				/*setting the slopes to zero for each iteration of the window*/
				slopeXl=slopeIl=slopeDiag=0; 
				sampPlusWin=samp-halfWin+k;		
				/*getting velocity at center of gather for each window*/
				sumVel+=velTr[halfTr][xNdx][sampPlusWin];
				++numVel;

				/*Plus-Plus is the upper right quad, Plus-Minus is the lower right quad, 
				* Minus-Minus is lower left quad, and Plus-Minus is the upper left quad
				* Plus-Zero goes right along center, Minus-Zero goes left along center
				* Zero-Plus goes upward, and Zero-Minus goes downward*/
				for(ilDirec = -1*1; ilDirec<=1; ++ilDirec){
				for(xlDirec = -1*1;xlDirec<=1; ++xlDirec){

					/* if both ilDirec and xlDirec are zero then the algorithim will sit at origin summing over the velocies
					* so, to avoid this I am skipping over ther case where both ilDirec and xlDirec are zero.*/
					if(ilDirec != 0 || xlDirec !=0){
				
					/*make sure im not runninf off the trace*/	
					if( boundCheck ==1){
						slopeDiag = nearest((ilDirec*dipTrI[halfTr][xNdx][sampPlusWin] + xlDirec*dipTrX[halfTr ][xNdx ][sampPlusWin])/((float)rate));	
					}else {slopeDiag =0;}
				
					slopeOld=slopeDiag;

					i=1;
					ilMove=i*ilDirec; 
					xlMove=i*xlDirec;
					while (i<=halfTr && (sampPlusWin-slopeDiag) >=0 && (sampPlusWin-slopeDiag) < ns && semTr[halfTr +ilMove][xNdx + xlMove][sampPlusWin-slopeDiag] < thresh && slopeChk(slopeOld,slopeDiag,maxDip)==1){							
						/*if for some reason i get a zero velocity, then do not include it*/
						if( fabs(velTr[halfTr + ilMove][xNdx + xlMove][sampPlusWin-slopeDiag])>0){
							/*getting the velocity along the diagonal*/
							sumVel+=velTr[halfTr + ilMove][xNdx + xlMove][sampPlusWin-slopeDiag];
							++numVel; 
						}
				
						/*going up/down*/				
						j=i; /*starting a t diagonal [i,i] and moving up/down*/
						ilMove=i*ilDirec; 
						xlMove=j*xlDirec;
						slopeOld=slopeIl=slopeXl=slopeDiag; 
						slopeXl+= nearest((xlDirec*dipTrX[halfTr + ilMove][xNdx + xlMove][sampPlusWin-slopeDiag])/((float)rate));
					
						while (j<halfTr && (sampPlusWin-slopeXl) >=0 && (sampPlusWin-slopeXl) < ns && semTr[halfTr+ilMove][xNdx +xlMove][sampPlusWin-slopeXl] < thresh &&  slopeChk(slopeOld,slopeXl,maxDip)==1 ){
							if(fabs( velTr[halfTr + ilMove][xNdx + xlMove + xlDirec][sampPlusWin-slopeXl])>0){
								sumVel+=velTr[halfTr + ilMove][xNdx + xlMove + xlDirec][sampPlusWin-slopeXl];
								++numVel; 	
							}
							
							slopeOld=slopeXl;
							/*getting next slope*/
							slopeXl+= nearest((xlDirec*dipTrX[halfTr + ilMove][xNdx + xlMove + xlDirec][sampPlusWin-slopeXl])/((float)rate));
							++j;
							xlMove=j*xlDirec;
						}
					
						/*going right/left*/
						ilMove=i*ilDirec; /*if ilDirec is positive then the ilMove will increase. If ilDirec is neg then ilMove will decrease*/
						xlMove=i*xlDirec;
					
						slopeOld=slopeIl=slopeXl=slopeDiag;
						slopeIl+= nearest((ilDirec*dipTrI[halfTr + ilMove][xNdx +xlMove][sampPlusWin-slopeDiag])/((float)rate));
					
						j=i;
						ilMove=j*ilDirec; 
						xlMove=i*xlDirec;
						while (j<halfTr && (sampPlusWin-slopeIl) >=0 && (sampPlusWin-slopeIl) < ns && semTr[halfTr + ilMove][xNdx +xlMove][sampPlusWin-slopeIl] < thresh && slopeChk(slopeOld,slopeIl,maxDip)==1 ){	
						//while (j<halfTr && boundCheck(sampPlusWin,slopeIl,ns)==1 && semTr[halfTr + ilMove][xNdx +xlMove][sampPlusWin-slopeIl] < thresh && slopeChk(slopeOld,slopeIl,maxDip)==1 ){		
							if( fabs(velTr[halfTr + ilMove + ilDirec][xNdx + xlMove][sampPlusWin-slopeIl])>0){
								sumVel+=velTr[halfTr + ilMove + ilDirec][xNdx + xlMove][sampPlusWin-slopeIl];
								++numVel; 
							}
						
							slopeOld=slopeIl;
							slopeIl+= nearest((ilDirec*dipTrI[halfTr + ilMove + ilDirec][xNdx +xlMove][sampPlusWin-slopeIl])/((float)rate));
							++j;
							ilMove=j*ilDirec; 
						}	
	
						/*getting next diagonal and adding it to previous diag slope, moving down diagonal!!*/	
						slopeOld=slopeDiag;
						
						ilMove=i*ilDirec;
						ilMove=i*xlDirec;
						if( sampPlusWin-slopeDiag >=0 && sampPlusWin -slopeDiag < ns){
						//if(boundCheck(sampPlusWin,slopeDiag,ns)==1){
							slopeDiag +=nearest((ilDirec*dipTrI[halfTr +ilMove][xNdx +xlMove][sampPlusWin - slopeDiag] + xlDirec*dipTrX[halfTr + ilMove][xNdx+xlMove][sampPlusWin - slopeDiag])/((float)rate));
						}else {slopeDiag =0;}
	
						++i;
						ilMove=i*ilDirec; 
						xlMove=i*xlDirec;
				
					}
				}
			}
		}

		}
		if (numVel >0){
			velOutTr[halfTr][xNdx][samp]=sumVel/(numVel);
		}else {
			velOutTr[halfTr][xNdx][samp]=0;
		}
	}		
}
}
}
/***********************************************************************
 * Function: Checks slope for conflicting dips
 * Input: slopeOld - previously calculated slope, slopeIXl - most 
 * 	      recently computed slope, maxDip - largest acceptable value
 * 	      for the absolute value of the difference in the dips 
 * 		  
 * 
 * Description:  If the the dips are moving in different directions and
 * 	             are larger than specified value, return 0 else 1.
***********************************************************************/
int slopeChk( float slopeOld, float slopeIXl, float maxDip){
	if(slopeOld*slopeIXl < 0 &&  fabs(slopeOld -slopeIXl) > maxDip){
		return 0;
	}else {
		return 1;
	}
}

/***********************************************************************
 * Function: Nerest - finds nearest interger
 * Input: input- input floating point value
 * 		  
 * 
 * Description:  Adds +.5 if input is positve, subtracts .5 if input is
 * 	             negative.  Then performs integer truncation!
***********************************************************************/
int nearest (float input){
	if(input < 0){
		input = (int)(input - 0.5000);
	}else{
		input = (int)(input + 0.5000);
	}

return input;
}

int boundCheck (int sampPlusWin, int slope, int ns){
	if ((sampPlusWin - slope) >=0 && (sampPlusWin - slope) < ns){
		return 1;
	}else{
		return 0;
	}
}

	

