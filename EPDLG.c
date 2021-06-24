/*
EPDLG.c

Shawn Stone Jan 18 2018

This program calculates the position of the galileo space craft in (IAU_MOON) coordinates as well as the look direction of the EPD instrument in (IAU_MOON) for the chosen Jovian moon encounter.  
It reads in data from two sources: Right acention, declination, and twist (EME50) of the rotor, and the step motor and sector position.  The first set is from the magteams PDS data set and the second is from the EPDteams PDS data set.  The time tags are not the same in these data sets so an interpolation function was written for this code to find ra,dec,twist for the EPD at its time stamp.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SpiceUsr.h"// we use the NAIF spice system heavely

float interpolate(int i,int a,double time1); // declare interpolating function described below
void R19502SC(double a,double d,double t,double (*r)[3]);// returns a 3x3 rotation matrix from eme50 to space craft (SC) coordinates given ra,dec, and twist

// doing this with a structures is overkill.  I was just practicing.  These are basically read in records for mag and epd data
struct in_mag
{
  char UTC[80]; //the utc character string
  char sclk[80];//the space craft clocktime of UTC
  float b[3];//the magnetic field in sc coordinates
  float bmag;//magnitude of magnetic field
  double ra;//eme50 right acention of rotor in radians
  double dec;//eme50 declination of rotor in radians
  double twist;//eme50 twist angle in radians fro mag and not epd
  float spin;//eme50 rotor spin angle? (no idea)
} inmag;

struct in_EPD
{
  char UTC[80];//The UTC time stamp
  int step;//the stepper motor position
  int sector;//the sector of the EPD look direction
  float pitch;//The pitch angle of the detector based on mag at that position
  float phase;//the phase anngle of the detector based on mag at that position
  float rates[64]; //the count rates for the EPD at 16 sector resolution (64 rate channels)
} inepd;

double angles[10000][4];// the read in EME50 angles.  this is global to use with the interpolate function
int max=0;//the max records of angles
//The spice kernel set needed at G7 (1997-04-05T06:44:40.32 start)
#define SCLKKER "mk00062a.tsc"
#define LEAP "naif0008.tls"
#define PLTEPH_SPK "gll_951120_021126_raj2007.bsp"
#define FK "gll_v0.tf"
#define PCK "pk96030a.tpc"
 
  

int main()
{
  FILE *in1; //Pointer to input file
  FILE *out; //The output file for the encounter
 char filename[80];
 char UTC[80];//Used to find decimal time later

 double ra,dec,twist;//the eme50 angles 
 int i=0,j=0;  //dummies
 SpiceInt sc=-77; //NAIF code for the galileo spacecraft
 SpiceInt sat=503;//NAIF code for Ganymede
 SpiceChar MOON[20];
 double  motor[10000][3]; //array used to store the needed EPD data
 SpiceDouble et;//ephemeris time
 //J1950 is the FK4 inertial system
 SpiceDouble RSC2J1950 [3][3];//rotation matrix from SC->J1950
 SpiceDouble RJ19502SC [3][3];//rotation matrix from J1950->SC
 SpiceDouble RJ19502J2000[3][3];//rotation matrix from J1950->J200
 SpiceDouble RJ20002GSII[3][3];//rotation matrix from J2000->GSII (IAU_GANYMEDE)
 SpiceDouble RSC2J2000[3][3];//rotation matrix from SC->J2000 (not used in this code at this time but RSC2J000=RSC2J1950*RJ19502J2000)
 double Mposition[8]={315.0,0.0,30.6,61.2,91.8,122.4,153.,183.6};// motor postions 0-7 and their angles on the space craft(0 is behind the sheild)
 double sector[17]={-1,0.0,22.5,45.,67.5,90.,112.5,135.,157.5,180.,202.5,225.,247.5,270,292.5,315.,337.5};// sector postions 1-16 and their angles on the space craft
 double look[3];// the look direction of the EPD instrument (center of bore site)
 double deg2rad=3.14159/180.; //used in trig functions for degree to radian conversion
 double phi0=3.14159/16.;// offset angle of epd from the spin axis
 int lenout=80;//used in et2utc_c
 char subbuff [6];//used to get substrings from UTC to calc the float value of UTC
 double time;//contains the float time of UTC
 double state[3];//The position of galileo relative to ganymede in GSII(IAU_GANYMEDE) coordinates
 SpiceDouble lt;//the light time correction in spkez (not used)
 SpiceDouble e[3]={0.0,-.397881202,.917436945},s[3],n[3]; //vectors to calculate the twist angle for epd  the e vector is the EME50 normal
 char encounters[19][6]={"G2","C3","E4","G7","G8","C9","C10","E11","E12","E14","E15","E19","I24","E26","I27","G28","G29","C30","I31"};//the encounter I have mag and EPD for
 int NAIFCODE[19]={503,504,502,503,503,504,504,502,502,502,502,502,501,502,501,503,503,504,501}; //naif id's
 float counts[10000]; //rate file info for future simulations
 
 char magfiles[19][80]={"ORB02_GAN_IRC.TAB","ORB03_CALL_IRC.TAB","ORB04_EUR_IRC.TAB","ORB07_GAN_IRC.TAB","ORB08_GAN_IRC.TAB","ORB09_CALL_IRC.TAB","ORB10_CALL_IRC.TAB","ORB11_EUR_IRC.TAB","ORB12_EUR_IRC.TAB","ORB14_EUR_IRC.TAB","ORB15_EUR_IRC.TAB","ORB19_EUR_IRC.TAB","ORB24_IO_IRC.TAB","ORB26_EUR_IRC.TAB","ORB27_IO_IRC.TAB","ORB28_GAN_IRC.TAB","ORB29_GAN_IRC.TAB","ORB30_CALL_IRC.TAB","ORB31_IO_IRC.TAB",};//the mag files with ra and dec angles in eme50
 
 char epdfiles[19][80]={"C_ORB_02_1996_250_GAN_16.TAB","C_ORB_03_1996_309_CALL_16.TAB","C_ORB_04_1996_354_EUR_16.TAB","C_ORB_07_1997_095_GAN_16.TAB","C_ORB_08_1997_127_GAN_16.TAB","C_ORB_09_1997_176_CALL_16.TAB","C_ORB_10_1997_259_CALL_16.TAB","C_ORB_11_1997_310_EUR_16.TAB","C_ORB_12_1997_350_EUR_16.TAB","C_ORB_14_1998_088_EUR_16.TAB","C_ORB_15_1998_151_EUR_16.TAB","C_ORB_19_1999_032_EUR_16.TAB","C_ORB_24_1999_284_IO_16.TAB","C_ORB_26_2000_003_EUR_16.LBL","C_ORB_27_2000_053_IO_16.TAB","C_ORB_28_2000_141_GAN_16.TAB","C_ORB_29_2000_363_GAN_16.TAB","C_ORB_30_2001_145_CALL_16.TAB","C_ORB_31_2001_218_IO_16.TAB"};//the epd files with time and step and sector info
 SpiceBoolean found;
 char body[80]="IAU_";
 
 //Load in the needed spice kernels
furnsh_c(SCLKKER);
furnsh_c(LEAP);
furnsh_c( PLTEPH_SPK);
furnsh_c(FK);
furnsh_c(PCK);
//pick a moon encounter
 printf("Enter the number for the encounter\n");
   for(int k=0;k<19;k++)
     printf("%d: %s\n",k,encounters[k]);

 int enc;
 scanf("%d",&enc);
 sat=NAIFCODE[enc]; //set the naif id for the satellite
 bodc2n_c ( NAIFCODE[enc], 80, MOON,   &found );//get the ascii designation for the moon
 strcat(body,MOON); //turn it into the proper NAIF id in ascii
 
 strcpy(filename,magfiles[enc] );//get the mag file
 in1=fopen(filename,"r");// open mag file first

 while(!feof(in1))        //read it in
  {

    
    fscanf(in1,"%s %s %f %f %f %f %lf %lf %lf %f",inmag.UTC,inmag.sclk,&inmag.b[0],&inmag.b[1],&inmag.b[2],&inmag.bmag,&inmag.ra,&inmag.dec,&inmag.twist,&inmag.spin);

    str2et_c ( inmag.UTC, &et );//convert UTC time tag to J2000 ephemeris time
    angles[i][0]=et; //ephemeris time of data
    angles[i][1]=inmag.ra;//EME50 right acention radians
    angles[i][2]=inmag.dec;//EME50 declination radians
    angles[i][3]=inmag.twist;//EME50 twist anle radians (not used as twist for epd unless you want to de-spin
    i++;
  }
 max=i;// max number of records in angles
 fclose(in1); //close it!

 //At this point we have the eme50 data for the g2 encounter
 //Next is to load the 16 sector EPD data for step and sector

 strcpy(filename,epdfiles[enc]);//get the epd file
 in1=fopen(filename,"r");// open epd

while(!feof(in1))       
  {

    
    fscanf(in1,"%s %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",inepd.UTC,&inepd.step,&inepd.sector,&inepd.pitch,&inepd.phase,&inepd.rates[0],&inepd.rates[1],&inepd.rates[2],&inepd.rates[3],&inepd.rates[4],&inepd.rates[5],&inepd.rates[6],&inepd.rates[7],&inepd.rates[8],&inepd.rates[9],&inepd.rates[10],&inepd.rates[11],&inepd.rates[12],&inepd.rates[13],&inepd.rates[14],&inepd.rates[15],&inepd.rates[16],&inepd.rates[17],&inepd.rates[18],&inepd.rates[19],&inepd.rates[20],&inepd.rates[21],&inepd.rates[22],&inepd.rates[23],&inepd.rates[24],&inepd.rates[25],&inepd.rates[26],&inepd.rates[27],&inepd.rates[28],&inepd.rates[29],&inepd.rates[30],&inepd.rates[31],&inepd.rates[32],&inepd.rates[33],&inepd.rates[34],&inepd.rates[35],&inepd.rates[36],&inepd.rates[37],&inepd.rates[38],&inepd.rates[39],&inepd.rates[40],&inepd.rates[41],&inepd.rates[42],&inepd.rates[43],&inepd.rates[44],&inepd.rates[45],&inepd.rates[46],&inepd.rates[47],&inepd.rates[48],&inepd.rates[49],&inepd.rates[50],&inepd.rates[51],&inepd.rates[52],&inepd.rates[53],&inepd.rates[54],&inepd.rates[55],&inepd.rates[56],&inepd.rates[57],&inepd.rates[58],&inepd.rates[59],&inepd.rates[60],&inepd.rates[61],&inepd.rates[62],&inepd.rates[63]);

    str2et_c ( inepd.UTC, &et );
    motor[j][0]=et;//time to J2000
    motor[j][1]=inepd.step;//step motor position
    motor[j][2]=inepd.sector;//look sector
    counts[j]=inepd.rates[1];
    j++;
  }
 
 fclose(in1);

 //at this point we have what we need to start interpolating eme50 angles
 //data proving linear interp is ok is in angles.txt

 
 strcpy(filename,"out.txt");
 out=fopen(filename,"w");// open out file
 for(int k=0;k<j;k++)
   {//find ra,dec,and twist at the EPD time
     ra=interpolate(1,0,motor[k][0]);
     dec=interpolate(1,1,motor[k][0]);
     // get twist angle for epd (outlined in my dissertation)
     s[0]=cos(ra)*cos(dec);
     s[1]=sin(ra)*cos(dec);
     s[2]=sin(dec);

     n[0]=s[1];
     n[1]=-s[0];
     n[2]=0;
     SpiceDouble p[3];
     
     vcrss_c (s,e,p );
     twist=vsep_c(n,p);
     if(vdot_c(e,n)<0.0)
       twist=-twist;
     
  


     
     //given the eme50 angles get the transform RJ19502SC
     R19502SC(ra,dec,twist,RJ19502SC);
     //transpose this matrix to RSC2J1950
     xpose_c(RJ19502SC,RSC2J1950);
   
     //get the rotation matrix from fk4 to j200: RJ19502J2000
     pxform_c ( "fk4", "J2000",   motor[k][0],    RJ19502J2000);
     //where is the galileo sc in relation to Ganymede
     spkezp_c(sc,motor[k][0],body,"NONE",sat,state,&lt);
     
    
     //find the transform from J2000 to IAU_MOON
   
     pxform_c ( "J2000", body, motor[k][0] , RJ20002GSII); 
    

     //at this point I have all the transforms
     //next calc look direction in SC if epd is not behind the shield
     if((int)motor[k][1]!=0){
       //look direction is SC coordinates based on motor and sector
     look[0]=cos(phi0+deg2rad*sector[(int)motor[k][2]])*sin(deg2rad*Mposition[(int)motor[k][1]]);
     look[1]=sin(phi0+deg2rad*sector[(int)motor[k][2]])*sin(deg2rad*Mposition[(int)motor[k][1]]);
     look[2]=cos(deg2rad*Mposition[(int)motor[k][1]]);

     //start the transforms
     mxv_c(RSC2J1950,look,look);
     mxv_c(RJ19502J2000,look,look);
     mxv_c(RJ20002GSII,look,look);
     //look is now in IAU_MOON

     //find the decimal time of utc for graphing purposes to compare to look directions calculated in 1996 and 1997
     et2utc_c ( motor[k][0] , "C", 3, lenout, UTC);
     memcpy( subbuff, &UTC[12], 2 );
     subbuff[2] = '\0';
     time=0.0;
     time=atoi(subbuff)*1.0;
      memcpy( subbuff, &UTC[15], 2 );
     subbuff[2] = '\0';
      time=time+(float)atoi(subbuff)/60. ;  
     memcpy( subbuff, &UTC[18], 2 );
     subbuff[2] = '\0';
     
     double temp=(double) atoi(subbuff);
      memcpy( subbuff, &UTC[21], 3 );
     subbuff[3] = '\0';
     temp=temp+(double)atoi(subbuff)/1000.0;
     time=time+temp/3600.;
     //end find the time
     
    

     //out put results
     fprintf(out," %f %d %d %f %f %f %f %f %f %f\n",time,(int)motor[k][1],(int)motor[k][2],look[0],look[1],look[2],state[0],state[1],state[2],counts[k]);
          
     }
   }
 fclose(out);
}


/*
This a simple linear interpolation function.  I graphed the angle data and it is clearly linear (barring a bit of bit noise for all angles).
WARNING: this is not a robust function and was written only for the G2 encounter

i is 0 for increasing angles and 1 is for decreasing (this takes care of the radian rollover from 2pi to 0.)
a is for the angle wanted (0->ra,1->dec,2->twist)
time1 is the epd ephemeris time
the interpolated angle is returned
twist in this context is for mag and not for epd
 */
float interpolate(int i,int a,double time1)
{
  //find indexes this time1 is between
  int top=0;
  double t,b,pi2=3.14159*2;
  //find the mag time interval that the epd time is between
  for (int k=0;k<max;k++)
    {
      t=angles[k][0];//top
      b=angles[k+1][0];//bottom
      if(time1>angles[k][0]&&time1<=angles[k+1][0])
	{//if epd time is inbetween then note it and break
	  top=k;
	  break;
	}
    }
  //lineraly interpolate between those times
  double v1,v2,t1,t2,m,vo;
  v1=angles[top][a+1];
  v2=angles[top+1][a+1];
  t1=angles[top][0];
  t2=angles[top+1][0];
 
  if(i==0)
    {
      if(v2<v1)
	{
	  v2=v2+pi2;
	}
    }
  m=(v2-v1);
 
  double m2=(t2-t1);
  m=m/m2;
  vo=v1+m*(t1-time1);
  if(i==0) 
    {
      if(vo>pi2) vo=vo-pi2;
    }
  return vo;

}
//given the right acension, declination, and twist angles return the transformation matrix J1950->SC 
void R19502SC(double a,double d,double t,double (*R)[3])
{
  

  R[0][0]=cos(t)*sin(d)*cos(a)-sin(t)*sin(a);
  R[0][1]=cos(t)*sin(d)*sin(a)+sin(t)*cos(a);
  R[0][2]=-cos(t)*cos(d);

  R[1][0]=-sin(t)*sin(d)*cos(a)-cos(t)*sin(a);
  R[1][1]=-sin(t)*sin(d)*sin(a)+cos(t)*cos(a);
  R[1][2]=sin(t)*cos(d);

  R[2][0]=cos(d)*cos(a);
  R[2][1]=cos(d)*sin(a);
  R[2][2]=sin(d);

  
}
