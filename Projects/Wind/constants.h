

#ifndef CONSTANTS_H
#define CONSTANTS_H



#define LERP(a,b,t) (1-t)*a + t*b



// Don't modify the values of these here.
// Modify the values of these in Constants.cpp instead.
extern const int theDim[3];
extern const double theCellSize;
extern const double theAirDensity;
extern const double theBuoyancyAlpha;
extern const double theBuoyancyBeta;	
extern const double theBuoyancyAmbientTemperature;
extern const double theVorticityEpsilon;

extern const double D_Methyl;
extern const double D_Chloro;
extern const double D_AGOC;
extern const double D_App;

extern const int framesPerHour;
extern const double unitConversion;

extern const int nFrames;

extern const bool writeActualTimes;

#endif