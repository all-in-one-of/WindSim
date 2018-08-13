

#include "constants.h"

#ifdef _DEBUG
const int theDim[3] = {4, 4, 1};
#else
const int theDim[3] = {200, 200, 1}; // Grid as specified by VAST Challenge, representing a 12 mi x 12 mi area
#endif

// Size of each cell, in grid units
const double theCellSize = 1.0;

const double theAirDensity = 1.0;

const double theBuoyancyAlpha = 0.08; // Gravity's effect on the smoke particles (irrelevant in wind sim)
const double theBuoyancyBeta = 0.37; // Buoyancy's effect due to temperature difference (irrelevant in wind sim)
const double theBuoyancyAmbientTemperature = 0.0; // Ambient temperature (irrelevant in wind sim)

const double theVorticityEpsilon = 0.0; // Coefficient for vorticity confinement (irrelevant in wind sim)

// Binary diffusion coefficient of each species into air
    // These chemicals don't appear to be real, so I'll assume they all diffuse into air at the same rate
const double D_Methyl = 1.0;
const double D_Chloro = 1.0;
const double D_AGOC = 1.0;
const double D_App = 1.0;

// Number of frames that represent 1 hour of data
    // Increasing this number makes the simulation slower but more accurate, as more particles will be rendered
const int framesPerHour = 1;

// Unit conversion for the grid
    // Wind speeds in m/s but grid is dimensioned in mi
    // 1 unit == 0.06 mi => 0.06 mi/unit (12 mi divided into 200 units)
    // 1 m/s == 2.237 mi/h
    // 2.237 mi/h * 100/6 unit/mi == 37.28333... unit/hr * (1/framesPerHour hr/frame) == 37.28333/fPH unit/frame
const double unitConversion = 37.28333 / framesPerHour;

// Number of frames to simulate
    // There are 2208 hours of data total; to run the full simulation, set nFrames = 2208 * framesPerHour
const int nFrames = 2208 * framesPerHour;

// Whether the writeTravelTimes() function should write the actual times (true) or estimated times (false)
const bool writeActualTimes = true;