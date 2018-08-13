#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "smoke_sim.h"
#include "constants.h"
#include "open_gl_headers.h" 
//#include "stb_image_write.h"
#include "custom_output.h"
#include "basic_math.h"
#include "mac_grid.h"
#include <fstream>

SmokeSim::SmokeSim() : mFrameNum(0), mTotalFrameNum(0)
{
   reset();
}

SmokeSim::~SmokeSim()
{
}

// Parse given concentration data
void SmokeSim::fillConcVec(std::vector<std::vector<vec4>>& c) {
  #ifdef TEST
  // Dummy loop that fills in random values (for testing only)
  for (int i = 0; i < nFrames; ++i) {
	    std::vector<vec4> internal;
	    for (int j = 0; j < 9; ++j) {
		    internal.push_back(vec4(rand() % 10 + 2, rand() % 10 + 2, rand() % 10 + 2, rand() % 10 + 2));
	    }
	    c.push_back(internal);
  }

  return;
  #endif // TEST

    std::ifstream inFile;
    inFile.open("../concData.csv");

    std::string currLine;
    std::string delim = ",";
    int lineNo = 0;

    std::vector<vec4> currConc = {vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0),
                                  vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0)};

    while (std::getline(inFile, currLine)) {
        if (lineNo != 0) {
            // Chemical name
            std::string token = currLine.substr(0, currLine.find(delim));
            int chemID;
            switch (token.length()) {
                // Methylosmolene
                case 14:
                    chemID = MACGrid::Methyl;
                    break;
                // Chlorodinine
                case 12:
                    chemID = MACGrid::Chloro;
                    break;
                // AGOC-3A
                case 7:
                    chemID = MACGrid::AGOC;
                    break;
                // Appluimonia
                case 11:
                    chemID = MACGrid::App;
                    break;
                // Unrecognized chemical -- ignore
                default:
                    throw std::invalid_argument("Received unrecognized chemical");
            }
            currLine.erase(0, currLine.find(delim) + delim.length());

            // Sensor number
            token = currLine.substr(0, currLine.find(delim));
            int sensorID = stoi(token) - 1;
            currLine.erase(0, currLine.find(delim) + delim.length());

            // Skip the date/time, as it's irrelevant
            currLine.erase(0, currLine.find(delim) + delim.length());

            // Write concentration to the appropriate location
            token = currLine.substr(0, currLine.find(delim));
            currConc[sensorID][chemID] = stof(token);

            // Only push the vector when all 36 entries for the timestamp have been parsed then reset the vector
            if (lineNo % 36 == 0) {
                c.push_back(currConc);
                currConc = {vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0),
                            vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0)};
            }
        }

        lineNo++;
    }
    inFile.close();
}

// Parse given wind data
void SmokeSim::fillWindVec(std::vector<vec3>& w) {
  #ifdef TEST
  // Dummy loop that fills in random values (for testing only)
  for (int i = 0; i < nFrames / 3; ++i) {
      w.push_back(vec3(rand() % 8 + 3, rand() % 360, 0));
  }

  return;
  #endif // TEST

    std::ifstream inFile;
    inFile.open("../windData.csv");

    std::string currLine;
    std::string delim = ",";
    int lineNo = 0;
    while (std::getline(inFile, currLine)) {
        if (lineNo != 0) {
            vec3 currWind;

            // Last entry is unused
            currWind[2] = 0;

            // Erase the date/time info, as it's irrelevant; the data is poorly taken
            currLine.erase(0, currLine.find(delim) + delim.length());

            // Direction
            std::string token = currLine.substr(0, currLine.find(delim));
            currWind[1] = stof(token);
            currLine.erase(0, currLine.find(delim) + delim.length());

            // Magnitude
            token = currLine.substr(0, currLine.find(delim) + delim.length());
            currWind[0] = stof(token);

            w.push_back(currWind);
        }

        lineNo++;
    }
    inFile.close();
}

void SmokeSim::reset()
{
   mGrid.reset();
   mTotalFrameNum = 0;
   fillConcVec(concVec);
   fillWindVec(windVec);
}

void SmokeSim::step()
{
    if (mFrameNum == 0) {
        reset();
    }

    // Time step
	double dt = 1.0 / framesPerHour;

   // Gather user forces
   vec3 windBack = windVec.back();
   std::vector<vec4> concBack = concVec.back();
   std::vector<vec4> concZero = {vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0),
                                 vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0), vec4(0,0,0,0)};

   // Apply time scale
   vec3 wind = (mFrameNum % (3 * framesPerHour) == 0) ? (windBack) : (vec3(0, 0, 0)); // Wind is updated every 3 hours
   std::vector<vec4> conc = (mFrameNum % framesPerHour == 0) ? (concBack) : (concZero); // Concentration is updated every hour

   // Main update call
    mGrid.updateSources(wind, conc, mFrameNum);

   // Adjust wind and concentration vectors if necessary
   if (mFrameNum % (3 * framesPerHour) == 0 && !windVec.empty()) {
       windVec.pop_back();
   }
   if (mFrameNum % framesPerHour == 0 && !concVec.empty()) {
       concVec.pop_back();
   }

   // Calculate new velocities
   mGrid.advectVelocity(dt);

   // Calculate new concentrations and apply mass diffusion
   mGrid.advectConcentration(dt);
   mGrid.diffuse(dt);

   // Pressure projection
   mGrid.project(dt);

   // Calculate new temperature
   mGrid.advectTemperature(dt);

   // Calculate new density
   mGrid.advectDensity(dt);

   // Advect rendering particles
   mGrid.advectRenderingParticles(dt);

   grabScreen();

   mTotalFrameNum++;
}

void SmokeSim::grabScreen()  // Code adapted from asst#1 . USING STB_IMAGE_WRITE INSTEAD OF DEVIL.
{
	if (mFrameNum > nFrames) {
        std::cout << "\n" << "All frames written. Writing output data..." << std::endl;
        mGrid.writeAmountProduced();
        mGrid.writeTravelTimes(writeActualTimes);
        exit(0);
    }

	// Write rendering particle data to .bgeo files
    std::cout << "Writing frame " << std::to_string(mFrameNum) << "..." << std::endl;
	std::string particleFile = "../records/frame" + std::to_string(mFrameNum) + ".bgeo";
    std::string fieldFile = "../records/FieldFrame" + std::to_string(mFrameNum) + ".bgeo";

	mGrid.saveParticle(particleFile, mFrameNum == 0);
    mGrid.saveVelocityField(fieldFile);

	mFrameNum++;
}