#include "mac_grid.h"
#include "open_gl_headers.h"
#include "custom_output.h" 
#include "constants.h" 
#include "basic_math.h"
#include <math.h>
#include <map>
#include <stdio.h>
#include <cstdlib>
#include <Eigen/Dense>
#undef max
#undef min 
#include <fstream>



// Globals
MACGrid target;
GridData fConfX, fConfY, fConfZ;

// NOTE: x -> cols, z -> rows, y -> stacks

#define FOR_EACH_CELL \
  for (int k = 0; k < theDim[MACGrid::Z]; k++) \
    for(int j = 0; j < theDim[MACGrid::Y]; j++) \
       for(int i = 0; i < theDim[MACGrid::X]; i++) 

#define FOR_EACH_CELL_REVERSE \
  for (int k = theDim[MACGrid::Z] - 1; k >= 0; k--) \
    for(int j = theDim[MACGrid::Y] - 1; j >= 0; j--) \
       for(int i = theDim[MACGrid::X] - 1; i >= 0; i--) 

#define FOR_EACH_FACE \
    for (int k = 0; k < theDim[MACGrid::Z]; k++) \
      for(int j = 0; j < theDim[MACGrid::Y]+1; j++) \
        for(int i = 0; i < theDim[MACGrid::X]+1; i++) 


MACGrid::MACGrid()
{
   initialize();
}

MACGrid::MACGrid(const MACGrid& orig)
{
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;
}

MACGrid& MACGrid::operator=(const MACGrid& orig)
{
   if (&orig == this)
   {
      return *this;
   }
   mU = orig.mU;
   mV = orig.mV;
   mW = orig.mW;
   mP = orig.mP;
   mD = orig.mD;
   mT = orig.mT;   

   return *this;
}

MACGrid::~MACGrid()
{
}

void MACGrid::reset()
{
   mU.initialize();
   mV.initialize();
   mW.initialize();
   mP.initialize();
   mD.initialize();
   mT.initialize(0.0);

   C_Methyl.initialize();
   C_Chloro.initialize();
   C_AGOC.initialize();
   C_App.initialize();

   nMethyl.initialize();
   nChloro.initialize();
   nAGOC.initialize();
   nApp.initialize();

   firstPartialX.initialize();
   firstPartialY.initialize();
   secondPartialXX.initialize();
   secondPartialYY.initialize();

   calculateAMatrix();
   calculatePreconditioner(AMatrix);
}

void MACGrid::initialize()
{
   reset();
}

// Private helper to refresh particles in each source
void MACGrid::refresh(int i, int j, int k, int chemical, int frameNum, int sensorIdx) {
    vec3 cell_center(theCellSize * (i + 0.5), theCellSize * (j + 0.5), 0);

    // Render more particles based on the concentration
    int numToRender;
    //double scale = 0.5;
    switch (chemical) {
        case MACGrid::Methyl :
            numToRender = int(std::ceil(getMethylConcentration(cell_center)));
            break;
        case MACGrid::Chloro :
            numToRender = int(std::ceil(getChloroConcentration(cell_center)));
            break;
        case MACGrid::AGOC :
            numToRender = int(std::ceil(getAGOCConcentration(cell_center)));
            break;
        case MACGrid::App :
            numToRender = int(std::ceil(getAppConcentration(cell_center)));
            break;
        default :
            // Only render factories once
            numToRender = 0;
            break;
    }
    for (int q = 0; q < numToRender; ++q) {
        double a = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
        double b = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
        //double c = ((float) rand() / RAND_MAX - 0.5) * theCellSize;
        vec3 shift(a, b, 0);
        vec3 xp = cell_center + shift;
        rendering_particles.push_back(xp);

        // Color Key: Methyl = cyan, Chloro = magenta, AGOC = yellow, App = white
        vec3 color(chemical != MACGrid::Methyl,
                   chemical != MACGrid::Chloro,
                   chemical != MACGrid::AGOC);
        rendering_particles_col.push_back(color);

        //  rendering_particles_rad.push_back(numToRender);

        // Store when the particle was produced
        rendering_particles_frame.push_back(vec4(frameNum, -1, sensorIdx, -1));
    }
}

// wind is a vec3: {magnitude, direction, 0}
// conc is a vec4: {Methyl, Chloro, AGOC, App}
void MACGrid::compoundSource(const vec3& wind, const vec4& conc, int i, int j, int frameNum, int sensorIdx) {
  // Direction is the angle that the wind *originates from*
    // That is, it's going in the opposite direction wrt North = 0 and CW+
  // Coordinates become (-sin, -cos) for vector pointing from the center of this circle (i.e. - the origin)
    // Flip vector to go backwards in time
  double theta = (270.0 - wind[1]) * BasicMath::PI / 180.0;
  mU(i + 1, j, 0) -= unitConversion * wind[0] * sin(theta);
  mU(i + 2, j, 0) -= unitConversion * wind[0] * sin(theta);
  mV(i, j + 1, 0) -= unitConversion * wind[0] * cos(theta);
  mV(i, j + 2, 0) -= unitConversion * wind[0] * cos(theta);

  mW(i, j, 0) = 0;

  mD(i, j, 0) = 1.0; //0.8
  mT(i, j, 0) = 0.0; //0.8

  C_Methyl(i, j, 0) += conc[0];
  C_Chloro(i, j, 0) += conc[1];
  C_AGOC(i, j, 0) += conc[2];
  C_App(i, j, 0) += conc[3];

  MACGrid::refresh(i, j, 0, MACGrid::Methyl, frameNum, sensorIdx);
  MACGrid::refresh(i, j, 0, MACGrid::Chloro, frameNum, sensorIdx);
  MACGrid::refresh(i, j, 0, MACGrid::AGOC, frameNum, sensorIdx);
  MACGrid::refresh(i, j, 0, MACGrid::App, frameNum, sensorIdx);

  // Check all particles for "collisions" with all factories at every frame
  for (int q = 0; q < rendering_particles_frame.size(); ++q) {
    // Check if the particle reached the factory; the second condition rules out the factories themselves
    int factoryId = reachesFactory(rendering_particles[q]);
    if (factoryId != -1 && rendering_particles_col[q] != vec3(1, 0.5333, 0)) {
        rendering_particles_frame[q][1] = frameNum;
        rendering_particles_frame[q][3] = factoryId;
    }
  }
}

// Helper to check if a given rendered particle reaches any factory
// Returns the index of that factory if so and -1 if not
    // Note that some factories may not be considered "reached" even if they have been predicted to output chemicals,
    // since this only makes use of *rendered* particles; thus, the chemical output is more accurate.
int MACGrid::reachesFactory(vec3 particle) {
    // Check the square surrounding each factory
    for (int k = 0; k < factoryPos.size(); ++k) {
        vec3 dist = particle - factoryPos[k];
        // n * theCellSize * sqrt(2) checks all points in an (n+1)x(n+1) square with the factory at the center
            // Note that n doesn't necessarily have to be an integer
        if (dist.Length() <= 2 * theCellSize * sqrt(2)) {
            return k;
        }
    }
    return -1;
}

// Treat sensors as sources and factories as sensors and run a reverse-time simulation
void MACGrid::updateSources(const vec3& wind, const std::vector<vec4>& conc, int frameNum) {
  // Update all sources
  for (int i = 0; i < sensorPos.size(); ++i) {
      MACGrid::compoundSource(wind, conc[i], int(sensorPos[i][0]), int(sensorPos[i][1]), frameNum, i);
  }
}

void MACGrid::advectVelocity(double dt)
{
    // Make sure the target values are initialized and retain this value in the case of invalid faces
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

    FOR_EACH_FACE {
        // Get previous positions of particles that are currently at the relevant x, y, z faces
        vec3 oldPosX, oldPosY;
        if (isValidFace(MACGrid::X, i, j, k)) {
            oldPosX = getRewoundPosition(getFacePosition(MACGrid::X, i, j, k), dt);

      			// Assign the velocity these particles used to have as the current velocity
      			target.mU(i, j, k) = getVelocityX(oldPosX);
        }
        if (isValidFace(MACGrid::Y, i, j, k)) {
            oldPosY = getRewoundPosition(getFacePosition(MACGrid::Y, i, j, k), dt);
			      target.mV(i, j, k) = getVelocityY(oldPosY);
        }
        target.mW(i, j, k) = 0;
    }

    // Then save the result to our object
    mU = target.mU;
    mV = target.mV;
    mW = target.mW;
}

void MACGrid::advectTemperature(double dt)
{
    target.mT = mT;

    FOR_EACH_CELL {
        // Find the old temperature of the particle at the center of this cell
        double oldTemp = 0;
        if (isValidCell(i, j, k)) {
            oldTemp = getTemperature(getRewoundPosition(getCenter(i, j, k), dt));

      			// Assign it as the current temperature
      			target.mT(i, j, k) = oldTemp;
        }
    }

    // Then save the result to our object
    mT = target.mT;
}

void MACGrid::advectConcentration(double dt) {
  target.C_Methyl = C_Methyl;
  target.C_Chloro = C_Chloro;
  target.C_AGOC = C_AGOC;
  target.C_App = C_App;

  FOR_EACH_CELL {
    double oldMethylConc = 0;
    double oldChloroConc = 0;
    double oldAGOCConc = 0;
    double oldAppConc = 0;

    if (isValidCell(i, j, k)) {
      oldMethylConc = getMethylConcentration(getRewoundPosition(getCenter(i, j, k), dt));
      oldChloroConc = getChloroConcentration(getRewoundPosition(getCenter(i, j, k), dt));
      oldAGOCConc = getAGOCConcentration(getRewoundPosition(getCenter(i, j, k), dt));
      oldAppConc = getAppConcentration(getRewoundPosition(getCenter(i, j, k), dt));

      target.C_Methyl(i, j, k) = oldMethylConc;
      target.C_Chloro(i, j, k) = oldChloroConc;
      target.C_AGOC(i, j, k) = oldAGOCConc;
      target.C_App(i, j, k) = oldAppConc;
    }
  }

  C_Methyl = target.C_Methyl;
  C_Chloro = target.C_Chloro;
  C_AGOC = target.C_AGOC;
  C_App = target.C_App;
}


void MACGrid::advectRenderingParticles(double dt) {
	rendering_particles_vel.resize(rendering_particles.size());
	for (size_t p = 0; p < rendering_particles.size(); p++) {
        // Prevent factory particles from advecting
        if (rendering_particles_col[p] == vec3(1, 0.5333, 0)) {
            continue;
        }
		vec3 currentPosition = rendering_particles[p];
        vec3 currentVelocity = getVelocity(currentPosition);
        vec3 nextPosition = currentPosition + currentVelocity * dt;
        vec3 clippedNextPosition = clipToGrid(nextPosition, currentPosition);
        // Keep going...
        vec3 nextVelocity = getVelocity(clippedNextPosition);
        vec3 averageVelocity = (currentVelocity + nextVelocity) / 2.0;
        vec3 betterNextPosition = currentPosition + averageVelocity * dt;
        vec3 clippedBetterNextPosition = clipToGrid(betterNextPosition, currentPosition);
        rendering_particles[p] = clippedBetterNextPosition;
		    rendering_particles_vel[p] = averageVelocity;
	}
}

void MACGrid::advectDensity(double dt)
{
    // Same process as with temperature
    FOR_EACH_CELL {
        double oldDensity = 0;
        if (isValidCell(i, j, k)) {
            oldDensity = getDensity(getRewoundPosition(getCenter(i, j, k), dt));
			      target.mD(i, j, k) = oldDensity;
        }
    }

    // Then save the result to our object
    mD = target.mD;
}

void MACGrid::computeBuoyancy(double dt)
{
    target.mV = mV;

    GridData buoyancy;
    buoyancy.initialize();
    FOR_EACH_CELL {
        if (isValidFace(MACGrid::Y, i, j, k) && j > 0) {
            buoyancy(i, j, k) = -theBuoyancyAlpha * mD(i, j, k) +
                                theBuoyancyBeta * (mT(i, j, k) - theBuoyancyAmbientTemperature);
        }
    }

    FOR_EACH_CELL {
        double fBuoy = 0;
        if (isValidFace(MACGrid::Y, i, j, k) && j > 0) {
            fBuoy = buoyancy.interpolate(getCenter(i, j, k));
        }

        target.mV(i, j, k) += fBuoy;
    }

    // and then save the result to our object
    mV = target.mV;
}

// Helper function to compute omega, the vorticity vector
vec3 MACGrid::computeOmega(int i, int j, int k) {
    //double oX = (mW(i, j + 1, k) - mW(i, j - 1, k)) / (2 * theCellSize) - (mV(i, j, k + 1) - mV(i, j, k - 1)) / (2 * theCellSize); // Should be 0
    //double oY = (mU(i, j, k + 1) - mU(i, j, k - 1)) / (2 * theCellSize) - (mW(i + 1, j, k) - mW(i - 1, j, k)) / (2 * theCellSize); // Should be 0
    double oZ = (mV(i + 1, j, k) - mV(i - 1, j, k)) / (2 * theCellSize) - (mU(i, j + 1, k) - mU(i, j - 1, k)) / (2 * theCellSize); // Nonzero

	return vec3(0, 0, oZ); // (0, 0, oZ)
}

void MACGrid::computeVorticityConfinement(double dt)
{
	// Calculate and store omega
	GridData omegaLength;
	omegaLength.initialize();

  FOR_EACH_CELL {
		omegaLength(i, j, k) = computeOmega(i, j, k).Length();
	}

	// Calculate and store the confinement force
	fConfX.initialize();
	fConfY.initialize();
    fConfZ.initialize();

	FOR_EACH_CELL {
		double gOX = (omegaLength(i + 1, j, k) - omegaLength(i - 1, j, k)) / (2 * theCellSize);
		double gOY = (omegaLength(i, j + 1, k) - omegaLength(i, j - 1, k)) / (2 * theCellSize);

        double oZ = computeOmega(i, j, k)[2];

    // Cross product is (gOX, gOY, 0) x (0, 0, oZ) = (gOY*oZ, -gOX*oZ, 0)
    // This makes sense, as there shouldn't be any confinement in the z direction
		fConfX(i, j, k) = theVorticityEpsilon * theCellSize * gOY * oZ; // Nonzero
		fConfY(i, j, k) = -theVorticityEpsilon * theCellSize * gOX * oZ; // Nonzero
        fConfZ(i, j, k) = 0; // Should be 0
	}

	// Update velocities
    // Important to ensure the values are initialized to the old values of mU, mV, mW so we're not adding to
    // garbage values.
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

	FOR_EACH_FACE {
		if (isValidFace(MACGrid::X, i, j, k)) {
			target.mU(i, j, k) += fConfX.interpolate(getFacePosition(MACGrid::X, i, j, k));
		}
		if (isValidFace(MACGrid::Y, i, j, k)) {
			target.mV(i, j, k) += fConfY.interpolate(getFacePosition(MACGrid::Y, i, j, k));
		}
        target.mW(i, j, k) = 0;
	}

   // Then save the result to our object
   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}

void MACGrid::addExternalForces(double dt)
{
   computeBuoyancy(dt);
   computeVorticityConfinement(dt);
}

void MACGrid::project(double dt)
{
    // Again, make sure we won't be subtracting from garbage values
    target.mU = mU;
    target.mV = mV;
    target.mW = mW;

	// Check border velocities
    FOR_EACH_CELL {
        if (i == 0) {
            mU(i, j, k) = 0;
        }
        if (!isValidCell(i + 1, j, k)) {
            mU(i + 1, j, k) = 0;
        }
        if (j == 0) {
            mV(i, j, k) = 0;
        }
        if (!isValidCell(i, j + 1, k)) {
            mV(i, j + 1, k) = 0;
        }

        mW(i, j, k) = 0;
    }

	// Construct b (divergences)
    GridData b;
    b.initialize();
	// Adapted from DEBUG block
	FOR_EACH_CELL {
        double div = 0;
        if (isValidCell(i, j, k)) {
            double velLowX = mU(i, j, k);
            double velHighX = mU(i + 1, j, k);
            double velLowY = mV(i, j, k);
            double velHighY = mV(i, j + 1, k);

            div = -((velHighX - velLowX) + (velHighY - velLowY)) / theCellSize;
        }

        b(i, j, k) = div;
	}

	// A is already calculated as AMatrix; solve eqn for P
	GridData P;
	P.initialize();
	preconditionedConjugateGradient(AMatrix, P, b, 100, 1e-5);

	// Set target.mP to an adjusted value of P
	double adjustP = theAirDensity * theCellSize * theCellSize / dt;
	FOR_EACH_CELL {
		target.mP(i, j, k) = adjustP * P(i, j, k);
	}

	// Subtract P from velocities, taking boundaries into account
	double adjustV = theAirDensity * theCellSize / dt;
	FOR_EACH_FACE {
		double p0, pF;

		if (isValidFace(MACGrid::X, i, j, k)) {
			if (i - 1 < 0) {
				p0 = target.mP(i, j, k) - adjustV * mU(i, j, k);
			}
			else {
				p0 = target.mP(i - 1, j, k);
			}

			if (i >= theDim[MACGrid::X]) {
				pF = target.mP(i - 1, j, k) + adjustV * mU(i, j, k);
			}
			else {
				pF = target.mP(i, j, k);
			}

			target.mU(i, j, k) -= (pF - p0) / adjustV;
		}

		if (isValidFace(MACGrid::Y, i, j, k)) {
			if (j - 1 < 0) {
				p0 = target.mP(i, j, k) - adjustV * mV(i, j, k);
			}
			else {
				p0 = target.mP(i, j - 1, k);
			}

			if (j >= theDim[MACGrid::Y]) {
				pF = target.mP(i, j - 1, k) + adjustV * mV(i, j, k);
			}
			else {
				pF = target.mP(i, j, k);
			}

			target.mV(i, j, k) -= (pF - p0) / adjustV;
		}

		if (isValidFace(MACGrid::Z, i, j, k)) {
			if (k - 1 < 0) {
				p0 = target.mP(i, j, k) - adjustV * mW(i, j, k);
			}
			else {
				p0 = target.mP(i, j, k - 1);
			}

			if (k >= theDim[MACGrid::Z]) {
				pF = target.mP(i, j, k - 1) + adjustV * mW(i, j, k);
			}
			else {
				pF = target.mP(i, j, k);
			}

			target.mW(i, j, k) -= (pF - p0) / adjustV;
		}
	}

  // Check border velocities again
  FOR_EACH_CELL {
      if (i == 0) {
          mU(i, j, k) = 0;
      }
      if (!isValidCell(i + 1, j, k)) {
          mU(i + 1, j, k) = 0;
      }
      if (j == 0) {
          mV(i, j, k) = 0;
      }
      if (!isValidCell(i, j + 1, k)) {
          mV(i, j + 1, k) = 0;
      }

      mW(i, j, k) = 0;
  }

   // Then save the result to our object
   mP = target.mP; 
   mU = target.mU;
   mV = target.mV;
   mW = target.mW;
}

// Approximates the first partial derivative with a difference method (only supports X and Y for purposes of wind sim)
// dir is the direction the derivative will be taken with respect to (i.e. - MACGrid::X or MACGrid::Y)
void MACGrid::computeFirstPartial(const GridData& C, int dir) {
  if (dir == MACGrid::X) {
    FOR_EACH_CELL {
      double cPlus = (isValidCell(i + 1, j, k)) ? (C(i + 1, j, k)) : (0.0);
      double cMinus = (isValidCell(i - 1, j, k)) ? (C(i - 1, j, k)) : (0.0);
      firstPartialX(i, j, k) = cPlus - cMinus;
    }
  }
  else if (dir == MACGrid::Y) {
    FOR_EACH_CELL {
      double cPlus = (isValidCell(i, j + 1, k)) ? (C(i, j + 1, k)) : (0.0);
      double cMinus = (isValidCell(i, j - 1, k)) ? (C(i, j - 1, k)) : (0.0);
      firstPartialY(i, j, k) = cPlus - cMinus;
    }
  }
}

// Approximates the second partial derivative using the above method
// Assumes no mixed partials will be taken
void MACGrid::computeSecondPartial(const GridData& C, int dir) {
  // Fill firstPartial GridData
  computeFirstPartial(C, dir);

  if (dir == MACGrid::X) {
    FOR_EACH_CELL {
      double dPlus = (isValidCell(i + 1, j, k)) ? (firstPartialX(i + 1, j, k)) : (0.0);
      double dMinus = (isValidCell(i - 1, j, k)) ? (firstPartialX(i - 1, j, k)) : (0.0);
      secondPartialXX(i, j, k) = dPlus - dMinus;
    }
  }
  else if (dir == MACGrid::Y) {
    FOR_EACH_CELL {
      double dPlus = (isValidCell(i, j + 1, k)) ? (firstPartialY(i, j + 1, k)) : (0.0);
      double dMinus = (isValidCell(i, j - 1, k)) ? (firstPartialY(i, j - 1, k)) : (0.0);
      secondPartialYY(i, j, k) = dPlus - dMinus;
    }
  }
}

// Computes chemical production using the law of conservation of species
// Stored at each cell's center
void MACGrid::applyMassDiffusion(double dt, const GridData& C, double D, GridData& n) {
    // Compute second partial derivatives of C, which in turn computes their first partial derivatives
  computeSecondPartial(C, MACGrid::X);
  computeSecondPartial(C, MACGrid::Y);

  FOR_EACH_CELL {
    // Rate of mass production is nDot = mU * dC/dx + mV * dC/dy - D * (d2C/dx2 + d2C/dy2)
    // Integral of nDot wrt time will be the total amount of <chemical> produced at each location
      // Do this by summing up nDot * dt then extract the values at each factory's location
    n(i, j, k) -= (mU(i, j, k) * firstPartialX(i, j, k) + mV(i, j, k) * firstPartialY(i, j, k) - D * (secondPartialXX(i, j, k) + secondPartialYY(i, j, k))) * dt;
  }
}

void MACGrid::diffuse(double dt) {
    applyMassDiffusion(dt, C_Methyl, D_Methyl, nMethyl);
    applyMassDiffusion(dt, C_Chloro, D_Chloro, nChloro);
    applyMassDiffusion(dt, C_AGOC, D_AGOC, nAGOC);
    applyMassDiffusion(dt, C_App, D_App, nApp);
}

vec3 MACGrid::getVelocity(const vec3& pt)
{
   vec3 vel;
   vel[0] = getVelocityX(pt); 
   vel[1] = getVelocityY(pt); 
   vel[2] = getVelocityZ(pt); 
   return vel;
}

double MACGrid::getVelocityX(const vec3& pt)
{
   return mU.interpolate(pt);
}

double MACGrid::getVelocityY(const vec3& pt)
{
   return mV.interpolate(pt);
}

double MACGrid::getVelocityZ(const vec3& pt)
{
   return mW.interpolate(pt);
}

double MACGrid::getTemperature(const vec3& pt)
{
   return mT.interpolate(pt);
}

double MACGrid::getDensity(const vec3& pt)
{
   return mD.interpolate(pt);
}

double MACGrid::getMethylConcentration(const vec3& pt) {
  return C_Methyl.interpolate(pt);
}

double MACGrid::getChloroConcentration(const vec3& pt) {
  return C_Chloro.interpolate(pt);
}

double MACGrid::getAGOCConcentration(const vec3& pt) {
  return C_AGOC.interpolate(pt);
}

double MACGrid::getAppConcentration(const vec3& pt) {
  return C_App.interpolate(pt);
}

vec3 MACGrid::getCenter(int i, int j, int k)
{
   double xstart = theCellSize/2.0;
   double ystart = theCellSize/2.0;
   double zstart = theCellSize/2.0;

   double x = xstart + i*theCellSize;
   double y = ystart + j*theCellSize;
   double z = zstart + k*theCellSize;
   return vec3(x, y, z);
}


vec3 MACGrid::getRewoundPosition(const vec3 & currentPosition, const double dt) {

	/*
	// EULER (RK1):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	return clippedRewoundPosition;
	*/

	// HEUN / MODIFIED EULER (RK2):
	vec3 currentVelocity = getVelocity(currentPosition);
	vec3 rewoundPosition = currentPosition - currentVelocity * dt;
	vec3 clippedRewoundPosition = clipToGrid(rewoundPosition, currentPosition);
	// Keep going...
	vec3 rewoundVelocity = getVelocity(clippedRewoundPosition);
	vec3 averageVelocity = (currentVelocity + rewoundVelocity) / 2.0;
	vec3 betterRewoundPosition = currentPosition - averageVelocity * dt;
	vec3 clippedBetterRewoundPosition = clipToGrid(betterRewoundPosition, currentPosition);
	return clippedBetterRewoundPosition;

}


vec3 MACGrid::clipToGrid(const vec3& outsidePoint, const vec3& insidePoint) {
	/*
	// OLD:
	vec3 rewindPosition = outsidePoint;
	if (rewindPosition[0] < 0) rewindPosition[0] = 0; // TEMP!
	if (rewindPosition[1] < 0) rewindPosition[1] = 0; // TEMP!
	if (rewindPosition[2] < 0) rewindPosition[2] = 0; // TEMP!
	if (rewindPosition[0] > theDim[MACGrid::X]) rewindPosition[0] = theDim[MACGrid::X]; // TEMP!
	if (rewindPosition[1] > theDim[MACGrid::Y]) rewindPosition[1] = theDim[MACGrid::Y]; // TEMP!
	if (rewindPosition[2] > theDim[MACGrid::Z]) rewindPosition[2] = theDim[MACGrid::Z]; // TEMP!
	return rewindPosition;
	*/

	vec3 clippedPoint = outsidePoint;

	for (int i = 0; i < 3; i++) {
		if (clippedPoint[i] < 0) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = 0 - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
		if (clippedPoint[i] > getSize(i)) {
			vec3 distance = clippedPoint - insidePoint;
			double newDistanceI = getSize(i) - insidePoint[i];
			double ratio = newDistanceI / distance[i];
			clippedPoint = insidePoint + distance * ratio;
		}
	}

#ifdef _DEBUG
	// Make sure the point is now in the grid:
	if (clippedPoint[0] < 0 || clippedPoint[1] < 0 || clippedPoint[2] < 0 || clippedPoint[0] > getSize(0) || clippedPoint[1] > getSize(1) || clippedPoint[2] > getSize(2)) {
		PRINT_LINE("WARNING: Clipped point is outside grid!");
	}
#endif

	return clippedPoint;

}


double MACGrid::getSize(int dimension) {
	return theDim[dimension] * theCellSize;
}


bool MACGrid::isValidCell(int i, int j, int k)
{
	if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
		return false;
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


bool MACGrid::isValidFace(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		if (i > theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 1) {
		if (i >= theDim[MACGrid::X] || j > theDim[MACGrid::Y] || k >= theDim[MACGrid::Z]) {
			return false;
		}
	} else if (dimension == 2) {
		if (i >= theDim[MACGrid::X] || j >= theDim[MACGrid::Y] || k > theDim[MACGrid::Z]) {
			return false;
		}
	}

	if (i < 0 || j < 0 || k < 0) {
		return false;
	}

	return true;
}


vec3 MACGrid::getFacePosition(int dimension, int i, int j, int k)
{
	if (dimension == 0) {
		return vec3(i * theCellSize, (j + 0.5) * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 1) {
		return vec3((i + 0.5) * theCellSize, j * theCellSize, (k + 0.5) * theCellSize);
	} else if (dimension == 2) {
		return vec3((i + 0.5) * theCellSize, (j + 0.5) * theCellSize, k * theCellSize);
	}

	return vec3(0,0,0); //???

}

void MACGrid::calculateAMatrix() {

	FOR_EACH_CELL {

		int numFluidNeighbors = 0;
		if (i-1 >= 0) {
			AMatrix.plusI(i-1,j,k) = -1;
			numFluidNeighbors++;
		}
		if (i+1 < theDim[MACGrid::X]) {
			AMatrix.plusI(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (j-1 >= 0) {
			AMatrix.plusJ(i,j-1,k) = -1;
			numFluidNeighbors++;
		}
		if (j+1 < theDim[MACGrid::Y]) {
			AMatrix.plusJ(i,j,k) = -1;
			numFluidNeighbors++;
		}
		if (k-1 >= 0) {
			AMatrix.plusK(i,j,k-1) = -1;
			numFluidNeighbors++;
		}
		if (k+1 < theDim[MACGrid::Z]) {
			AMatrix.plusK(i,j,k) = -1;
			numFluidNeighbors++;
		}
		// Set the diagonal:
		AMatrix.diag(i,j,k) = numFluidNeighbors;
	}
}

bool MACGrid::preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance) {
	// Solves Ap = d for p.

	FOR_EACH_CELL {
		p(i,j,k) = 0.0; // Initial guess p = 0.	
	}

	GridData r = d; // Residual vector.

	/*
	PRINT_LINE("r: ");
	FOR_EACH_CELL {
		PRINT_LINE(r(i,j,k));
	}
	*/
	GridData z; z.initialize();
	applyPreconditioner(r, A, z); // Auxillary vector.
	/*
	PRINT_LINE("z: ");
	FOR_EACH_CELL {
		PRINT_LINE(z(i,j,k));
	}
	*/

	GridData s = z; // Search vector;

	double sigma = dotProduct(z, r);

	for (int iteration = 0; iteration < maxIterations; iteration++) {

		double rho = sigma; // According to TA. Here???

		apply(A, s, z); // z = applyA(s);

		double alpha = rho/dotProduct(z, s);

		GridData alphaTimesS; alphaTimesS.initialize();
		multiply(alpha, s, alphaTimesS);
		add(p, alphaTimesS, p);
		//p += alpha * s;

		GridData alphaTimesZ; alphaTimesZ.initialize();
		multiply(alpha, z, alphaTimesZ);
		subtract(r, alphaTimesZ, r);
		//r -= alpha * z;

		if (maxMagnitude(r) <= tolerance) {
			//PRINT_LINE("PCG converged in " << (iteration + 1) << " iterations.");
			return true; //return p;
		}

		applyPreconditioner(r, A, z); // z = applyPreconditioner(r);

		double sigmaNew = dotProduct(z, r);

		double beta = sigmaNew / rho;

		GridData betaTimesS; betaTimesS.initialize();
		multiply(beta, s, betaTimesS);
		add(z, betaTimesS, s);
		//s = z + beta * s;

		sigma = sigmaNew;
	}

	PRINT_LINE( "PCG didn't converge!" );
	return false;

}


void MACGrid::calculatePreconditioner(const GridDataMatrix & A) {

	precon.initialize();

    // Filling in precon(i,j,k) for all cells
    FOR_EACH_CELL {
        double tau = 0.97;
        double e = 0;
        if (A.diag(i, j, k) != 0.0) {
            e = A.diag(i, j, k) -
                    (A.plusI(i - 1, j, k) * precon(i - 1, j, k)) * (A.plusI(i - 1, j, k) * precon(i - 1, j, k)) -
                    (A.plusJ(i, j - 1, k) * precon(i, j - 1, k)) * (A.plusJ(i, j - 1, k) * precon(i, j - 1, k)) -
                    (A.plusK(i, j, k - 1) * precon(i, j, k - 1)) * (A.plusK(i, j, k - 1) * precon(i, j, k - 1)) -
                    tau * (A.plusI(i - 1, j, k) * (A.plusJ(i - 1, j, k) + A.plusK(i - 1, j, k)) * precon(i - 1, j, k) * precon(i - 1, j, k) +
                          A.plusJ(i, j - 1, k) * (A.plusI(i, j - 1, k) + A.plusK(i, j - 1, k)) * precon(i, j - 1, k) * precon(i, j - 1, k) +
                          A.plusK(i, j, k - 1) * (A.plusI(i, j, k - 1) + A.plusJ(i, j, k - 1)) * precon(i, j, k - 1) * precon(i, j, k - 1));
        }

        precon(i, j, k) = 1 / sqrt(e + 1e-30);
    }

}


void MACGrid::applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z) {
    if(1) {

        // APPLY THE PRECONDITIONER:
        // Solve Lq = r for q:
        GridData q;
        q.initialize();
        FOR_EACH_CELL {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = r(i, j, k) - A.plusI(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k)
                               - A.plusJ(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k)
                               - A.plusK(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                    q(i, j, k) = t * precon(i, j, k);
                    //}
                }
        // Solve L^Tz = q for z:
        FOR_EACH_CELL_REVERSE {
                    //if (A.diag(i,j,k) != 0.0) { // If cell is a fluid.
                    double t = q(i, j, k) - A.plusI(i, j, k) * precon(i, j, k) * z(i + 1, j, k)
                               - A.plusJ(i, j, k) * precon(i, j, k) * z(i, j + 1, k)
                               - A.plusK(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                    z(i, j, k) = t * precon(i, j, k);
                    //}
                }
    }
    else{
        // Unpreconditioned CG: Bypass preconditioner:
        z = r;
        return;
    }

}



double MACGrid::dotProduct(const GridData & vector1, const GridData & vector2) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		result += vector1(i,j,k) * vector2(i,j,k);
	}

	return result;
}


void MACGrid::add(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) + vector2(i,j,k);
	}

}


void MACGrid::subtract(const GridData & vector1, const GridData & vector2, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = vector1(i,j,k) - vector2(i,j,k);
	}

}


void MACGrid::multiply(const double scalar, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL {
		result(i,j,k) = scalar * vector(i,j,k);
	}

}


double MACGrid::maxMagnitude(const GridData & vector) {
	
	double result = 0.0;

	FOR_EACH_CELL {
		if (abs(vector(i,j,k)) > result) result = abs(vector(i,j,k));
	}

	return result;
}


void MACGrid::apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result) {
	
	FOR_EACH_CELL { // For each row of the matrix.

		double diag = 0;
		double plusI = 0;
		double plusJ = 0;
		double plusK = 0;
		double minusI = 0;
		double minusJ = 0;
		double minusK = 0;

		diag = matrix.diag(i,j,k) * vector(i,j,k);
		if (isValidCell(i+1,j,k)) plusI = matrix.plusI(i,j,k) * vector(i+1,j,k);
		if (isValidCell(i,j+1,k)) plusJ = matrix.plusJ(i,j,k) * vector(i,j+1,k);
		if (isValidCell(i,j,k+1)) plusK = matrix.plusK(i,j,k) * vector(i,j,k+1);
		if (isValidCell(i-1,j,k)) minusI = matrix.plusI(i-1,j,k) * vector(i-1,j,k);
		if (isValidCell(i,j-1,k)) minusJ = matrix.plusJ(i,j-1,k) * vector(i,j-1,k);
		if (isValidCell(i,j,k-1)) minusK = matrix.plusK(i,j,k-1) * vector(i,j,k-1);

		result(i,j,k) = diag + plusI + plusJ + plusK + minusI + minusJ + minusK;
	}

}

void MACGrid::saveParticle(std::string filename, bool save){
	Partio::ParticlesDataMutable *parts = Partio::create();
	Partio::ParticleAttribute posH, vH, cH, rH;
	posH = parts->addAttribute("position", Partio::VECTOR, 3);
	vH = parts->addAttribute("v", Partio::VECTOR, 3);
    cH = parts->addAttribute("color", Partio::VECTOR, 3);
    rH = parts->addAttribute("pscale", Partio::FLOAT, 1);

    // Create factories as stationary orange particles
    for (int i = 0; i < 4 && save; ++i) {
        rendering_particles.push_back(factoryPos[i]);
        rendering_particles_vel.push_back(vec3(0, 0, 0));
        rendering_particles_col.push_back(vec3(1, 0.5333, 0));
    }

	for (unsigned int i = 0; i < rendering_particles.size(); i++)
	{
		int idx = parts->addParticle();
		float *p = parts->dataWrite<float>(posH, idx);
		float *v = parts->dataWrite<float>(vH, idx);
        float *c = parts->dataWrite<float>(cH, idx);
        float *r = parts->dataWrite<float>(rH, idx);
		for (int k = 0; k < 3; k++)
		{
			p[k] = rendering_particles[i][k];
			v[k] = rendering_particles_vel[i][k];
            c[k] = rendering_particles_col[i][k];
//            if (k == 0) {
//                r[k] = rendering_particles_rad[i];
//            }
		}
	}
	
	Partio::write(filename.c_str(), *parts);
	parts->release();
}

void MACGrid::saveVelocityField(std::string filename) {
    Partio::ParticlesDataMutable* field = Partio::create();
    Partio::ParticleAttribute posH, vcH, vH;
    posH = field->addAttribute("position", Partio::VECTOR, 3);
    vcH = field->addAttribute("centerV", Partio::VECTOR, 3);
    vH = field->addAttribute("v", Partio::VECTOR, 3);

    FOR_EACH_CELL {
        int idx = field->addParticle();
        float* p = field->dataWrite<float>(posH, idx);
        float* vC = field->dataWrite<float>(vcH, idx);
        float* v = field->dataWrite<float>(vH, idx);

        // The velocity field vectors represent the velocity at the *center* of each cell
        vec3 center = getCenter(i, j, k);
        vec3 centerVel = getVelocity(center);
        for (int l = 0; l < 3; ++l) {
            p[l] = center[l];
            vC[l] = centerVel[l];
            v[l] = 0;
        }
    }

    Partio::write(filename.c_str(), *field);
    field->release();
}

// Outputs the amount predicted to have been produced at each factory to a text file
    // In other words, it writes the cumulative value of N for each chemical at each factory location
void MACGrid::writeAmountProduced() {
    std::ofstream outFile;
    outFile.open("../records/chemicalOutput.txt");

    // Write headers
    outFile << "\t" << "METHYLOSMOLENE \t CHLORODININE \t AGOC-3A \t APPLUIMONIA" << std::endl;

    // Write data of amount produced
    for (int i = 0; i < factoryPos.size(); ++i) {
        auto x = int(factoryPos[i][0]);
        auto y = int(factoryPos[i][1]);
        double methylOut = 0;
        double chloroOut = 0;
        double agocOut = 0;
        double appOut = 0;

        // Treat each factory as occupying a square surrounding its given location
        for (int a = -1; a < 2; ++a) {
            for (int b = -1; b < 2; ++b) {
                methylOut += nMethyl(x + a, y + b, 0);
                chloroOut += nChloro(x + a, y + b, 0);
                agocOut += nAGOC(x + a, y + b, 0);
                appOut += nApp(x + a, y + b, 0);
            }
        }

        // Output the amounts produced -- these must be positive numbers, so take the absolute value
        outFile << factoryName[i] << "\t" << std::to_string(std::abs(methylOut))
                << "\t" << std::to_string(std::abs(chloroOut)) << "\t" << std::to_string(std::abs(agocOut))
                << "\t" << std::to_string(std::abs(appOut)) << std::endl;
    }

    outFile.close();
    std::cout << "Estimated amounts produced, in PPM, have been written to ../records/chemicalOutput.txt as TSV" << std::endl;
}

// Private helper to calculate simple Pythagorean distance between a factory and a sensor
double MACGrid::calculateDistance(int factoryIdx, int sensorIdx) {
    double dx = factoryPos[factoryIdx][0] - sensorPos[sensorIdx][0];
    double dy = factoryPos[factoryIdx][1] - sensorPos[sensorIdx][1];

    return std::sqrt(dx * dx + dy * dy);
}

// Output function to write travel times
void MACGrid::writeTravelTimes(bool calculateActual) {
    std::ofstream outFile;

    // Write estimated travel times (using 5 m/s as velocity ceiling)
    if (!calculateActual) {
        outFile.open("../records/estimatedTravelTimes.txt");

        // Write headers
        outFile << "\t"
                << "Sensor 1 \t Sensor 2 \t Sensor 3 \t Sensor 4 \t Sensor 5 \t Sensor 6 \t Sensor 7 \t Sensor 8 \t Sensor 9"
                << std::endl;

        // Write estimated travel times between factories and sensors
        for (int i = 0; i < factoryPos.size(); ++i) {
            outFile << factoryName[i];
            for (int j = 0; j < sensorPos.size(); ++j) {
                // Note that calculateDistance() returns a distance in *grid units*
                // Multiply function output by 0.06 mi/unit
                // Also multiply by 1609.34 m/mi then divide by 5.0 m/s to get travel time in s
                outFile << "\t" << std::to_string(calculateDistance(i, j) * 0.06 * 1609.34 / 5.0);
            }
            outFile << std::endl;
        }
        outFile.close();
        std::cout << "Estimated travel times, in seconds, have been written to ../records/estimatedTravelTimes.txt as TSV; note" <<
                  " that this assumes that all particles travel in a straight line path from factory to sensor at a " <<
                  "constant velocity of 5 m/s." << std::endl;
    }
    // Calculate the travel times based on rendered particles
    else {
        outFile.open("../records/calcTravelTimes.txt");

        // Matrix of all travel times; rows are factories, columns are sensors
        // NOTE: The dimensions of this matrix *must* be passed as constexpr (I've hardcoded them for simplicity)
        Eigen::Matrix<double, 4, 9> allTimes;
        allTimes.setZero();

        // Number of particles that have traveled between each factory-sensor pair
        // NOTE: This matrix must have the same dimensions as allTimes!
        Eigen::Matrix<int, 4, 9> numDetected;
        numDetected.setZero();

        for (int i = 0; i < rendering_particles_frame.size(); ++i) {
            // Row and column indices (zero-indexed)
            auto col = int(rendering_particles_frame[i][2]); // Sensor number
            auto row = int(rendering_particles_frame[i][3]); // Factory number

            // Factories themselves will have their second/fourth entries as -1 since they're never captured; skip them
            if (rendering_particles_frame[i][1] != -1 && row != -1) {
                // Travel time of this particle in hours
                allTimes(row, col) += (rendering_particles_frame[i][1] - rendering_particles_frame[i][0]) / framesPerHour;
                numDetected(row, col)++;

                // Delete captured particles
                rendering_particles.erase(rendering_particles.begin() + i);
                rendering_particles_vel.erase(rendering_particles_vel.begin() + i);
                rendering_particles_col.erase(rendering_particles_col.begin() + i);
                rendering_particles_frame.erase(rendering_particles_frame.begin() + i);
//                rendering_particles_rad.erase(rendering_particles_rad.begin() + i);r
            }
        }

        // Compute average travel time between each factory-sensor pair
        for (int i = 0; i < allTimes.rows(); ++i) {
            for (int j = 0; j < allTimes.cols(); ++j) {
                if (numDetected(i, j) > 0) {
                    allTimes(i, j) /= numDetected(i, j);
                }
                else {
                    // A travel time of -1 indicates no particles reached this factory from the current sensor
                    allTimes(i, j) = -1;
                }
            }
        }

        // Write to file
        outFile << "\t"
                << "Sensor 1 \t Sensor 2 \t Sensor 3 \t Sensor 4 \t Sensor 5 \t Sensor 6 \t Sensor 7 \t Sensor 8 \t Sensor 9"
                << std::endl;
        for (int i = 0; i < factoryPos.size(); ++i) {
            outFile << factoryName[i];
            for (int j = 0; j < sensorPos.size(); ++j) {
                if (allTimes(i, j) == -1) {
                    outFile << "\t" << "N/A";
                }
                else {
                    outFile << "\t" << std::to_string(allTimes(i, j));
                }
            }
            outFile << std::endl;
        }
        outFile.close();
        std::cout << "Calculated travel times, in hours, have been written to ../records/calcTravelTimes.txt as TSV; note that " <<
                  "these are based on rendered particles only" << std::endl;
    }
}
