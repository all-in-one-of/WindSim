#ifndef MACGrid_H_
#define MACGrid_H_

#pragma warning(disable: 4244 4267 4996)

#include "open_gl_headers.h" 
#include "vec.h"
#include "grid_data.h"
#include "grid_data_matrix.h" 
#include <Partio.h>
class Camera;

class MACGrid
{

public:
	MACGrid();
	~MACGrid();
	MACGrid(const MACGrid& orig);
	MACGrid& operator=(const MACGrid& orig);

	void reset();

	void advectVelocity(double dt);
	void addExternalForces(double dt); // Not used in wind sim; useful in 3D simulations
	void project(double dt);
	void advectTemperature(double dt);
	void advectDensity(double dt);
	void advectRenderingParticles(double dt);
	void advectConcentration(double dt);

    // Public function that applies mass diffusion
    void diffuse(double dt);

protected:

	// Setup
	void initialize();

	// Simulation
	void computeBuoyancy(double dt);
	void computeVorticityConfinement(double dt);

	// GridData accessors
	enum Direction { X, Y, Z };
	vec3 getVelocity(const vec3& pt);
	double getVelocityX(const vec3& pt);
	double getVelocityY(const vec3& pt);
	double getVelocityZ(const vec3& pt);
	double getTemperature(const vec3& pt);
	double getDensity(const vec3& pt);
	double getMethylConcentration(const vec3& pt);
	double getChloroConcentration(const vec3& pt);
	double getAGOCConcentration(const vec3& pt);
	double getAppConcentration(const vec3& pt);
	vec3 getCenter(int i, int j, int k);

	
	vec3 getRewoundPosition(const vec3 & currentPosition, const double dt);
	vec3 clipToGrid(const vec3& outsidePoint, const vec3& insidePoint);
	double getSize(int dimension);
	bool isValidCell(int i, int j, int k);
	bool isValidFace(int dimension, int i, int j, int k);
	vec3 getFacePosition(int dimension, int i, int j, int k);
	void calculateAMatrix();
	bool preconditionedConjugateGradient(const GridDataMatrix & A, GridData & p, const GridData & d, int maxIterations, double tolerance);
	void calculatePreconditioner(const GridDataMatrix & A);
	void applyPreconditioner(const GridData & r, const GridDataMatrix & A, GridData & z);
	double dotProduct(const GridData & vector1, const GridData & vector2);
	void add(const GridData & vector1, const GridData & vector2, GridData & result);
	void subtract(const GridData & vector1, const GridData & vector2, GridData & result);
	void multiply(const double scalar, const GridData & vector, GridData & result);
	double maxMagnitude(const GridData & vector);
	void apply(const GridDataMatrix & matrix, const GridData & vector, GridData & result);


	GridDataX mU; // X component of velocity, stored on X faces, size is (dimX+1)*dimY*dimZ
	GridDataY mV; // Y component of velocity, stored on Y faces, size is dimX*(dimY+1)*dimZ
	GridDataZ mW; // Z component of velocity, stored on Z faces, size is dimX*dimY*(dimZ+1)
	GridData mP;  // Pressure, stored at grid centers, size is dimX*dimY*dimZ
	GridData mD;  // Density, stored at grid centers, size is dimX*dimY*dimZ
	GridData mT;  // Temperature, stored at grid centers, size is dimX*dimY*dimZ

	// Amount of each species produced (stored at cell centers)
	GridData nMethyl;
	GridData nChloro;
	GridData nAGOC;
	GridData nApp;

	// Species concentrations (stored at cell centers)
	GridData C_Methyl;
	GridData C_Chloro;
	GridData C_AGOC;
	GridData C_App;

	// Partial derivatives of concentrations
	GridData firstPartialX;
	GridData firstPartialY;
	GridData secondPartialXX;
	GridData secondPartialYY;

	// Matrices	
	GridDataMatrix AMatrix;
	GridData precon;

	// Data from VAST sheets -- sensors listed in numerical order
	const std::vector<vec3> sensorPos = {vec3(62, 21, 0), vec3(66, 35, 0), vec3(76, 41, 0), vec3(88, 45, 0),
                                         vec3(103, 43, 0), vec3(102, 22, 0), vec3(89, 3, 0), vec3(74, 7, 0), vec3(119, 42, 0)};
	const std::vector<vec3> factoryPos = {vec3(89, 27, 0), vec3(90, 21, 0), vec3(109, 26, 0), vec3(120, 22, 0)};
    const std::vector<std::string> factoryName = {"RFE", "KOF", "RCT", "ISB"};

public:
    enum Chemical {Methyl, Chloro, AGOC, App};

	// Rendering particles
	std::vector<vec3> rendering_particles;
	std::vector<vec3> rendering_particles_vel;
	std::vector<vec3> rendering_particles_col;
    std::vector<int> rendering_particles_rad;
    std::vector<vec4> rendering_particles_frame; // Info stored as (frameReleased, frameCaptured, sensor, factory)

	void saveParticle(std::string filename, bool save);
    void saveVelocityField(std::string filename);

    void writeAmountProduced();
    void writeTravelTimes(bool calculateActual);

    // Main particle update call
    void updateSources(const vec3 &wind, const std::vector<vec4> &conc, int frameNum);

private:
	// Particle update helpers
    void refresh(int i, int j, int k, int chemical, int frameNum, int sensorIdx);
	void compoundSource(const vec3 &wind, const vec4 &conc, int i, int j, int frameNum, int sensorIdx);

    // Computes vorticity vector (not used in wind sim)
	vec3 computeOmega(int i, int j, int k);

	// Partial derivative helpers
	void computeFirstPartial(const GridData& C, int dir);
	void computeSecondPartial(const GridData& C, int dir);

    // Apply Conservation of Species
    void applyMassDiffusion(double dt, const GridData& C, double D, GridData& n);

    // Distance between a factory and sensor
    double calculateDistance(int factoryIdx, int sensorIdx);

	// Check for a particle in a factory's surrounding square
    int reachesFactory(vec3 particle);

};

#endif