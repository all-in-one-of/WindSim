
#ifndef smokeSim_H_
#define smokeSim_H_

#include "mac_grid.h"
#include "vec.h"
#include <Partio.h>

class Camera;
class SmokeSim
{
public:
   SmokeSim();
   virtual ~SmokeSim();

   virtual void reset();
   virtual void step();

protected:
   virtual void grabScreen();

protected:
	MACGrid mGrid;
	int mFrameNum;
	int mTotalFrameNum;

	// All wind data
	std::vector<vec3> windVec;
	// This vec<vec<vec4>> represents: {{conc of all chemicals at sensor 1 at time=0, conc of all chemicals at sensor 2 at time=0, ...}, {conc of all at sensor 1 at time=1, conc of all at sensor 2 at time=1, ...}, ...}
	std::vector<std::vector<vec4>> concVec;

private:
	void fillConcVec(std::vector<std::vector<vec4>>& c);
	void fillWindVec(std::vector<vec3>& w);
};

#endif