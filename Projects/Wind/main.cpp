#include "smoke_sim.h"
#include "camera.h"
#include "fps.h"
#include "constants.h"
#include <stdio.h>
#include <cmath>
#include "open_gl_headers.h"
#include "basic_math.h"
#include <string.h>
#include <iostream>

// Geometry and whatnot
SmokeSim theSmokeSim;

// Simple infinite loop for running the simulation (has internal checks to stop itself when appropriate)
int main(int argc, char **argv) {
    while (true) {
        theSmokeSim.step();
    }
}
