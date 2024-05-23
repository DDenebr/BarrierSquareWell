#ifndef INITIALIZATION_INCLUDED
#define INITIALIZATION_INCLUDED

#include "simulation_system.h"

System InitializeRandomHardSphereECMC(const Parameters& env);

System InitializeFromDump(const string& dump_filename, const SquareWellCoefficient& coef);

#endif