#include <cassert>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <fstream>

#include <omp.h>

#include "event_driven_molecular_dynamics.h"
#include "initialization.h"

using std::cout;
using std::endl;
using std::filesystem::path;

int main (int argc, char* argv[]){

    auto log_output = std::freopen("/home/raop0002/BarrierSquareWell/log/output.o", "w", stdout);
    auto log_error = std::freopen("/home/raop0002/BarrierSquareWell/log/error.e", "w", stderr);

    omp_set_num_threads(8);

    Unit unit;
    PeriodicalSquareBox box(47, 31);
    SquareWellCoefficient coef(1.5, -1, 0.5);
    ParticleData pdata("0", 1, 1, 20000);

    Parameters env(box, coef);
    env.CreateAtom(pdata);

    System sys(InitializeRandomHardSphereECMC(env));

    EDMD edmd(sys);

    path dump_path("/home/raop0002/BarrierSquareWell/dump/T0.5_phi0.1");

    edmd.SetSamplingParameters(10 * (std::pow(10, 0.1) - 1), dump_path);
    edmd.InitializeParticleVelocity(0.5);
    edmd.ExecuteEDMD(1e5);

    return 0;
}