#include <cassert>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <fstream>

#include <omp.h>

#include "event_driven_molecular_dynamics.h"

using std::cout;
using std::endl;
using std::filesystem::path;

int main (int argc, char* argv[]){

    Unit unit;
    PeriodicalSquareBox box(16, 10);
    SquareWellCoefficient coef(1.5, -1, 0.5);
    ParticleData pdata("a", 1, 1, 1);

    Parameters env(box, coef);
    env.CreateAtom(pdata);
    env.CreateAtom("b", 1, 1, 1);

    System sys(env);

    auto p1 = std::make_shared<System::Particle>("a", 1, 1, 0);
    auto p2 = std::make_shared<System::Particle>("b", 1, 1, 0.1);

    p1->location(array<int, 3>{0, 0, 0}, array<double, 3>{0, 0, 0});
    p2->location(array<int, 3>{1, 2, 3}, array<double, 3>{0.1, 0.1, 0.1});

    p2->velocity(array<double, 3>{-1, -2, -3});

    sys.AddParticle(p1);
    sys.AddParticle(p2);

    double dist = sys.distance(*p1, *p2, 0);
    auto time = sys.contacttime(*p1, *p2, 1);

    EDMD edmd(sys);

    edmd.ExecuteEDMD(10);

    return 0;
}