#ifndef PARAMETERS_INCLUDED
#define PARAMETERS_INCLUDED

#include <chrono>
#include <cstring>
#include <execution>
#include <map>
#include <numeric>
#include <random>

#include "vector_operation.h"

using std::map;
using std::pair;
using std::string;

//Unit
class Unit
{
    public:
    Unit() :
        mass(1.0), sigma(1.0), epsilon(1.0), kB(1.0) {
    };
    //Swap (Copy & Swap Idiom)
    void swap(Unit& unit) {
        std::swap (mass, unit.mass);
        std::swap (sigma, unit.sigma);
        std::swap (epsilon, unit.epsilon);
        std::swap (kB, unit.kB);
    };

    protected:
    double mass;
    double sigma;
    double epsilon;
    double kB;
};

//Periodical Square Box
class PeriodicalSquareBox
{
    public:
    PeriodicalSquareBox(const double& box_length, const int& cell_number_per_edge){
        L = box_length;
        n = cell_number_per_edge; 
        l = L / static_cast<double>(n); 
    };
    //Copy & Swap Idiom
    void swap(PeriodicalSquareBox& box) {
        std::swap (n, box.n);
        std::swap (L, box.L);
        std::swap (l, box.l);
    };
    
    protected:
    int n;                 //cell number per edge
    double L;              //box length
    double l;              //cell length
};

//Particle parameters and number
struct ParticleData
{
    ParticleData() = default;
    ParticleData(const string& type, const double& mass, const double& diameter, const unsigned& number) :
        type(type), mass(mass), diameter(diameter), number(number){
    };
    //Copy & Swap Idiom
    void swap(ParticleData& pdata){
        std::swap (type, pdata.type);
        std::swap (mass, pdata.mass);
        std::swap (diameter, pdata.diameter);
        std::swap (number, pdata.number);
    }

    string type;           //particle type identifier
    double mass;           //particle type mass 
    double diameter;       //particle type diameter
    unsigned number;       //particle type number in the system
};

class SquareWellCoefficient
{
    public:
    SquareWellCoefficient(const double& width, const double& depth, const double& barrier) :
        squarewell_width(width), squarewell_depth(depth), squarewell_barrier(barrier){
    };
    //Copy & Swap Idiom
    void swap(SquareWellCoefficient& coef){
        std::swap (squarewell_width, coef.squarewell_width);
        std::swap (squarewell_depth, coef.squarewell_depth);
        std::swap (squarewell_barrier, coef.squarewell_barrier);
    }
    
    protected:             
    double squarewell_width;          //attraction range of the square well
    double squarewell_depth;          //attraction energy of the square well
    double squarewell_barrier;        //barrier of the square well
};

class Parameters : public Unit, public PeriodicalSquareBox, public SquareWellCoefficient
{
    public:
    //Pseudorandom number generator
    std::mt19937 gen;                                           //Mersenne Twister pseudo-random generator of 32-bit numbers
    unsigned int seed;                                          //Seed of generator

    //Friend classes
    friend class Cell;
    friend class Particle;

    //Constructor
    Parameters(const double& box_length, const int& cell_number_per_edge, const double& width, const double& depth, const double& barrier) : 
        Unit(), 
        PeriodicalSquareBox(box_length, cell_number_per_edge), 
        SquareWellCoefficient(width, depth, barrier) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        gen.seed(seed);
    };
    Parameters(const PeriodicalSquareBox& simulation_box, const SquareWellCoefficient& square_well_coefficient) :
        Unit(), 
        PeriodicalSquareBox(simulation_box), 
        SquareWellCoefficient(square_well_coefficient) {
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        gen.seed(seed);
    };
    //Copy constructor
    Parameters(const Parameters&) = default;
    //Move constructor
    Parameters(Parameters&&) = default;
    //Copy & Swap Idiom
    void swap(Parameters& env) {
        env.Unit::swap(*this);
        env.PeriodicalSquareBox::swap(*this);
        env.SquareWellCoefficient::swap(*this);
        std::swap(gen, env.gen);
        std::swap(seed, env.seed);
        std::swap(particle_data, env.particle_data);
    }

    //Add particle
    void CreateAtom(const string& type, const double& mass, const double& diameter, const unsigned& number){
        particle_data.emplace(type, ParticleData(type, mass, diameter, number));
    };
    void CreateAtom(const ParticleData& pdata){
        particle_data.emplace(pdata.type, pdata);
    };

    //Total particle number of the system
    unsigned Ntot() const{
        return std::accumulate(particle_data.begin(), particle_data.end(), 0, 
            [](const int& n_sum, const std::pair<string, ParticleData>& pdata)
            { return n_sum + pdata.second.number;}); 
    };

    //Packing fraction of the system
    double eta() const{
        return std::accumulate(particle_data.begin(), particle_data.end(), 0.0, 
            [](const double& v_sum, const std::pair<string, ParticleData>& pdata)
            {return v_sum + pdata.second.diameter * pdata.second.diameter * pdata.second.diameter * pdata.second.number;})
            * M_PI / (L * L * L * 6.0);        
    }

    //Max particle diameter
    double dmax() const{
        return particle_data.size() == 0 ? 0 : 
            std::max_element(particle_data.begin(), particle_data.end(), 
            [](const pair<string, ParticleData>& largest, const pair<string, ParticleData>& pdata){
                return largest.second.diameter < pdata.second.diameter;
            })->second.diameter;
    }

    protected:
    //Particle information
    map<string, ParticleData> particle_data;
};

#endif