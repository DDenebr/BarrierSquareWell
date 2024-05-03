#ifndef SIMULATION_SYSTEM_INCLUDED
#define SIMULATION_SYSTEM_INCLUDED

#include <filesystem>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "parameters.h"
#include "utilities.h"
#include "vector_operation.h"

using std::cout, std::cerr, std::endl;
using std::filesystem::path;
using std::shared_ptr;
using std::make_shared;
using std::vector;

#define MAX_CELL_PARTICLE_LIST_SIZE 64
#define MAX_SYSTEM_PARTICLE_POOL_SIZE 131072

class System : public Parameters
{
    public:

    //Friend class
    friend System InitializeRandomHardSphereECMC(const Parameters& env);

    //Nested struct & template
    struct Cell;
    struct Particle;

    //Typedefs
    typedef shared_ptr<Particle> PtrParticle;

    System() = delete;
    //Default constructor, system without particle
    System(const Parameters& env);
    //Copy constructor
    System(const System& system);
    //Move constructor
    System(System&& system);
    //Swap (Copy & Swap Idiom)
    void swap(System& system);
    //Unify copy assignment
    System& operator=(System system);
    //Clear particle pool and cell grid
    void clear();

    //Add particle in cell particle list
    void AddParticleInCell(PtrParticle& particle_add);
    //Remove particle from cell particle list
    void RemoveParticleFromCell(PtrParticle& particle_remove);

    //Calculate the distance vector bewteen two particles
    array<double, 3> distancevector(const Particle& p1, const Particle& p2, const double& t); 
    //Calculate the distance between two particles
    double distance(const Particle& p1, const Particle& p2, const double& t);   
    //Calculate the contact time between two particles                                             
    array<double, 2> contacttime(const Particle& p1, const Particle& p2, const double& len);
    
    //Swap two particle in the particle list
    void SwapParticle(PtrParticle& p1, PtrParticle& p2);
    //Random add a particle to the system, return true when successfully added
    bool AddParticle(const PtrParticle& particle_add);
    //Move a particle in the system
    void MoveParticle(PtrParticle& particle_move, const array<int, 3>& new_cell, const array<double, 3>& new_coordinate);
    //Remove a particle from the system
    void RemoveParticle(PtrParticle& particle_remove);

    //Initialize particle velocity according to a initial kinetic temperature
    void InitializeParticleVelocity(const double& kinetic_temperature);

    //Export system configuration to dat file, LAMMPS style
    void Dump(const path& dump_directory, const string& dump_filename, const unsigned& epoch, const double& current_time);

    //TEST
    bool CheckOverlap();

    protected:
    vector<PtrParticle> particle_pool;
    vector<vector<vector<Cell>>> cell_grid;             //cell list
    // vector<Cluster> cluster_list;

    private:
    PtrParticle const* PARTICLE_POOL_ANCHOR_;
};

struct System::Cell
{
    //Default constructor
    Cell(){
        particle_list.reserve(MAX_CELL_PARTICLE_LIST_SIZE);
        PARTICLE_LIST_ANCHOR_ = &particle_list.front();
    };
    //Copy constructor
    Cell(const Cell& cell) : 
        particle_list(cell.particle_list){
        particle_list.reserve(MAX_CELL_PARTICLE_LIST_SIZE);
        PARTICLE_LIST_ANCHOR_ = &particle_list.front();
    };
    //Move constructor
    Cell(Cell&& cell) : 
        particle_list(std::move(cell.particle_list)),
        PARTICLE_LIST_ANCHOR_(std::exchange(cell.PARTICLE_LIST_ANCHOR_, nullptr)){
    };
    //Swap (Copy & Swap Idiom)
    void swap(Cell& cell){
       std::swap(particle_list, cell.particle_list);
       std::swap(PARTICLE_LIST_ANCHOR_, cell.PARTICLE_LIST_ANCHOR_);
    };
    //Unifying copy assignment
    Cell& operator=(Cell cell){
        cell.swap(*this);
        return *this;
    }

    //Number of particles in the cell
    int N(){
        return particle_list.size();
    };
    //Add a particle, return a pointer to the added position in p                                                
    int* append(const int& new_particle){
        return ::append(new_particle, particle_list, PARTICLE_LIST_ANCHOR_);
    };
    //Erase a particle, return a pointer to the erased position in p                      
    int* erase(int* const& ptr_old_particle){
        return ::erase(ptr_old_particle, particle_list);
    };       
    //Clear the cell               
    void clear(){
        particle_list.clear();
    }; 

    vector<int> particle_list;                          //cell particle list                                    

    private:
    int* PARTICLE_LIST_ANCHOR_;                         //cell particle list address (avoid reallocation)
};


struct System::Particle
{
    //Default Constructor
    Particle(const string& type = "", const double& diameter = 0, const double& mass = 0, const double& time = 0) : 
        type(type), 
        d(diameter), m(mass), t(time),
        b{}, c{}, v{}, r{},
        ptr_cell_particle_list(nullptr) {
    };
    //Copy constructor
    Particle(const Particle&) = default;
    //Move constructor
    Particle(Particle&&) = default;
    //Virtual copy constructor
    virtual PtrParticle clone() const{
        return std::make_shared<Particle>(Particle(*this));
    };
    //Swap (Copy & Swap Idiom)
    void swap(Particle& particle){
       std::swap(type, particle.type); 
       std::swap(d, particle.d); std::swap(m, particle.m); std::swap(t, particle.t);
       std::swap(b, particle.b); std::swap(c, particle.c); std::swap(r, particle.r); std::swap(v, particle.v);
       std::swap(ptr_cell_particle_list, particle.ptr_cell_particle_list);
    };
    //Unifying copy assignment
    Particle& operator=(Particle particle){
        particle.swap(*this);
        return *this;
    }

    void location(const array<int, 3>& cell, const array<double, 3>& coordinate){
        c = cell; r = coordinate;
    };
    void velocity(const array<double, 3>& velocity){
        v = velocity;
    }
    //Clear the particle
    virtual void clear(){
       Particle().swap(*this);
    };    

    //Basic information
    string type;                                        //particle type
    double d;                                           //particle diameter
    double t;                                           //particle time updated
    double m;                                           //particle mass

    array<int, 3> b;                                    //particle box coordinate (default in (0,0,0) box, record of box crossing)
    array<int, 3> c;                                    //particle cell coordinate in the box [0, m]
    array<double, 3> r;                                 //particle coordinate in the cell, [0, 1]
    array<double, 3> v;                                 //particle velocity

    //links to other class
    int* ptr_cell_particle_list;                        //pointer to the particle information stored in the cell particle list
};

inline System::System(const Parameters& env) :
    Parameters(env),
    cell_grid(n, vector<vector<Cell>>(n, vector<Cell>(n, Cell()))) {
    particle_pool.reserve(MAX_SYSTEM_PARTICLE_POOL_SIZE);
    PARTICLE_POOL_ANCHOR_ = &particle_pool.front();
};

inline System::System(const System& system) : 
    Parameters(system),
    cell_grid(system.cell_grid) {
    particle_pool.reserve(MAX_SYSTEM_PARTICLE_POOL_SIZE);
    PARTICLE_POOL_ANCHOR_ = &particle_pool.front();
    for (const auto& pi : system.particle_pool)
        particle_pool.emplace_back(pi->clone());
    //Adjust ptr_cell_particle_list (pointer to particle list element in cell)
    for (auto& ci : cell_grid)
    for (auto& cj : ci)
    for (auto& ck : cj)
    for (int& i : ck.particle_list)
        particle_pool.at(i)->ptr_cell_particle_list = &i;
};

inline System::System(System&& system) : 
    Parameters(std::move(system)),
    particle_pool(std::move(system.particle_pool)),
    cell_grid(std::move(system.cell_grid)),
    PARTICLE_POOL_ANCHOR_(std::exchange(system.PARTICLE_POOL_ANCHOR_, nullptr)){
};

inline void System::swap(System& system){
    system.Parameters::swap(*this);
    std::swap(particle_pool, system.particle_pool);
    std::swap(cell_grid, system.cell_grid);
    std::swap(PARTICLE_POOL_ANCHOR_, system.PARTICLE_POOL_ANCHOR_);
};

inline System& System::operator=(System system){
    system.swap(*this);
    return *this;
};

inline void System::clear(){
    particle_pool.clear();
    for (auto& ci : cell_grid)
    for (auto& cj : ci)
    for (auto& ck : cj)
        ck.clear();
};

inline void System::AddParticleInCell(PtrParticle& particle_add){  
    Cell& c = cell_grid[particle_add->c[0]][particle_add->c[1]][particle_add->c[2]];
    particle_add->ptr_cell_particle_list = c.append(index(&particle_add, particle_pool));  
};

inline void System::RemoveParticleFromCell(PtrParticle& particle_remove){
    Cell& c = cell_grid[particle_remove->c[0]][particle_remove->c[1]][particle_remove->c[2]];
    //Erase the old particle from the cell particle list.
    int* ptr_erased = c.erase(particle_remove->ptr_cell_particle_list);
    //Adjust the pointer of the last particle in the cell particle list
    if (ptr_erased != &*c.particle_list.end())
        particle_pool.at(*ptr_erased)->ptr_cell_particle_list = ptr_erased;
    //Set pointer to nullptr
    particle_remove->ptr_cell_particle_list = nullptr;
};

#endif