#ifndef EVENT_DRIVEN_MOLECULAR_DYNAMICS_INCLUDED
#define EVENT_DRIVEN_MOLECULAR_DYNAMICS_INCLUDED

#include <list>

#include "simulation_system.h"

using std::list;

#define MAX_PARTICLE_COLLIDE_NODE_LIST_SIZE 64
#define MAX_PARTICLE_SQUAREWELL_NODE_LIST_SIZE 64
#define MAX_NEIGHBOR_SIZE 32
#define MAX_VIRIAL_SIZE 65536

class EDMD : public System
{
    public:

    //Friend classes

    //Nested classes
    //Events 
    class Event;
    class CrossEvent;
    class CollideEvent;
    class ResetEvent;
    class SampleEvent;
    class AndersenThermostatEvent;

    //Specialized to square well
    class SquareWellEvent;

    //Particle
    class ParticleEDMD;

    //Typedefs
    typedef shared_ptr<Event> PtrEvent;
    typedef shared_ptr<ParticleEDMD> PtrParticleEDMD;
    typedef map<double, PtrEvent> Tree;
    typedef typename Tree::iterator Node; 

    //Default constructor deleted
    EDMD() = delete;
    //Copy constructor
    EDMD(const EDMD&) = default;
    //Move constructor
    EDMD(EDMD&&) = default;
    //Copy constructor from System
    EDMD(const System& system);
    //Destructor
    virtual ~EDMD() = default;

    //Initialize kinetic temperature data output filestream
    void SetSamplingParameters(const double& sample_interval, const path& dump_directory);

    //Event Functions
    //Add event in particle event list
    void AddCrossInParticle(const Node& node_cross);
    void AddCollideInParticle(const Node& node_collide);
    void AddSquareWellInParticle(const Node& node_squarewell);

    //Remove event in particle event list
    void RemoveCrossFromParticle(const Node& node_cross);
    void RemoveCollideFromParticle(const Node& node_collide);
    void RemoveSquareWellFromParticle(const Node& node_squarewell);

    //Event initialize
    void InitializeCrossNode(PtrParticleEDMD& particle_cross);
    void InitializeCollideNode(PtrParticleEDMD& particle_collide1, PtrParticleEDMD& particle_collide2);
    void InitializeSquareWellNode(PtrParticleEDMD& particle_squarewell1, PtrParticleEDMD& particle_squarewell2);
    void InitializeAndersenNode(const double& last_thermostat_time, const double& T_bath, const double& lambda);
    void InitializeResetNode();
    void InitializeSampleNode(const double& sampling_time);

    //Event erase
    void EraseCrossNode(const PtrParticleEDMD& particle_cross);
    void EraseCollideNode(const PtrParticleEDMD& particle_collide);
    void EraseSquareWellNode(const PtrParticleEDMD& particle_squarewell);
    void EraseAndersenNode(const Node& node_andersen);
    void EraseResetNode(const Node& node_reset);
    void EraseSampleNode(const Node& node_sample);

    //Event execution
    void ExecuteCrossNode(const Node& node);
    void ExecuteCollideNode(const Node& node);
    void ExecuteSquareWellNode(const Node& node);
    void ExecuteAndersenNode(const Node& node);
    void ExecuteResetNode(const Node& node);
    void ExecuteSampleNode(const Node& node);

    //Tree execution
    void ExecuteEventTree();

    //EDMD
    void ExecuteEDMD(const double& termination_time);

    protected:
    vector<PtrParticleEDMD> particleEDMD_pool;
    Tree event_tree;
    list<pair<double, double>> virial;

 
    double diameter_max;
    int collision_search_range;
    int squarewell_search_range;
};

class EDMD::ParticleEDMD : public Particle
{
    public:

    //Friend classes
    friend class EDMD;
    friend class CrossEvent;
    
    //Default constructor
    ParticleEDMD() :
        Particle(),
        node_cross(NULL_TREE.end()){
        node_collide_list.reserve(MAX_PARTICLE_COLLIDE_NODE_LIST_SIZE);
        node_squarewell_list.reserve(MAX_PARTICLE_SQUAREWELL_NODE_LIST_SIZE);
        neighbor_list.reserve(MAX_NEIGHBOR_SIZE);
        COLLIDE_LIST_ANCHOR_ = &node_collide_list.front();
        SQUAREWELL_LIST_ANCHOR_ = &node_squarewell_list.front();
    };
    //Copy constructor
    ParticleEDMD(const ParticleEDMD& particle_edmd) :
        Particle(particle_edmd),
        node_cross(particle_edmd.node_cross),
        node_collide_list(particle_edmd.node_collide_list),
        node_squarewell_list(particle_edmd.node_squarewell_list),
        neighbor_list(particle_edmd.neighbor_list){
        node_collide_list.reserve(MAX_PARTICLE_COLLIDE_NODE_LIST_SIZE);
        node_squarewell_list.reserve(MAX_PARTICLE_SQUAREWELL_NODE_LIST_SIZE);
        neighbor_list.reserve(MAX_NEIGHBOR_SIZE);
        COLLIDE_LIST_ANCHOR_ = &node_collide_list.front();
        SQUAREWELL_LIST_ANCHOR_ = &node_squarewell_list.front();
    }
    //Move constructor
    ParticleEDMD(ParticleEDMD&& particle_edmd) :
        Particle(std::move(particle_edmd)),
        node_cross(std::move(particle_edmd.node_cross)),
        node_collide_list(std::move(particle_edmd.node_collide_list)),
        node_squarewell_list(std::move(particle_edmd.node_squarewell_list)),
        neighbor_list(std::move(particle_edmd.neighbor_list)),
        COLLIDE_LIST_ANCHOR_(std::exchange(particle_edmd.COLLIDE_LIST_ANCHOR_, nullptr)),
        SQUAREWELL_LIST_ANCHOR_(std::exchange(particle_edmd.SQUAREWELL_LIST_ANCHOR_, nullptr))
        {
    };
    //Destructor
    virtual ~ParticleEDMD() = default;
    //Copy construct from System::Particle
    ParticleEDMD(const Particle& particle) :
        Particle(particle),
        node_cross(NULL_TREE.end()){
        node_collide_list.reserve(MAX_PARTICLE_COLLIDE_NODE_LIST_SIZE);
        node_squarewell_list.reserve(MAX_PARTICLE_SQUAREWELL_NODE_LIST_SIZE);
        neighbor_list.reserve(MAX_NEIGHBOR_SIZE);
        COLLIDE_LIST_ANCHOR_ = &node_collide_list.front();
        SQUAREWELL_LIST_ANCHOR_ = &node_squarewell_list.front();
    };
    //Virtual copy constructor
    virtual System::PtrParticle clone() const override{
        return System::PtrParticle(new ParticleEDMD(*this));
    };
    //Swap (copy and swap idiom)
    void swap(ParticleEDMD& particle_edmd){
        static_cast<Particle&>(particle_edmd).swap(static_cast<Particle&>(*this));
        std::swap(node_cross, particle_edmd.node_cross);
        std::swap(node_collide_list, particle_edmd.node_collide_list);
        std::swap(node_squarewell_list, particle_edmd.node_squarewell_list);
        std::swap(neighbor_list, particle_edmd.neighbor_list);
        std::swap(COLLIDE_LIST_ANCHOR_, particle_edmd.COLLIDE_LIST_ANCHOR_);
        std::swap(SQUAREWELL_LIST_ANCHOR_, particle_edmd.SQUAREWELL_LIST_ANCHOR_);
    };
    //Unify copy assignment
    ParticleEDMD& operator=(ParticleEDMD particle_edmd){
        particle_edmd.swap(*this);
        return *this;
    };
    virtual void clear(){
        ParticleEDMD().swap(*this);
    };

    protected:
    // Particle* base;

    //Event data
    //Cross event
    Node node_cross;
    //Collide
    vector<Node> node_collide_list;                                         //All possible collide event involved
    //Squarewell
    vector<Node> node_squarewell_list;                                      //All possible square-well event involved
    //Bond
    vector<PtrParticleEDMD*> neighbor_list;                                 //List of bonded neighbors
    //Add collision index
    Node* AppendCollideNode(const Node& new_node_collide){
        return append(new_node_collide, node_collide_list, COLLIDE_LIST_ANCHOR_);
    };
    //Erase collision index            
    Node* EraseCollideNode(Node* const& ptr_old_node_collide){
        return erase(ptr_old_node_collide, node_collide_list);
    };
    //Add squarewell index
    Node* AppendSquareWellNode(const Node& new_node_squarewell){
        return append(new_node_squarewell, node_squarewell_list, SQUAREWELL_LIST_ANCHOR_);
    };
    //Erase squarewell index            
    Node* EraseSquareWellNode(Node* const& ptr_old_node_squarewell){
        return erase(ptr_old_node_squarewell, node_squarewell_list);
    };         

    private:
    static Tree NULL_TREE;
    Node* COLLIDE_LIST_ANCHOR_;
    Node* SQUAREWELL_LIST_ANCHOR_;
};

#endif