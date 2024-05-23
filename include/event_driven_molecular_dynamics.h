#ifndef EVENT_DRIVEN_MOLECULAR_DYNAMICS_INCLUDED
#define EVENT_DRIVEN_MOLECULAR_DYNAMICS_INCLUDED

#include <list>
#include <functional>

#include "simulation_system.h"

using std::list;
using std::function;

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
    class TerminateEvent;
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
    void SetSamplingParameters(const string& sample_mode, const double& sample_interval, const vector<string>& sample_properties, const path& dump_directory);
    void SetResetParameters(const double& time_per_epoch);
    void SetResetEpoch(const double& total_time);
    void SetSamplingOperation(function<void(const double&, EDMD&)> sample_operation);
    void SetTerminateOperation(function<void(const double&, EDMD&)> terminate_operation);
    int GetEpoch();
    double GetTimePerEpoch();
    double GetSamplingTime();


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
    void InitializeTerminateNode(const double& terminate_time);
    void InitializeSampleNode(const double& sampling_time);

    //Event erase
    void EraseCrossNode(const PtrParticleEDMD& particle_cross);
    void EraseCollideNode(const PtrParticleEDMD& particle_collide);
    void EraseSquareWellNode(const PtrParticleEDMD& particle_squarewell);
    void EraseAndersenNode(const Node& node_andersen);
    void EraseResetNode(const Node& node_reset);
    void EraseTerminateNode(const Node& node_terminate);
    void EraseSampleNode(const Node& node_sample);

    //Event execution
    void ExecuteCrossNode(const Node& node);
    void ExecuteCollideNode(const Node& node);
    void ExecuteSquareWellNode(const Node& node);
    void ExecuteAndersenNode(const Node& node);
    void ExecuteResetNode(const Node& node);
    void ExecuteTerminateNode(const Node& node);
    void ExecuteSampleNode(const Node& node);

    //Tree execution
    void ExecuteEventTree();

    //EDMD
    void ExecuteEDMD(const double& termination_time, const double& bath_temperature, const double& thermostat_frequency, const double& initial_sample_time);

    //Various properties
    double Pressure(const double& t);
    double MaximumBubbleVolume(const double& t);

    protected:
    vector<PtrParticleEDMD> particleEDMD_pool;
    Tree event_tree;
    list<pair<double, double>> virial;

    double diameter_max;
    int collision_search_range;
    int squarewell_search_range;

    //Reset
    int epoch = 0;
    double TIME_PER_EPOCH_ = 100;
    int max_epoch = 0;
    double time_remainder = 0;

    //Sampling
    int sampling_time = 0;
    double sampling_interval = 0.1;
    string sampling_mode = "Arithmetic";
    vector<string> sampling_property_list{"dump"};
    path dump_directory;

    function<void(const double&, EDMD&)> sampling_operation = [&](const double& t, EDMD& edmd){

        const auto& properties = sampling_property_list;
        map<string, double> property_map;

        if (std::find(properties.begin(), properties.end(), "dump") != properties.end()){
            assert(std::filesystem::exists(dump_directory));
            const int& epoch = GetEpoch();
            const double& time_per_epoch = GetTimePerEpoch();
            std::stringstream dump_filename;
            dump_filename << "dump_" << sampling_time << ".dat";
            // if (sample_event->kinetic_temperature != 0)
            edmd.Dump(dump_directory, dump_filename.str(), epoch, time_per_epoch, t);
            property_map.emplace("time", t + epoch * time_per_epoch);
        }

        if (std::find(properties.begin(), properties.end(), "time") != properties.end()){
            property_map.emplace("time", t + epoch * TIME_PER_EPOCH_);
        }

        if (std::find(properties.begin(), properties.end(), "temperature") != properties.end()){
            const double& kinetic_temperature = edmd.KineticTemperature(t);
            property_map.emplace("temperature", kinetic_temperature);
        }

        if (std::find(properties.begin(), properties.end(), "pressure") != properties.end()){
            const double& pressure = edmd.Pressure(t);
            property_map.emplace("pressure", pressure);
        }

        if (std::find(properties.begin(), properties.end(), "bubble") != properties.end()){
            const double& max_bubble_volume = edmd.MaximumBubbleVolume(t);
            property_map.emplace("bubble", max_bubble_volume);
        }

        for (const auto& [name, value] : property_map)
            cout << std::setprecision(4) << value << ' ';
        cout << endl;
    };
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