#ifndef EVENT_INCLUDED
#define EVENT_INCLUDED

#include <fstream>

#include "event_driven_molecular_dynamics.h"

using std::ofstream;

class EDMD::Event 
{
    public:
    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const = 0;

    //Default constructor
    Event(const double& time = 0) : t(time) {};
    //Copy constructor
    Event(const Event&) = default;
    //Move constructor
    Event(Event&&) = default;
    //Virtual copy cosntructor
    virtual PtrEvent clone() const = 0;
    //Destructor
    virtual ~Event() {};
    //Clear
    virtual void clear() = 0;

    protected:
    double t;                                           //Event time                                       
};

class EDMD::CrossEvent : public Event
{
    public:

    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "cross";
    }

    //Default constructor
    CrossEvent() : 
        Event(), 
        ptr_particle_cross(nullptr), 
        cross_direction{}{
    };
    //Constructor
    CrossEvent(const double& time, PtrParticleEDMD& particle_cross, const array<int, 2>& cross_direction) : 
        Event(time),
        ptr_particle_cross(&particle_cross),
        cross_direction(cross_direction){  
    };
    //Copy construcotr
    CrossEvent(const CrossEvent&) = default;
    //Move constructor
    CrossEvent(CrossEvent&&) = default;
    //Destructor
    virtual ~CrossEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new CrossEvent(*this));
    };
    //Swap (copy and swap idiom)
    void swap(CrossEvent& cross_event){
        std::swap(t, cross_event.t);
        std::swap(ptr_particle_cross, cross_event.ptr_particle_cross);
        std::swap(cross_direction, cross_event.cross_direction);
    };
    //Unify copy assignment
    CrossEvent& operator= (CrossEvent cross_event){
        cross_event.swap(*this);
        return *this;
    };
    //Clear
    virtual void clear() override{
        CrossEvent().swap(*this);
    };

    protected:
    PtrParticleEDMD* ptr_particle_cross;                       //Particle to cross
    array<int, 2>cross_direction;                   //[0], xyz directions, [1], +- directions
};

class EDMD::CollideEvent : public Event
{
    public:
    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "collide";
    }

    //Default constructor
    CollideEvent() :
        Event(),
        ptr_particle_collide{nullptr, nullptr},
        ptr_node_collide{nullptr, nullptr}{
    };
    //Constructor
    CollideEvent(const double& time, PtrParticleEDMD& particle_collide1, PtrParticleEDMD& particle_collide2) :
        Event(time),
        ptr_particle_collide{&particle_collide1, &particle_collide2},
        ptr_node_collide{}{ 
    };
    //Copy constructor
    CollideEvent(const CollideEvent&) = default;
    //Move constructor
    CollideEvent(CollideEvent&&) = default;
    //Destructor
    virtual ~CollideEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new CollideEvent(*this));
    };
    //Swap (copy and swap idiom) 
    void swap(CollideEvent& collide_event){
        std::swap(t, collide_event.t);
        std::swap(ptr_particle_collide, collide_event.ptr_particle_collide);
        std::swap(ptr_node_collide, collide_event.ptr_node_collide);
    };
    //Unify copy assignment
    CollideEvent& operator= (CollideEvent collide_event){
        collide_event.swap(*this);
        return *this;
    };
    //Clear
    virtual void clear() override{
        CollideEvent().swap(*this);
    };

    protected:
        array<PtrParticleEDMD*, 2> ptr_particle_collide;
        array<Node*, 2> ptr_node_collide;
};

class EDMD::ResetEvent : public Event
{
    public:
    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "reset";
    }

    //Default constructor
    ResetEvent() : 
        Event(){
    };
    ResetEvent(const double& time_per_epoch) : 
        Event(time_per_epoch){
    };
    //Copy constructor
    ResetEvent(const ResetEvent&) = default;
    //Move constructor
    ResetEvent(ResetEvent&&) = default;
    //Destructor
    virtual ~ResetEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new ResetEvent(*this));
    };
    //Swap (copy and swap idiom)
    void swap(ResetEvent& rest_event){
        std::swap(t, rest_event.t);
    }
    //Unify copy assignment
    ResetEvent& operator=(ResetEvent rest_event){
        rest_event.swap(*this);
        return *this;
    }
    //Clear
    virtual void clear() override{
        ResetEvent().swap(*this);
    };

    protected:
};

class EDMD::TerminateEvent : public Event
{
    public:
    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "terminate";
    }

    //Default constructor
    TerminateEvent() = default; 
    //Constructor
    TerminateEvent(const double& time) :
        Event(time) {
    };
    //Copy constructor
    TerminateEvent(const TerminateEvent&) = default;
    //Move constructor
    TerminateEvent(TerminateEvent&&) = default;
    //Destructor
    virtual ~TerminateEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new TerminateEvent(*this));
    };
    //Swap (copy and swap idiom)
    void swap(TerminateEvent& terminate_event){
        std::swap(t, terminate_event.t);
    }
    //Unify copy assignment
    TerminateEvent& operator=(TerminateEvent terminate_event){
        terminate_event.swap(*this);
        return *this;
    }
    //Clear
    virtual void clear() override{
        TerminateEvent().swap(*this);
    };

    protected:
    inline static function<void(const double&, EDMD&)> terminate_operation;
};

class EDMD::SampleEvent : public Event
{
    public:
    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "sample";
    }

    //Default constructor
    SampleEvent() = default;
    //Constructor
    SampleEvent(const double& sampling_time) :
        Event(sampling_time){
        // kinetic_temperature(0.0){
    };
    //Copy constructor
    SampleEvent(const SampleEvent&) = default;
    //Move constructor
    SampleEvent(SampleEvent&&) = default;
    //Destructor
    virtual ~SampleEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new SampleEvent(*this));
    };
    //Swap (copy and swap idiom)
    void swap(SampleEvent& sample_event){
        std::swap(t, sample_event.t);
    }
    //Unify copy assignment
    SampleEvent& operator=(SampleEvent sample_event){
        sample_event.swap(*this);
        return *this;
    }
    //Clear
    virtual void clear() override{
        SampleEvent().swap(*this);
    };

    //Print sampled physical properties
    //Print kinetic temperature
    // void PrintKineticTemperature(){
    //     kinetic_termperature_output << std::setprecision(7) << kinetic_temperature << endl;
    // };

    protected:

    // double kinetic_temperature;

};

class EDMD::AndersenThermostatEvent : public Event
{
    public:

    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "andersen";
    }

    //Default constructor
    AndersenThermostatEvent() : 
        Event(),
        T_bath(0),
        lambda(0){
    };
    //Constructor
    AndersenThermostatEvent(const double& time, const double& T_bath, const double& lambda) : 
        Event(time),
        T_bath(T_bath),
        lambda(lambda){  
    };
    //Copy construcotr
    AndersenThermostatEvent(const AndersenThermostatEvent&) = default;
    //Move constructor
    AndersenThermostatEvent(AndersenThermostatEvent&&) = default;
    //Destructor
    virtual ~AndersenThermostatEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new AndersenThermostatEvent(*this));
    };
    //Swap (copy and swap idiom)
    void swap(AndersenThermostatEvent& andersen_event){
        std::swap(t, andersen_event.t);
        std::swap(T_bath, andersen_event.T_bath);
        std::swap(lambda, andersen_event.lambda);
    };
    //Unify copy assignment
    AndersenThermostatEvent& operator= (AndersenThermostatEvent andersen_event){
        andersen_event.swap(*this);
        return *this;
    };
    //Clear
    virtual void clear() override{
        AndersenThermostatEvent().swap(*this);
    };

    protected:
    double T_bath;                                                //Thermostat temperature;
    double lambda;                                                //Poisson parameter, rate of thermostat > 0
};

class EDMD::SquareWellEvent : public Event
{
    public:

    //Friend class
    friend class EDMD;
    friend class ParticleEDMD;

    //Type
    virtual string Type() const override{
        return "squarewell";
    }

    //Default constructor
    SquareWellEvent() :
        Event(),
        ptr_particle_squarewell{nullptr, nullptr},
        ptr_node_squarewell{nullptr, nullptr}{
    };
    //Constructor
    SquareWellEvent(const double& time, PtrParticleEDMD& particle_squarewell1, PtrParticleEDMD& particle_squarewell2) :
        Event(time),
        ptr_particle_squarewell{&particle_squarewell1, &particle_squarewell2},
        ptr_node_squarewell{}{ 
    };
    //Copy constructor
    SquareWellEvent(const SquareWellEvent&) = default;
    //Move constructor
    SquareWellEvent(SquareWellEvent&&) = default;
    //Destructor
    virtual ~SquareWellEvent() = default;
    //Virtual copy constructor
    virtual PtrEvent clone() const override {
        return PtrEvent(new SquareWellEvent(*this));
    };
    //Swap (copy and swap idiom) 
    void swap(SquareWellEvent& squarewell_event){
        std::swap(t, squarewell_event.t);
        std::swap(ptr_particle_squarewell, squarewell_event.ptr_particle_squarewell);
        std::swap(ptr_node_squarewell, squarewell_event.ptr_node_squarewell);
    };
    //Unify copy assignment
    SquareWellEvent& operator= (SquareWellEvent squarewell_event){
        squarewell_event.swap(*this);
        return *this;
    };
    //Clear
    virtual void clear() override{
        SquareWellEvent().swap(*this);
    };

    protected:
    array<PtrParticleEDMD*, 2> ptr_particle_squarewell;
    array<Node*, 2> ptr_node_squarewell;
};



#endif