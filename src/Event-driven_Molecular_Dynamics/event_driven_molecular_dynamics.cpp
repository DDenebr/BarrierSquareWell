#include <sstream>

#include "event_driven_molecular_dynamics.h"
#include "event.h"

EDMD::Tree EDMD::ParticleEDMD::NULL_TREE;

EDMD::EDMD(const System& system) : 
    System(system),
    event_tree(),
    virial(){
    particleEDMD_pool.reserve(MAX_SYSTEM_PARTICLE_POOL_SIZE);
    for (auto& i : particle_pool){
        particleEDMD_pool.emplace_back(make_shared<ParticleEDMD>(ParticleEDMD(*i)));
        i = particleEDMD_pool.back();
    }
    diameter_max = dmax();
    collision_search_range = static_cast<int>(std::ceil(diameter_max));
    squarewell_search_range = static_cast<int>(std::ceil(squarewell_width));
};

void EDMD::SetSamplingParameters(const string& sample_mode, const double& sample_interval, const vector<string>& sample_properties, const path& dump_directory){
    this->sampling_time = 0;
    this->sampling_mode = sample_mode;
    this->sampling_interval = sample_interval;
    this->sampling_property_list = sample_properties;
    this->dump_directory = dump_directory;
    // SampleEvent::kinetic_termperature_output.open(filename_kinetic_temperature, std::ios::trunc);
    // assert(SampleEvent::kinetic_termperature_output.is_open());
};

void EDMD::SetSamplingOperation(function<void(const double&, EDMD&)> sample_operation){
    sampling_operation = sample_operation;
}

void EDMD::SetTerminateOperation(function<void(const double&, EDMD&)> terminate_operation){
    TerminateEvent::terminate_operation = terminate_operation;
}

void EDMD::SetResetParameters(const double& time_per_epoch){
    TIME_PER_EPOCH_ = time_per_epoch;
}

void EDMD::SetResetEpoch(const double& total_time){
    epoch = 0;
    time_remainder = std::remquo(total_time, TIME_PER_EPOCH_, &max_epoch);
};

int EDMD::GetEpoch() {
    return epoch; 
};

double EDMD::GetTimePerEpoch(){
    return TIME_PER_EPOCH_;
}

double EDMD::GetSamplingTime(){
    return sampling_time;
}

void EDMD::AddCrossInParticle(const Node& node_cross){
    const shared_ptr<CrossEvent>& cross_event = std::dynamic_pointer_cast<CrossEvent>(node_cross->second);
    assert(cross_event != shared_ptr<CrossEvent>());
    //Update pointers to event node of cross event
    const PtrParticleEDMD& particle_cross = *cross_event->ptr_particle_cross;
    particle_cross->node_cross = node_cross;
};

void EDMD::AddCollideInParticle(const Node& node_collide){
    const shared_ptr<CollideEvent>& collide_event = std::dynamic_pointer_cast<CollideEvent>(node_collide->second);
    assert(collide_event != shared_ptr<CollideEvent>());
    //Update pointers to event node list in particles
    collide_event->ptr_node_collide[0] = (*collide_event->ptr_particle_collide[0])->AppendCollideNode(node_collide);
    collide_event->ptr_node_collide[1] = (*collide_event->ptr_particle_collide[1])->AppendCollideNode(node_collide);
};

void EDMD::AddSquareWellInParticle(const Node& node_squarewell){
    const shared_ptr<SquareWellEvent>& squarewell_event = std::dynamic_pointer_cast<SquareWellEvent>(node_squarewell->second);
    assert(squarewell_event != shared_ptr<SquareWellEvent>());
    //Update pointers to event node list in particles
    squarewell_event->ptr_node_squarewell[0] = (*squarewell_event->ptr_particle_squarewell[0])->AppendSquareWellNode(node_squarewell);
    squarewell_event->ptr_node_squarewell[1] = (*squarewell_event->ptr_particle_squarewell[1])->AppendSquareWellNode(node_squarewell);
};

void EDMD::RemoveCrossFromParticle(const Node& node_cross){
    const shared_ptr<CrossEvent>& cross_event = std::dynamic_pointer_cast<CrossEvent>(node_cross->second);
    assert(cross_event != shared_ptr<CrossEvent>());
    const int& pi = index(cross_event->ptr_particle_cross, particleEDMD_pool);
    //Set the cross node stored in the cross particle to be null
    particleEDMD_pool.at(pi)->node_cross = ParticleEDMD::NULL_TREE.end();
}

void EDMD::RemoveCollideFromParticle(const Node& node_collide){
    const shared_ptr<CollideEvent> collide_event = std::dynamic_pointer_cast<CollideEvent>(node_collide->second);
    assert(collide_event != shared_ptr<CollideEvent>());

    for(int i = 0; i < 2; ++ i){
        const PtrParticleEDMD& particle_collide = *collide_event->ptr_particle_collide[i];
        Node *const & ptr_node_collide = collide_event->ptr_node_collide[i];

        //Adjust pointer of particle in collide event of erased position
        Node* ptr_erased = particle_collide->EraseCollideNode(ptr_node_collide);
        if (ptr_erased != &*particle_collide->node_collide_list.end()){
            //Collide event originally at the back of the collide event list
            const shared_ptr<CollideEvent> collide_event_at_list_back = std::dynamic_pointer_cast<CollideEvent>((*ptr_erased)->second);
            //Adjust pointer from the collide event list end to the new position
            collide_event_at_list_back->ptr_node_collide[0] == &*particle_collide->node_collide_list.end() ? 
                collide_event_at_list_back->ptr_node_collide[0] = ptr_erased :
                collide_event_at_list_back->ptr_node_collide[1] = ptr_erased;
        }
    }
};

void EDMD::RemoveSquareWellFromParticle(const Node& node_squarewell){
    const shared_ptr<SquareWellEvent> squarewell_event = std::dynamic_pointer_cast<SquareWellEvent>(node_squarewell->second);
    assert(squarewell_event != shared_ptr<SquareWellEvent>());

    for(int i = 0; i < 2; ++ i){
        const PtrParticleEDMD& particle_squarewell = *squarewell_event->ptr_particle_squarewell[i];
        Node *const & ptr_node_squarewell = squarewell_event->ptr_node_squarewell[i];

        //Adjust pointer of particle in collide event of erased position
        Node* ptr_erased = particle_squarewell->EraseSquareWellNode(ptr_node_squarewell);
        if (ptr_erased != &*particle_squarewell->node_squarewell_list.end()){
            //Collide event originally at the back of the collide event list
            const shared_ptr<SquareWellEvent> squarewell_event_at_list_back = std::dynamic_pointer_cast<SquareWellEvent>((*ptr_erased)->second);
            //Adjust pointer from the collide event list end to the new position
            squarewell_event_at_list_back->ptr_node_squarewell[0] == &*particle_squarewell->node_squarewell_list.end() ? 
                squarewell_event_at_list_back->ptr_node_squarewell[0] = ptr_erased :
                squarewell_event_at_list_back->ptr_node_squarewell[1] = ptr_erased;
        }
    }
}

//Initialize cross event node
void EDMD::InitializeCrossNode(PtrParticleEDMD& particle_cross){
    const int& pi = index(&particle_cross, particleEDMD_pool);

    double t = 0;
    array<int, 2> cross_direction{};
    double Delta_t_prime = std::numeric_limits<double>::max();

    //Particle at rest
    if (particle_cross->v == array<double, 3>{})
        return;

    //Calculate the event time (crs[].t) and cross direction (crs[].crs_dir[])
    for (int i = 0; i < 2; ++i)                                                         //i=0, cross the bottom; i=1, cross the top
    for (int j = 0; j < 3; ++j)                                                         //x,y,z direction
    {
        const double& t = (i * l - particle_cross->r[j]) / (particle_cross->v[j]);      //Calculate crossing time (t) in 6 directions (+-)(xyz)
        if (t > 0 && t < Delta_t_prime){                                                //Choose the smallest positive t{
            Delta_t_prime = t;
            cross_direction[1] = -1 + 2*i; 
            cross_direction[0] = j;
        }
    }

    //Rest
    t = Delta_t_prime + particle_cross->t;
    
    // Emplace node on event tree
    PtrEvent cross_event = std::dynamic_pointer_cast<Event>(make_shared<CrossEvent>(CrossEvent(t, particle_cross, cross_direction)));
    auto itr_bool = event_tree.emplace(t, cross_event);
    while(!itr_bool.second){
        t = std::nextafter(t, std::numeric_limits<double>::max());
        cross_event = std::dynamic_pointer_cast<Event>(make_shared<CrossEvent>(CrossEvent(t, particle_cross, cross_direction)));
        itr_bool = event_tree.emplace(t, cross_event);
    }
    assert(itr_bool.second);
    AddCrossInParticle(itr_bool.first);
};

// Initialize collide event node
void EDMD::InitializeCollideNode(PtrParticleEDMD& particle_collide1, PtrParticleEDMD& particle_collide2){
    const int& pi1 = index(&particle_collide1, particleEDMD_pool);
    const int& pi2 = index(&particle_collide2, particleEDMD_pool);

    const PtrParticle& particle_collide_base1 = std::dynamic_pointer_cast<Particle>(particle_collide1);
    const PtrParticle& particle_collide_base2 = std::dynamic_pointer_cast<Particle>(particle_collide2);

    if (pi1 == pi2)
        return;
    
    const double& contact_length = 0.5 * (particle_collide1->d + particle_collide2->d);

    //Collision time, no collision when {-1, -1}
    array<double, 2> t{};
    
    t = contacttime(*particle_collide_base1, *particle_collide_base2, contact_length);

    //No collision
    if (t == array<double, 2>{-1, -1})
        return;
    
    const double& t1 = particle_collide1->t;
    const double& t2 = particle_collide2->t;

    //Overlap check
    if ((t[0] - t1 > 0 && t[1] - t1 < 0) || (t[0] - t1 < 0 && t[1] - t1 > 0) || 
        (t[0] - t2 > 0 && t[1] - t2 < 0) || (t[0] - t2 < 0 && t[1] - t2 > 0))
    if (!((t[0] - t1 < 0 && t[1] - t1 < 0) || (t[0] - t2 < 0 && t[1] - t2 < 0))){
        double t = (t1 > t2) ? t1 : t2;
        double l = distance(*particle_collide_base1, *particle_collide_base2, t);
        cerr << "Collide_Event::Initialize: Particle overlapping"<<std::endl;
        cerr << "At " << t << ", contact length = " << contact_length << ", l= "<< l <<", offset="<< contact_length - l;
        assert(false);
    }

    //No collision
    // if (t[1] == std::numeric_limits<double>::infinity())
    //     return;

    //Collision before current time point
    if ((t[1] - t1) < 0 || (t[1] - t2) < 0)
        return;

    // Emplace node on event tree
    PtrEvent collide_event = std::dynamic_pointer_cast<Event>(make_shared<CollideEvent>(CollideEvent(t[1], particle_collide1, particle_collide2)));
    auto itr_bool = event_tree.emplace(t[1], collide_event);
    while(!itr_bool.second){
        t[1] = std::nextafter(t[1], std::numeric_limits<double>::max());
        collide_event = std::dynamic_pointer_cast<Event>(make_shared<CollideEvent>(CollideEvent(t[1], particle_collide1, particle_collide2)));
        itr_bool = event_tree.emplace(t[1], collide_event);
    }
    assert(itr_bool.second);
    AddCollideInParticle(itr_bool.first);
}

//Initialize square well event node
void EDMD::InitializeSquareWellNode(PtrParticleEDMD& particle_squarewell1, PtrParticleEDMD& particle_squarewell2){
    const int& pi1 = index(&particle_squarewell1, particleEDMD_pool);
    const int& pi2 = index(&particle_squarewell2, particleEDMD_pool);

    const PtrParticle& particle_squarewell_base1 = std::dynamic_pointer_cast<Particle>(particle_squarewell1);
    const PtrParticle& particle_squarewell_base2 = std::dynamic_pointer_cast<Particle>(particle_squarewell2);

    if (pi1 == pi2)
        return;
    
    const double& contact_length = squarewell_width;

    //Square-well contact time, no contact when {-1, -1}
    array<double, 2> t{};
    
    t = contacttime(*particle_squarewell_base1, *particle_squarewell_base2, contact_length);

    //No contact
    if (t == array<double, 2>{-1, -1})
        return;
    
    const double& t1 = particle_squarewell1->t;
    const double& t2 = particle_squarewell2->t;

    const auto& sgn1 = sgn_array(t - t1);
    const auto& sgn2 = sgn_array(t - t2);

    //Contact before current time point
    if (sgn1 == array<int, 2>{-1, -1} || sgn2 == array<int, 2>{-1, -1})
        return;

    //Identify bond connection status
    double event_t;
    if (sgn1 == array<int, 2>{1, -1} || sgn2 == array<int, 2>{1, -1})
        event_t = t[0];
    else if (sgn1 == array<int, 2>{1, 1} && sgn2 == array<int, 2>{1, 1})
        event_t = t[1];
    else
        assert(("Unexpected square well time calculation", false));
    
    //Emplace node on event tree
    PtrEvent squarewell_event = std::dynamic_pointer_cast<Event>(make_shared<SquareWellEvent>(SquareWellEvent(event_t, particle_squarewell1, particle_squarewell2)));
    auto itr_bool = event_tree.emplace(event_t, squarewell_event);
    while(!itr_bool.second){
        event_t = std::nextafter(event_t, std::numeric_limits<double>::max());
        squarewell_event = std::dynamic_pointer_cast<Event>(make_shared<SquareWellEvent>(SquareWellEvent(event_t, particle_squarewell1, particle_squarewell2)));
        itr_bool = event_tree.emplace(event_t, squarewell_event);
    }
    assert(itr_bool.second);
    AddSquareWellInParticle(itr_bool.first);
}

//Initialize andersen event node
void EDMD::InitializeAndersenNode(const double& last_thermostat_time, const double& T_bath, const double& lambda){
    
    std::uniform_real_distribution<double> Ur_dice(0, 1);
    double t = last_thermostat_time - std::log(Ur_dice(gen)) / lambda;

    // Emplace node on event tree
    PtrEvent andersen_event = std::dynamic_pointer_cast<Event>(make_shared<AndersenThermostatEvent>(AndersenThermostatEvent(t, T_bath, lambda)));
    auto itr_bool = event_tree.emplace(t, andersen_event);
    while(!itr_bool.second){
        t = std::nextafter(t, std::numeric_limits<double>::max());
        andersen_event = std::dynamic_pointer_cast<Event>(make_shared<AndersenThermostatEvent>(AndersenThermostatEvent(t, T_bath, lambda)));
        itr_bool = event_tree.emplace(t, andersen_event);
    }
    assert(itr_bool.second);
}

//Initialize reset event node
void EDMD::InitializeResetNode(){
    //If reach the max epoch, initialize terminate event
    if (epoch >= max_epoch)
        InitializeTerminateNode(time_remainder);
    else{
        // Emplace node on event tree
        PtrEvent reset_event = std::dynamic_pointer_cast<Event>(make_shared<ResetEvent>(ResetEvent(TIME_PER_EPOCH_)));
        auto itr_bool = event_tree.emplace(TIME_PER_EPOCH_, reset_event);
        assert(itr_bool.second);
    }
}

//Initialize terminate event node
void EDMD::InitializeTerminateNode(const double& terminate_time){
    // Emplace node on event tree
    PtrEvent terminate_event = std::dynamic_pointer_cast<Event>(make_shared<TerminateEvent>(TerminateEvent(terminate_time)));
    auto itr_bool = event_tree.emplace(terminate_time, terminate_event);
    double t = terminate_time;
    while(!itr_bool.second){
        // t += std::numeric_limits<double>::epsilon() * std::exp2(std::floor(std::log2(t)) - 1.0) * (1.0 + std::numeric_limits<double>::epsilon());
        t = std::nextafter(t, std::numeric_limits<double>::max());
        terminate_event = std::dynamic_pointer_cast<Event>(make_shared<TerminateEvent>(TerminateEvent(t)));
        itr_bool = event_tree.emplace(t, terminate_event);
    }
    assert(itr_bool.second);
}

//Initialize sample event
void EDMD::InitializeSampleNode(const double& sampling_time){
    // Emplace node on event tree
    PtrEvent sample_event = std::dynamic_pointer_cast<Event>(make_shared<SampleEvent>(SampleEvent(sampling_time)));
    auto itr_bool = event_tree.emplace(sampling_time, sample_event);
    double t = sampling_time;
    while(!itr_bool.second){
        // t += std::numeric_limits<double>::epsilon() * std::exp2(std::floor(std::log2(t)) - 1.0) * (1.0 + std::numeric_limits<double>::epsilon());
        t = std::nextafter(t, std::numeric_limits<double>::max());
        sample_event = std::dynamic_pointer_cast<Event>(make_shared<SampleEvent>(SampleEvent(t)));
        itr_bool = event_tree.emplace(t, sample_event);
    }
    assert(itr_bool.second);
}

//Erase cross event node
void EDMD::EraseCrossNode(const PtrParticleEDMD& particle_cross){
    if (particle_cross->node_cross != ParticleEDMD::NULL_TREE.end()){
        Node node_cross(particle_cross->node_cross);
        RemoveCrossFromParticle(node_cross);
        event_tree.erase(node_cross);
    }
};

//Erase collide event node
void EDMD::EraseCollideNode(const PtrParticleEDMD& particle_collide){
    for (auto ritr = particle_collide->node_collide_list.rbegin(); ritr < particle_collide->node_collide_list.rend(); ++ ritr){
        Node node_collide(*ritr);
        RemoveCollideFromParticle(node_collide);
        event_tree.erase(node_collide);
    }
}

//Erase squarewell event node
void EDMD::EraseSquareWellNode(const PtrParticleEDMD& particle_squarewell){
    for (auto ritr = particle_squarewell->node_squarewell_list.rbegin(); ritr < particle_squarewell->node_squarewell_list.rend(); ++ ritr){
        Node node_squarewell(*ritr);
        RemoveSquareWellFromParticle(node_squarewell);
        event_tree.erase(node_squarewell);
    }
}

//Erase andersen thermostat event node
void EDMD::EraseAndersenNode(const Node& node_andersen){
    event_tree.erase(node_andersen);
}

//Erase reset event node
void EDMD::EraseResetNode(const Node& node_reset){
    event_tree.erase(node_reset);
}

//Erase terminate node
void EDMD::EraseTerminateNode(const Node& node_terminate){
    event_tree.erase(node_terminate);
}

//Erase sample event node
void EDMD::EraseSampleNode(const Node& node_sample){
    event_tree.erase(node_sample);
}

//Execute cross event node
void EDMD::ExecuteCrossNode(const Node& node){
    const shared_ptr<CrossEvent>& cross_event = std::dynamic_pointer_cast<CrossEvent>(node->second);
    assert(cross_event != shared_ptr<CrossEvent>());
    const int& pi = index(cross_event->ptr_particle_cross, particleEDMD_pool);

    PtrParticleEDMD& particle_cross = particleEDMD_pool.at(pi);
    PtrParticle& particle_cross_base = particle_pool.at(pi);

    //Remove particle from its original cell
    RemoveParticleFromCell(particle_cross_base);

    const double& Delta_t = cross_event->t - particle_cross->t;
    // const double& velocity_decay = std::exp(- gamma / m * Delta_t);

    // Update accumulated time stored in particle cross
    particle_cross->t = cross_event->t; 

    // Calculate the new particle coordinate
    for (int i = 0; i < 3; ++ i)
    if (cross_event->cross_direction[0] == i)
        particle_cross->r[i] = static_cast<double>((1 - cross_event->cross_direction[1]) / 2) * l;
    else 
        particle_cross->r[i] += particle_cross->v[i] * Delta_t;

    // Update cell information stored in particle cross
    for (int i = 0; i < 3; ++ i)
    if (cross_event->cross_direction[0] == i){
        // Update cross_record
        if (particle_cross->c[i] + cross_event->cross_direction[1] > n - 1)
            // ++ particle_cross->cross_record[i];
            ++ particle_cross->b[i];
        if (particle_cross->c[i] + cross_event->cross_direction[1] < 0)
            // -- particle_cross->cross_record[i];
            -- particle_cross->b[i];

        // Update new cell
        particle_cross->c[i] = (particle_cross->c[i] + cross_event->cross_direction[1] + n) % n;
    }

    //Add particle to the new cell
    AddParticleInCell(particle_cross_base);

    //Update squarewell event
    array<array<int, 3>, 2> search_range;
    search_range[0].fill(-squarewell_search_range);
    search_range[1].fill(squarewell_search_range);

    //Narrow the search range to the x9 cells in moving direction
    for (auto& i : search_range)
        i.at(cross_event->cross_direction[0]) = cross_event->cross_direction[1] * squarewell_search_range;
    
    //Get the current squarewell cell list of cross particle
    vector<PtrParticleEDMD> squarewell_cell_list;
    squarewell_cell_list.reserve(particle_cross->node_squarewell_list.size());
    for (const auto& i : particle_cross->node_squarewell_list){
        const shared_ptr<SquareWellEvent>& squarewell_event = std::dynamic_pointer_cast<SquareWellEvent>(i->second);
        assert(squarewell_event != shared_ptr<SquareWellEvent>());
        squarewell_cell_list.emplace_back(*squarewell_event->ptr_particle_squarewell[0] == particle_cross ? 
            *squarewell_event->ptr_particle_squarewell[1] : *squarewell_event->ptr_particle_squarewell[0]);
    }

    //Add possible new squarewell events to event tree
    for (int i = search_range[0][0]; i <= search_range[1][0]; ++ i)
    for (int j = search_range[0][1]; j <= search_range[1][1]; ++ j)
    for (int k = search_range[0][2]; k <= search_range[1][2]; ++ k)
    {
        const Cell& c = cell_grid[(particle_cross->c[0] + i + n) % n]
                                 [(particle_cross->c[1] + j + n) % n]
                                 [(particle_cross->c[2] + k + n) % n];

        for (const int& l : c.particle_list){
            PtrParticleEDMD& p1 = particle_cross;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (std::find(squarewell_cell_list.begin(), squarewell_cell_list.end(), p2) == squarewell_cell_list.end())
                InitializeSquareWellNode(p1, p2);
        }
    }

    // Reinitialize cross event
    EraseCrossNode(particle_cross);
    InitializeCrossNode(particle_cross);
};

//Execute collide node
void EDMD::ExecuteCollideNode(const Node& node){
    const shared_ptr<CollideEvent>& collide_event = std::dynamic_pointer_cast<CollideEvent>(node->second);
    assert(collide_event != shared_ptr<CollideEvent>());
    const int& pi1 = index(collide_event->ptr_particle_collide[0], particleEDMD_pool);
    const int& pi2 = index(collide_event->ptr_particle_collide[1], particleEDMD_pool);

    PtrParticleEDMD& particle_collide1 = particleEDMD_pool.at(pi1);
    PtrParticleEDMD& particle_collide2 = particleEDMD_pool.at(pi2);

    array<double, 3> dr{};
    array<double, 3> dv{};
    double drdv = 0;
    const double& contact_length = 0.5 * (particle_collide1->d + particle_collide2->d);

    for (int i = 0; i < 3; ++ i){
        // Update particles' new coordinate
        particle_collide1->r[i] += particle_collide1->v[i] * (collide_event->t - particle_collide1->t);
        particle_collide2->r[i] += particle_collide2->v[i] * (collide_event->t - particle_collide2->t);

        // Calculate Delta r and Delta v
        dr[i] = static_cast<double>(particle_collide1->c[i] - particle_collide2->c[i]) * l + 
                particle_collide1->r[i] - particle_collide2->r[i];
        if (2 * dr[i] > L)
            dr[i] = std::fmod(dr[i], L) - L;
        if (2 * dr[i] < -L)
            dr[i] = std::fmod(dr[i], L) + L;
        dv[i] = particle_collide1->v[i] - particle_collide2->v[i];

        // Calculate Delta r · Delta v
        drdv += dr[i] * dv[i];
    }

    assert(drdv < 0);

    //Normal collision with additional kinetic energy 
    assert(particle_collide1->m == particle_collide2->m);
    const double& scalar_coefficient = - drdv / contact_length / contact_length; 

    // Calculate particles' new velocities
    for (int i = 0; i < 3; ++ i){
        particle_collide1->v[i] += scalar_coefficient * dr[i];
        particle_collide2->v[i] -= scalar_coefficient * dr[i];
    }

    //Update particles' time
    particle_collide1->t = collide_event->t;
    particle_collide2->t = collide_event->t;

    //Update virial
    if (virial.size() >= MAX_VIRIAL_SIZE)
        virial.pop_front();
    virial.emplace_back(std::make_pair(collide_event->t, scalar_coefficient * dot(dr, dr)));
        
    //Update all Cross/Collide/Rest events for particles
    //Update cross event
    EraseCrossNode(particle_collide1);
    EraseCrossNode(particle_collide2);

    InitializeCrossNode(particle_collide1);
    InitializeCrossNode(particle_collide2);

    // Update collide event
    EraseCollideNode(particle_collide1);
    EraseCollideNode(particle_collide2);

    for (const auto& neighbor_ptr : particle_collide1->neighbor_list){
        PtrParticleEDMD& p1 = particle_collide1;
        PtrParticleEDMD& p2 = *neighbor_ptr;
        if (p2 != particle_collide2)
            InitializeCollideNode(p1, p2);  
    }

    for (const auto& neighbor_ptr : particle_collide2->neighbor_list){
        PtrParticleEDMD& p1 = particle_collide2;
        PtrParticleEDMD& p2 = *neighbor_ptr;
        if (p2 != particle_collide1)
            InitializeCollideNode(p1, p2);  
    }

    // Update squarewell event
    EraseSquareWellNode(particle_collide1);
    EraseSquareWellNode(particle_collide2);

    InitializeSquareWellNode(particle_collide1, particle_collide2);
    for (int i = -squarewell_search_range; i <= squarewell_search_range; ++ i)
    for (int j = -squarewell_search_range; j <= squarewell_search_range; ++ j)
    for (int k = -squarewell_search_range; k <= squarewell_search_range; ++ k){
        const Cell& c1 = cell_grid[(particle_collide1->c[0] + i + n) % n]
                                [(particle_collide1->c[1] + j + n) % n]
                                [(particle_collide1->c[2] + k + n) % n];

        //New possible collision of collide particle 1 
        for (const int& l : c1.particle_list){
            PtrParticleEDMD& p1 = particle_collide1;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (p2 != particle_collide1 && p2 != particle_collide2)
                InitializeSquareWellNode(p1, p2);    
        }
        
        const Cell& c2 = cell_grid[(particle_collide2->c[0] + i + n) % n]
                                [(particle_collide2->c[1] + j + n) % n]
                                [(particle_collide2->c[2] + k + n) % n];
        
        //New possible collision of collide particle 2 
        for (const int& l : c2.particle_list){
            PtrParticleEDMD& p1 = particle_collide2;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (p2 != particle_collide1 && p2 != particle_collide2)
                InitializeSquareWellNode(p1, p2);
        }    
    }
}

void EDMD::ExecuteSquareWellNode(const Node& node){
    const shared_ptr<SquareWellEvent>& squarewell_event = std::dynamic_pointer_cast<SquareWellEvent>(node->second);
    assert(squarewell_event != shared_ptr<SquareWellEvent>());
    const int& pi1 = index(squarewell_event->ptr_particle_squarewell[0], particleEDMD_pool);
    const int& pi2 = index(squarewell_event->ptr_particle_squarewell[1], particleEDMD_pool);

    PtrParticleEDMD& particle_squarewell1 = particleEDMD_pool.at(pi1);
    PtrParticleEDMD& particle_squarewell2 = particleEDMD_pool.at(pi2);

    array<double, 3> dr{};
    array<double, 3> dv{};
    double drdv = 0;
    const double& contact_length = squarewell_width;

    for (int i = 0; i < 3; ++ i){
        // Update particles' new coordinate
        particle_squarewell1->r[i] += particle_squarewell1->v[i] * (squarewell_event->t - particle_squarewell1->t);
        particle_squarewell2->r[i] += particle_squarewell2->v[i] * (squarewell_event->t - particle_squarewell2->t);

        // Calculate Delta r and Delta v
        dr[i] = static_cast<double>(particle_squarewell1->c[i] - particle_squarewell2->c[i]) * l + 
                particle_squarewell1->r[i] - particle_squarewell2->r[i];
        if (2 * dr[i] > L)
            dr[i] = std::fmod(dr[i], L) - L;
        if (2 * dr[i] < -L)
            dr[i] = std::fmod(dr[i], L) + L;
        dv[i] = particle_squarewell1->v[i] - particle_squarewell2->v[i];

        // Calculate Delta r · Delta v
        drdv += dr[i] * dv[i];
    }

    //Discrimiant for barrier cross
    assert(particle_squarewell1->m == particle_squarewell2->m);
    const double& m = particle_squarewell1->m;
    const int& sign1 = sgn(drdv);                               //cross in -,cross out +

    const double& Delta_barrier = drdv * drdv - 4 * contact_length * contact_length * (squarewell_barrier - (sign1 > 0) * squarewell_depth) / m;
    const int& sign2 = sgn(Delta_barrier);                      //not cross -, cross +

    const double& Delta = drdv * drdv + 4 * sign1 * (sign2 > 0) * contact_length * contact_length * squarewell_depth / m;

    const double& scalar_coefficient = (- drdv + sign1 * sign2 * std::sqrt(Delta)) / (2.0 * contact_length * contact_length);

    
    // in (bond == -1, barrier == 1)
    // (- drdv - std::sqrt(Delta)) / 2.0 / contact_length / contact_length; -
    // out (bond == 1, barrier == 1)
    // (- drdv + std::sqrt(Delta)) / 2.0 / contact_length / contact_length; +

    // repel (bond == -1, barrier == 0)
    // - drdv / contact_length / contact_length; +
    // rebound (bond == 1, barrier == 0)
    // - drdv / contact_length / contact_length; -

    // Calculate particles' new velocities
    for (int i = 0; i < 3; ++ i){
        particle_squarewell1->v[i] += scalar_coefficient * dr[i];
        particle_squarewell2->v[i] -= scalar_coefficient * dr[i];
    }
    
    //Update particles' time
    particle_squarewell1->t = squarewell_event->t;
    particle_squarewell2->t = squarewell_event->t;

    //Update virial
    if (virial.size() >= MAX_VIRIAL_SIZE)
        virial.pop_front();
    virial.emplace_back(std::make_pair(squarewell_event->t, scalar_coefficient * dot(dr, dr)));

    //Update neighbor list
    if (sign2 == 1)
    //Add neighbors
    if (sign1 == -1){
        vector<PtrParticleEDMD*>& n1 = particle_squarewell1->neighbor_list;
        vector<PtrParticleEDMD*>& n2 = particle_squarewell2->neighbor_list;
        assert(std::find(n1.begin(), n1.end(), &particle_squarewell2) == n1.end());
        assert(std::find(n2.begin(), n2.end(), &particle_squarewell1) == n2.end());

        n1.emplace_back(&particle_squarewell2);
        n2.emplace_back(&particle_squarewell1);
    }
    //Erase neighbor
    else{
        vector<PtrParticleEDMD*>& n1 = particle_squarewell1->neighbor_list;
        vector<PtrParticleEDMD*>& n2 = particle_squarewell2->neighbor_list;
        assert(std::find(n1.begin(), n1.end(), &particle_squarewell2) != n1.end());
        assert(std::find(n2.begin(), n2.end(), &particle_squarewell1) != n2.end());

        n1.erase(std::find(n1.begin(), n1.end(), &particle_squarewell2));
        n2.erase(std::find(n2.begin(), n2.end(), &particle_squarewell1));
    }
        
    //Update all Cross/Collide/Rest events for particles
    //Update cross event
    EraseCrossNode(particle_squarewell1);
    EraseCrossNode(particle_squarewell2);

    InitializeCrossNode(particle_squarewell1);
    InitializeCrossNode(particle_squarewell2);

    //Update collide event
    EraseCollideNode(particle_squarewell1);
    EraseCollideNode(particle_squarewell2);

    InitializeCollideNode(particle_squarewell1, particle_squarewell2);
    for (const auto& neighbor_ptr : particle_squarewell1->neighbor_list){
        PtrParticleEDMD& p1 = particle_squarewell1;
        PtrParticleEDMD& p2 = *neighbor_ptr;
        if (p2 != particle_squarewell2)
        InitializeCollideNode(p1, p2);  
    }
    for (const auto& neighbor_ptr : particle_squarewell2->neighbor_list){
        PtrParticleEDMD& p1 = particle_squarewell2;
        PtrParticleEDMD& p2 = *neighbor_ptr;
        if (p2 != particle_squarewell1)
        InitializeCollideNode(p1, p2);  
    }

    //Update squarewell event
    EraseSquareWellNode(particle_squarewell1);
    EraseSquareWellNode(particle_squarewell2);

    for (int i = -squarewell_search_range; i <= squarewell_search_range; ++ i)
    for (int j = -squarewell_search_range; j <= squarewell_search_range; ++ j)
    for (int k = -squarewell_search_range; k <= squarewell_search_range; ++ k){
        const Cell& c1 = cell_grid[(particle_squarewell1->c[0] + i + n) % n]
                                [(particle_squarewell1->c[1] + j + n) % n]
                                [(particle_squarewell1->c[2] + k + n) % n];

        //New possible collision of collide particle 1 
        for (const int& l : c1.particle_list){
            PtrParticleEDMD& p1 = particle_squarewell1;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (p2 != particle_squarewell1 && p2 != particle_squarewell2)
                InitializeSquareWellNode(p1, p2);    
        }
        
        const Cell& c2 = cell_grid[(particle_squarewell2->c[0] + i + n) % n]
                                [(particle_squarewell2->c[1] + j + n) % n]
                                [(particle_squarewell2->c[2] + k + n) % n];
        
        //New possible collision of collide particle 2 
        for (const int& l : c2.particle_list){
            PtrParticleEDMD& p1 = particle_squarewell2;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (p2 != particle_squarewell1 && p2 != particle_squarewell2)
                InitializeSquareWellNode(p1, p2);
        }    
    }

    //Mannually update boundary scenario
    if (sign1 * sign2 == -1){
        array<double, 2> t = contacttime(*particle_squarewell1, *particle_squarewell2, squarewell_width);
        assert(t[0] != -1);

        //Emplace node on event tree
        PtrEvent squarewell_event = std::dynamic_pointer_cast<Event>(make_shared<SquareWellEvent>(SquareWellEvent(t[0], particle_squarewell1, particle_squarewell2)));
        auto itr_bool = event_tree.emplace(t[0], squarewell_event);
        while(!itr_bool.second){
            t[0] = std::nextafter(t[0], std::numeric_limits<double>::max());
            squarewell_event = std::dynamic_pointer_cast<Event>(make_shared<SquareWellEvent>(SquareWellEvent(t[0], particle_squarewell1, particle_squarewell2)));
            itr_bool = event_tree.emplace(t[0], squarewell_event);
        }
        assert(itr_bool.second);
        AddSquareWellInParticle(itr_bool.first);
    }
}

//Execute andersen thermostat event
void EDMD::ExecuteAndersenNode(const Node& node){
    const shared_ptr<AndersenThermostatEvent>& andersen_event = std::dynamic_pointer_cast<AndersenThermostatEvent>(node->second);
    assert(andersen_event != shared_ptr<AndersenThermostatEvent>());

    std::uniform_int_distribution<int> Ui_dice(0, particleEDMD_pool.size() - 1);

    //Choose a thermostat particle
    const int& pi = Ui_dice(gen);
    PtrParticleEDMD& particle_andersen = particleEDMD_pool.at(pi);

    //Update the position of the particle
    particle_andersen->r = particle_andersen->r + particle_andersen->v * (andersen_event->t - particle_andersen->t);

    //Re-assign velocity
    const double& T_bath = andersen_event->T_bath;
    const double& lambda = andersen_event->lambda;
    const double& m = particle_andersen->m;

    std::normal_distribution<double> Nr_dice(0, std::sqrt(kB * T_bath / m));

    particle_andersen->velocity(array<double, 3>{Nr_dice(gen), Nr_dice(gen), Nr_dice(gen)});

    //Update particles' time
    particle_andersen->t = andersen_event->t;
        
    //Update all Cross/Collide/Squarewell events for particles
    //Update cross event
    EraseCrossNode(particle_andersen);
    InitializeCrossNode(particle_andersen);

    //Update collide event
    EraseCollideNode(particle_andersen);

    for (const auto& neighbor_ptr : particle_andersen->neighbor_list){
        PtrParticleEDMD& p1 = particle_andersen;
        PtrParticleEDMD& p2 = *neighbor_ptr;
        InitializeCollideNode(p1, p2);  
    }

    //Update squarewell event
    EraseSquareWellNode(particle_andersen);
    
    for (int i = -squarewell_search_range; i <= squarewell_search_range; ++ i)
    for (int j = -squarewell_search_range; j <= squarewell_search_range; ++ j)
    for (int k = -squarewell_search_range; k <= squarewell_search_range; ++ k){
        const Cell& c = cell_grid[(particle_andersen->c[0] + i + n) % n]
                                 [(particle_andersen->c[1] + j + n) % n]
                                 [(particle_andersen->c[2] + k + n) % n];

        //New possible collision of andersen particle
        for (const int& l : c.particle_list){
            PtrParticleEDMD& p1 = particle_andersen;
            PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
            if (p2 != particle_andersen)
                InitializeSquareWellNode(p1, p2);    
        }
    }

    //Re-initialize Andersen thermostat
    EraseAndersenNode(node);
    InitializeAndersenNode(andersen_event->t, T_bath, lambda);
}

//Execute reset event
void EDMD::ExecuteResetNode(const Node& node){
    const shared_ptr<ResetEvent>& reset_event = std::dynamic_pointer_cast<ResetEvent>(node->second);
    assert(reset_event != shared_ptr<ResetEvent>());

    //Modify the particle
    for (auto& pi : particleEDMD_pool){
        const double& Delta_t = reset_event->t - pi->t;

        // Update accumulated time stored in particle cross
        pi->t = 0; 

        for (int i = 0; i < 3; ++ i){
            // Calculate the new particle coordinate
            pi->r[i] += pi->v[i] * Delta_t;
        }            
    }
    
    //Reset the time to decrease numeric error
    if (epoch < max_epoch){
        //Modify the event tree
        for (auto i = event_tree.begin(); i != event_tree.end(); ++ i){
            auto node = event_tree.extract(i);
            node.key() -= TIME_PER_EPOCH_;
            node.mapped()->t -= TIME_PER_EPOCH_;
            event_tree.insert(std::move(node));
        }
        //Modify the virial list
        for (auto& [t, v] : virial)
            t -= TIME_PER_EPOCH_;

        EraseResetNode(node);
        ++ epoch;
        InitializeResetNode();
    }
    else{
        //Clear the event tree
        event_tree.clear();
        //Close sample data output file
        // SampleEvent::kinetic_termperature_output.close();
    }
}

void EDMD::ExecuteTerminateNode(const Node& node){
    const shared_ptr<TerminateEvent>& terminate_event = std::dynamic_pointer_cast<TerminateEvent>(node->second);
    assert(terminate_event != shared_ptr<TerminateEvent>());
    
    //Terminate operation
    TerminateEvent::terminate_operation(terminate_event->t, *this);

    event_tree.clear();
}

//Execute sample event
void EDMD::ExecuteSampleNode(const Node& node){
    const shared_ptr<SampleEvent>& sample_event = std::dynamic_pointer_cast<SampleEvent>(node->second);
    assert(sample_event != shared_ptr<SampleEvent>());

    /*-------------------------Sampling Desired Data-------------------------*/
    ++ sampling_time;
    sampling_operation(sample_event->t, *this);
    // const auto& properties = SampleEvent::sampling_property_list;
    // map<string, double> property_map;
    
    /*-----------------------------------------------------------------------*/
    const double& last_sampling_time = sample_event->t;
    EraseSampleNode(node);

    if (sampling_mode == "Arithmetic")
        //Arithmetic sequences sampling time
        InitializeSampleNode(last_sampling_time + sampling_interval);
    else if (sampling_mode == "Geometric"){
        //Geometric sequences sampling time
        InitializeSampleNode(last_sampling_time + sampling_interval);
        sampling_interval *= std::pow(10.0, 0.01);
    }
    else 
        assert(("unknown sampling mode", false));
        
}


//Execute event tree
void EDMD::ExecuteEventTree(){
    //Check if the tree is empty
    if (event_tree.empty())
        return;
    const Node& node_execute = event_tree.begin();
    const PtrEvent& event_execute = node_execute->second;

    // cout << event_execute->t << ' ' << event_execute->Type() << endl;

    if (event_execute->Type() == "andersen")
        ExecuteAndersenNode(node_execute);
    else if (event_execute->Type() == "cross")
        ExecuteCrossNode(node_execute);
    else if (event_execute->Type() == "collide")
        ExecuteCollideNode(node_execute);
    else if (event_execute->Type() == "reset")
        ExecuteResetNode(node_execute);
    else if (event_execute->Type() == "terminate")
        ExecuteTerminateNode(node_execute);
    else if (event_execute->Type() == "sample")
        ExecuteSampleNode(node_execute);
    else if (event_execute->Type() == "squarewell")
        ExecuteSquareWellNode(node_execute);
    else
        assert(("Unexpected event type", false));
}

void EDMD::ExecuteEDMD(const double& termination_time, const double& bath_temperature, const double& thermostat_frequency, const double& initial_sample_time){
    
    //Initialize event tree
    event_tree.clear();
    //Andersen thermostat
    InitializeAndersenNode(0, bath_temperature, thermostat_frequency);
    //Reset event
    SetResetEpoch(termination_time);
    InitializeResetNode();
    //Sample event
    InitializeSampleNode(initial_sample_time);

    for (PtrParticleEDMD& pi : particleEDMD_pool){
        //Rest events
        // InitializeRestNode(pi);
        //Cross events
        InitializeCrossNode(pi);
        //Collide events (neighbor list)
        // #pragma omp parallel for collapse(3) shared(particleEDMD_pool, cell_grid, event_tree)
        for (int i =- collision_search_range; i <= collision_search_range; ++ i)
        for (int j =- collision_search_range; j <= collision_search_range; ++ j)
        for (int k =- collision_search_range; k <= collision_search_range; ++ k){
            const System::Cell& c = cell_grid[(pi->c[0] + i + n) % n][(pi->c[1] + j + n) % n][(pi->c[2] + k + n) % n];

            for (const int& l : c.particle_list)
            if (l > index(&pi, particleEDMD_pool))
            if (distance(*pi, *particleEDMD_pool.at(l), 0) < squarewell_width){
                pi->neighbor_list.emplace_back(&particleEDMD_pool.at(l));
                particleEDMD_pool.at(l)->neighbor_list.emplace_back(&pi);

                InitializeCollideNode(pi, particleEDMD_pool.at(l));
            }
        }

        //Squarewell events
        for (int i =- squarewell_search_range; i <= squarewell_search_range; ++ i)
        for (int j =- squarewell_search_range; j <= squarewell_search_range; ++ j)
        for (int k =- squarewell_search_range; k <= squarewell_search_range; ++ k){
            const System::Cell& c = cell_grid[(pi->c[0] + i + n) % n][(pi->c[1] + j + n) % n][(pi->c[2] + k + n) % n];

            for (const int& l : c.particle_list)
            if (l > index(&pi, particleEDMD_pool)){
                InitializeSquareWellNode(pi, particleEDMD_pool.at(l));
            }
        }
    }

    // cout << event_tree.size() << endl;

    int node_execute_times = 0;

    //Execute event tree
    while(!event_tree.empty()){
        ExecuteEventTree();
        ++ node_execute_times;
    }
};

double EDMD::Pressure(const double& t){
    //Pressure
    //Kinetic temperature
    const double& kinetic_temperature = KineticTemperature(t);

    //Virial
    long double virial_part(0);
    for (auto& [t, v] : virial)
        virial_part += v;
    virial_part /= 3 * (virial.back().first - virial.front().first);

    return (kinetic_temperature * particleEDMD_pool.size() * kB + virial_part) / (L * L * L); 
}


double EDMD::MaximumBubbleVolume(const double& t){

    //resolution
    const int& nb = static_cast<int>(std::floor(L / diameter_max) * 2);
    const double& lb = L / nb;
    
    const int& vapor_neighbor_threshold = 5;
    const int& bubble_neighbor_threshold = 7;

    struct LocalInfo
    {
        //0 liquid, 1 vapor, 2 vapor after 1 neighbor check, 3 vapor afte 2 neighbor check
        int status;
        //number of vapor neighbor
        int vapor_neighbor;
        //index of bubble
        int bubble_index;

        LocalInfo(const int& s, const int& n, const int& i){
            status = s; vapor_neighbor = n; bubble_index = i;
        };
    };
    
    //Initialize LocalInfo as vapor
    vector<vector<vector<LocalInfo>>> bubble_grid(nb, vector<vector<LocalInfo>>(nb, vector<LocalInfo>(nb, LocalInfo(1, -1, 0))));

    //Preliminary calculate vapor bubble
    // #pragma omp parallel for collapse(3) num_threads(8)
    for (int bi = 0; bi < nb; ++ bi)
    for (int bj = 0; bj < nb; ++ bj)
    for (int bk = 0; bk < nb; ++ bk){

        const array<double, 3>& position{bi * lb, bj * lb, bk * lb};
        const array<int, 3>& c{static_cast<int>(bi * lb / l), static_cast<int>(bj * lb / l), static_cast<int>(bk * lb / l)};
        const array<double, 3> r{bi * lb / l - c[0], bj * lb / l - c[1], bk * lb / l - c[2]};
        Particle test_point("test_point", 0, 0, t);
        test_point.location(c, r);

        bool flag = true;
        for (int i =- 1; i <= 1 && flag; ++ i)
        for (int j =- 1; j <= 1 && flag; ++ j)
        for (int k =- 1; k <= 1 && flag; ++ k){
            const auto& nc = cell_grid[(c[0] + i + n) % n][(c[1] + j + n) % n][(c[2] + k + n) % n];

            for (const auto& l : nc.particle_list){
                const auto& np = particleEDMD_pool[l];
                const auto& d = distance(test_point, *np, t);
                //Determine if the test point is vapor or liquid
                if (d < np->d && np->neighbor_list.size() >= vapor_neighbor_threshold){
                    
                    // #pragma omp critical
                    bubble_grid[bi][bj][bk].status = 0;

                    flag = false;
                    break;
                }
            }
        }
    }

    //Calculate neighbor
    for (int bi = 0; bi < nb; ++ bi)
    for (int bj = 0; bj < nb; ++ bj)
    for (int bk = 0; bk < nb; ++ bk)
    if (bubble_grid[bi][bj][bk].status == 1)
        for (int i =- 1; i <= 1; ++ i)
        for (int j =- 1; j <= 1; ++ j)
        for (int k =- 1; k <= 1; ++ k)
            bubble_grid[bi][bj][bk].vapor_neighbor += bubble_grid[(bi + i + nb) % nb][(bj + j + nb) % nb][(bk + k + nb) % nb].status;

    //Neighbor test (>threshold)
    for (int bi = 0; bi < nb; ++ bi)
    for (int bj = 0; bj < nb; ++ bj)
    for (int bk = 0; bk < nb; ++ bk)
    if (bubble_grid[bi][bj][bk].status > 0){
        if (bubble_grid[bi][bj][bk].vapor_neighbor > bubble_neighbor_threshold)
            ++ bubble_grid[bi][bj][bk].status;
        
        int second_neighbor_counts = 0;
        for (int i =- 1; i <= 1; ++ i)
        for (int j =- 1; j <= 1; ++ j)
        for (int k =- 1; k <= 1; ++ k)
        if (bubble_grid[(bi + i + nb) % nb][(bj + j + nb) % nb][(bk + k + nb) % nb].vapor_neighbor > bubble_neighbor_threshold)
            ++ second_neighbor_counts;
        if (second_neighbor_counts > bubble_neighbor_threshold)
            ++ bubble_grid[bi][bj][bk].status;
    }

    //Calculate bubble volume
    // vector<vector<vector<int>>> bubble_list(nb, vector<vector<int>>(nb, vector<int>(nb, 0)));
    vector<int> bubble_volume_list;
    bubble_volume_list.reserve(32);

    int bubble_index = 1;
    for (int bi = 0; bi < nb; ++ bi)
    for (int bj = 0; bj < nb; ++ bj)
    for (int bk = 0; bk < nb; ++ bk)
    if (bubble_grid[bi][bj][bk].status == 3)
    if (bubble_grid[bi][bj][bk].bubble_index == 0){
        
        //Breadth first search
        list<array<int, 3>> queue;
        //Origin
        queue.emplace_back(array<int, 3>{bi, bj, bk});
        bubble_grid[bi][bj][bk].bubble_index = bubble_index;
        auto cursor = queue.begin();
        //Search
        while (cursor != queue.end()){
            const auto& current_vertice = *cursor;
            bubble_grid[current_vertice[0]][current_vertice[1]][current_vertice[2]].bubble_index = bubble_index;
            for (int i =- 1; i <= 1; ++ i)
            for (int j =- 1; j <= 1; ++ j)
            for (int k =- 1; k <= 1; ++ k){
                const int& cx = (current_vertice[0] + i + nb) % nb;
                const int& cy = (current_vertice[1] + j + nb) % nb;
                const int& cz = (current_vertice[2] + k + nb) % nb;
                if (bubble_grid[cx][cy][cz].status == 3)
                if (std::find(queue.begin(), queue.end(), array<int, 3>{cx, cy, cz}) == queue.end()){
                    queue.emplace_back(array<int, 3>{cx, cy, cz});
                }
            }
            ++ cursor;
        }
        ++ bubble_index;
        bubble_volume_list.emplace_back(queue.size());
    }

    return *std::max_element(bubble_volume_list.begin(), bubble_volume_list.end()) * lb * lb * lb;
}
