#include <sstream>

#include "event_driven_molecular_dynamics.h"
#include "event.h"

EDMD::Tree EDMD::ParticleEDMD::NULL_TREE;

EDMD::EDMD(const System& system) : 
    System(system),
    event_tree(){
    particleEDMD_pool.reserve(MAX_SYSTEM_PARTICLE_POOL_SIZE);
    for (auto& i : particle_pool){
        particleEDMD_pool.emplace_back(make_shared<ParticleEDMD>(ParticleEDMD(*i)));
        i = particleEDMD_pool.back();
    }
    diameter_max = dmax();
    collision_search_range = static_cast<int>(std::ceil(diameter_max));
    squarewell_search_range = static_cast<int>(std::ceil(squarewell_width));
};

void EDMD::SetSamplingParameters(const double& sample_interval, const path& dump_directory, const string& filename_kinetic_temperature){
    SampleEvent::sampling_time = 0;
    SampleEvent::sampling_interval = sample_interval;
    SampleEvent::dump_directory = dump_directory;
    // SampleEvent::kinetic_termperature_output.open(filename_kinetic_temperature, std::ios::trunc);
    // assert(SampleEvent::kinetic_termperature_output.is_open());
};

void EDMD::AddCrossInParticle(const Node& node_cross){
    const shared_ptr<CrossEvent>& cross_event = std::dynamic_pointer_cast<CrossEvent>(node_cross->second);
    assert(cross_event != shared_ptr<CrossEvent>());
    //Update pointers to event node of cross event
    const PtrParticleEDMD& particle_cross = *cross_event->ptr_particle_cross;
    particle_cross->node_cross = node_cross;
};

// void EDMD::AddRestInParticle(const Node& node_rest){
//     const shared_ptr<RestEvent>& rest_event = std::dynamic_pointer_cast<RestEvent>(node_rest->second);
//     assert(rest_event != shared_ptr<RestEvent>());  
//     //Update pointers to event node of cross event
//     const PtrParticleEDMD& particle_rest = *rest_event->ptr_particle_rest;
//     particle_rest->node_rest = node_rest; 
// }

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

// void EDMD::RemoveRestFromParticle(const Node& node_rest){
//     const shared_ptr<RestEvent>& rest_event = std::dynamic_pointer_cast<RestEvent>(node_rest->second);
//     assert(rest_event != shared_ptr<RestEvent>());
//     const int& pi = index(rest_event->ptr_particle_rest, particleEDMD_pool);
//     //Set the rest node stored in the rest particle to be null
//     particleEDMD_pool.at(pi)->node_rest = ParticleEDMD::NULL_TREE.end();
// }

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

//Initialize rest event node
// void EDMD::InitializeRestNode(PtrParticleEDMD& particle_rest){
//     //if no damping
//     if (gamma == 0)
//         return;

//     const int& pi = index(&particle_rest, particleEDMD_pool);

//     //Velocity to the 10^-15 of the original (double accuracy limit)
//     const double& v = std::sqrt(std::accumulate(particle_rest->v.begin(), particle_rest->v.end(), 0.0,
//         [](const double& sum_of_square, const double& i){ return sum_of_square + i * i;}));
//     //Particle already at rest
//     if (v == 0)
//         return;

//     const double& Delta_t = m / gamma * (std::log(v * m / gamma) + 15.0 * std::log(10.0));
//     double t = Delta_t + particle_rest->t;

//     // Emplace node on event tree
//     PtrEvent rest_event = std::dynamic_pointer_cast<Event>(make_shared<RestEvent>(RestEvent(t, particle_rest)));
//     auto itr_bool = event_tree.emplace(t, rest_event);
//     while(!itr_bool.second){
//         t = std::nextafter(t, std::numeric_limits<double>::max());
//         rest_event = std::dynamic_pointer_cast<Event>(make_shared<RestEvent>(RestEvent(t, particle_rest)));
//         itr_bool = event_tree.emplace(t, rest_event);
//     }
//     assert(itr_bool.second);
//     AddRestInParticle(itr_bool.first);
// };

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

    //Mistake check
    assert(("Omit particle entering the square well",
            !((sgn1 == array<int, 2>{1, 1} && sgn2 == array<int, 2>{1, -1}) || (sgn1 == array<int, 2>{1, -1} && sgn2 == array<int, 2>{1, 1}))));

    //Overlap check
    // if ((t[0] - t1 > 0 && t[1] - t1 < 0) || (t[0] - t1 < 0 && t[1] - t1 > 0) || 
    //     (t[0] - t2 > 0 && t[1] - t2 < 0) || (t[0] - t2 < 0 && t[1] - t2 > 0))
    // if (!((t[0] - t1 < 0 && t[1] - t1 < 0) || (t[0] - t2 < 0 && t[1] - t2 < 0))){
    //     double t = (t1 > t2) ? t1 : t2;
    //     double l = distance(*particle_collide_base1, *particle_collide_base2, t);
    //     cerr << "Collide_Event::Initialize: Particle overlapping"<<std::endl;
    //     cerr << "At " << t << ", contact length = " << contact_length << ", l= "<< l <<", offset="<< contact_length - l;
    //     assert(false);
    // }

    //No collision
    // if (t[1] == std::numeric_limits<double>::infinity())
    //     return;

    //Contact before current time point
    if (sgn1 == array<int, 2>{-1, -1} || sgn2 == array<int, 2>{-1, -1})
        return;

    //Identify bond connection status
    double event_t;
    if (sgn1 == array<int, 2>{1, -1} && sgn2 == array<int, 2>{1, -1})
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
void EDMD::InitializeAndersenNode(PtrParticleEDMD& particle_andersen, const double& last_thermostat_time, const double& T_bath, const double& lambda){
    
    std::uniform_real_distribution<double> Ur_dice(0, 1);
    double t = last_thermostat_time - std::log(Ur_dice(gen)) / lambda;

    // Emplace node on event tree
    PtrEvent andersen_event = std::dynamic_pointer_cast<Event>(make_shared<AndersenThermostatEvent>(AndersenThermostatEvent(t, T_bath, lambda, particle_andersen)));
    auto itr_bool = event_tree.emplace(t, andersen_event);
    while(!itr_bool.second){
        t = std::nextafter(t, std::numeric_limits<double>::max());
        andersen_event = std::dynamic_pointer_cast<Event>(make_shared<AndersenThermostatEvent>(AndersenThermostatEvent(t, T_bath, lambda, particle_andersen)));
        itr_bool = event_tree.emplace(t, andersen_event);
    }
    assert(itr_bool.second);
}

//Initialize reset event node
void EDMD::InitializeResetNode(){
    // Emplace node on event tree
    const PtrEvent& reset_event = std::dynamic_pointer_cast<Event>(make_shared<ResetEvent>(ResetEvent()));
    auto itr_bool = event_tree.emplace(ResetEvent::TIME_PER_EPOCH_, reset_event);
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

//Erase rest event node
// void EDMD::EraseRestNode(const PtrParticleEDMD& particle_rest){
//     if (particle_rest->node_rest != ParticleEDMD::NULL_TREE.end()){
//         Node node_rest(particle_rest->node_rest);
//         RemoveRestFromParticle(node_rest);
//         event_tree.erase(node_rest);
//     }
// };

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

    // Calculate the new particle velocity
    // for (int i = 0; i < 3; ++ i)
    //     particle_cross->v[i] *= velocity_decay;

    //Add particle to the new cell
    AddParticleInCell(particle_cross_base);

    // array<array<int, 3>, 2> search_range;
    // search_range[0].fill(-collision_search_range);
    // search_range[1].fill(collision_search_range);

    // //Narrow the search range to the x9 cells in moving direction
    // for (auto& i : search_range)
    //     i.at(cross_event->cross_direction[0]) = cross_event->cross_direction[1] * collision_search_range;
    
    // //Get the current collision cell list of cross particle
    // vector<PtrParticleEDMD> collision_cell_list;
    // collision_cell_list.reserve(particle_cross->node_collide_list.size());
    // for (const auto& i : particle_cross->node_collide_list){
    //     const shared_ptr<CollideEvent>& collide_event = std::dynamic_pointer_cast<CollideEvent>(i->second);
    //     collision_cell_list.emplace_back(*collide_event->ptr_particle_collide[0] == particle_cross ? 
    //         *collide_event->ptr_particle_collide[1] : *collide_event->ptr_particle_collide[0]);
    // }

    // //Add possible new collide events to event tree
    // for (int i = search_range[0][0]; i <= search_range[1][0]; ++ i)
    // for (int j = search_range[0][1]; j <= search_range[1][1]; ++ j)
    // for (int k = search_range[0][2]; k <= search_range[1][2]; ++ k)
    // {
    //     const Cell& c = cell_grid[(particle_cross->c[0] + i + n) % n]
    //                              [(particle_cross->c[1] + j + n) % n]
    //                              [(particle_cross->c[2] + k + n) % n];

    //     for (const int& l : c.particle_list){
    //         PtrParticleEDMD& p1 = particle_cross;
    //         PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
    //         if (std::find(collision_cell_list.begin(), collision_cell_list.end(), p2) == collision_cell_list.end())
    //             InitializeCollideNode(p1, p2);
    //     }
    // }

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
    for (const auto& i : particle_cross->node_collide_list){
        const shared_ptr<SquareWellEvent>& squarewell_event = std::dynamic_pointer_cast<SquareWellEvent>(i->second);
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

//Execute rest node
// void EDMD::ExecuteRestNode(const Node& node){
//     const shared_ptr<RestEvent>& rest_event = std::dynamic_pointer_cast<RestEvent>(node->second);
//     assert(rest_event != shared_ptr<RestEvent>());
//     const int& pi = index(rest_event->ptr_particle_rest, particleEDMD_pool);

//     PtrParticleEDMD& particle_rest = particleEDMD_pool.at(pi);
//     PtrParticle& particle_rest_base = particle_pool.at(pi);

//     // Update accumulated time stored in particle cross
//     particle_rest->t = rest_event->t; 

//     // Calculate the new particle coordinate
//     for (int i = 0; i < 3; ++ i)
//         particle_rest->r[i] += particle_rest->v[i] * m / gamma;

//     // Calculate the new particle velocity
//     for (int i = 0; i < 3; ++ i)
//         particle_rest->v[i] = 0;
    
//     //Erase rest node
//     EraseRestNode(particle_rest);
// }

//Execute collide node
void EDMD::ExecuteCollideNode(const Node& node){
    const shared_ptr<CollideEvent>& collide_event = std::dynamic_pointer_cast<CollideEvent>(node->second);
    assert(collide_event != shared_ptr<CollideEvent>());
    const int& pi1 = index(collide_event->ptr_particle_collide[0], particleEDMD_pool);
    const int& pi2 = index(collide_event->ptr_particle_collide[1], particleEDMD_pool);

    PtrParticleEDMD& particle_collide1 = particleEDMD_pool.at(pi1);
    PtrParticleEDMD& particle_collide2 = particleEDMD_pool.at(pi2);

    // assert(particle_collide1->mobility || particle_collide2->mobility);

    // auto t_prime = [this](const double& t){ 
    //     if (gamma != 0)
    //         return m / gamma * (1 - std::exp(- gamma / m * t));
    //     else
    //         return t;
    // };

    // auto v_decay = [this](const double& t){
    //     if (gamma != 0)
    //         return std::exp(- gamma / m * t);
    //     else
    //         return 1.0;
    // };

    array<double, 3> dr{};
    array<double, 3> dv{};
    double drdv = 0;
    const double& contact_length = 0.5 * (particle_collide1->d + particle_collide2->d);

    for (int i = 0; i < 3; ++ i){
        // Update particles' new coordinate
        particle_collide1->r[i] += particle_collide1->v[i] * (collide_event->t - particle_collide1->t);
        particle_collide2->r[i] += particle_collide2->v[i] * (collide_event->t - particle_collide2->t);

        // Calculate particles' velocity just before collision
        // particle_collide1->v[i] *= v_decay(collide_event->t - particle_collide1->t);
        // particle_collide2->v[i] *= v_decay(collide_event->t - particle_collide2->t);

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

    //Collision with immobile particle, elastic rebounce
    // else if (particle_collide1->mobility)
    //     for (int i = 0; i < 3; ++ i)
    //             particle_collide1->v[i] -= 2.0 * drdv / contact_length * dr[i];
    //     // particle_collide1->velocity(- particle_collide1->v);
    // else if (particle_collide2->mobility)
    //     for (int i = 0; i < 3; ++ i)
    //             particle_collide2->v[i] += 2.0 * drdv / contact_length * dr[i];
    //     // particle_collide2->velocity(- particle_collide2->v);
    // else
    //     assert(false);
    

    // Update particles' time
    particle_collide1->t = collide_event->t;
    particle_collide2->t = collide_event->t;
        
    //Update all Cross/Collide/Rest events for particles
    //Update cross event
    EraseCrossNode(particle_collide1);
    EraseCrossNode(particle_collide2);

    InitializeCrossNode(particle_collide1);
    InitializeCrossNode(particle_collide2);

    //Update rest event
    // EraseRestNode(particle_collide1);
    // EraseRestNode(particle_collide2);

    // InitializeRestNode(particle_collide1);
    // InitializeRestNode(particle_collide2);

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


    // for (int i = -collision_search_range; i <= collision_search_range; ++ i)
    // for (int j = -collision_search_range; j <= collision_search_range; ++ j)
    // for (int k = -collision_search_range; k <= collision_search_range; ++ k){
    //     const Cell& c1 = cell_grid[(particle_collide1->c[0] + i + n) % n]
    //                             [(particle_collide1->c[1] + j + n) % n]
    //                             [(particle_collide1->c[2] + k + n) % n];

    //     //New possible collision of collide particle 1 
    //     for (const int& l : c1.particle_list){
    //         PtrParticleEDMD& p1 = particle_collide1;
    //         PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
    //         if (p2 != particle_collide1 && p2 != particle_collide2)
    //             InitializeCollideNode(p1, p2);    
    //     }
        
    //     const Cell& c2 = cell_grid[(particle_collide2->c[0] + i + n) % n]
    //                             [(particle_collide2->c[1] + j + n) % n]
    //                             [(particle_collide2->c[2] + k + n) % n];
        
    //     //New possible collision of collide particle 2 
    //     for (const int& l : c2.particle_list){
    //         PtrParticleEDMD& p1 = particle_collide2;
    //         PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
    //         if (p2 != particle_collide1 && p2 != particle_collide2)
    //             InitializeCollideNode(p1, p2);
    //     }    
    // }
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
    const double& Delta_barrier = drdv * drdv - 4 * contact_length * contact_length * squarewell_barrier / m;
    const int& sign1 = sgn(drdv);                               //cross in -,cross out +
    const int& sign2 = sgn(Delta_barrier);                      //not cross -, cross +
    // const int& bond_identifier = squarewell_event->bond_status;
    // const int& barrier_identifier = Delta_barrier > 0;

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

    // if(abs(r2-d*d)>1e-8)
    // {
    //     std::cerr<<"Attract_Exe: Numerical Accuracy Loss "<<abs(r2-sig*sig)<<' '<<atr->t;
    //     exit(1);
    // }

    //Square well interaction with additional kinetic energy 
    // const double& scalar_coefficient = (- drdv + bond_identifier * std::sqrt(Delta)) / 2.0 / contact_length / contact_length;
    // const double& scalar_coefficient = (- drdv + std::sqrt(Delta)) / 2.0 / contact_length / contact_length;

    // Calculate particles' new velocities
    for (int i = 0; i < 3; ++ i){
        particle_squarewell1->v[i] += scalar_coefficient * dr[i];
        particle_squarewell2->v[i] -= scalar_coefficient * dr[i];
    }

    //Collision with immobile particle, elastic rebounce
    // else if (particle_collide1->mobility)
    //     for (int i = 0; i < 3; ++ i)
    //             particle_collide1->v[i] -= 2.0 * drdv / contact_length * dr[i];
    //     // particle_collide1->velocity(- particle_collide1->v);
    // else if (particle_collide2->mobility)
    //     for (int i = 0; i < 3; ++ i)
    //             particle_collide2->v[i] += 2.0 * drdv / contact_length * dr[i];
    //     // particle_collide2->velocity(- particle_collide2->v);
    // else
    //     assert(false);
    

    //Update particles' time
    particle_squarewell1->t = squarewell_event->t;
    particle_squarewell2->t = squarewell_event->t;

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

    //Update rest event
    // EraseRestNode(particle_collide1);
    // EraseRestNode(particle_collide2);

    // InitializeRestNode(particle_collide1);
    // InitializeRestNode(particle_collide2);

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
    

    

    // for (int i = -squarewell_search_range; i <= squarewell_search_range; ++ i)
    // for (int j = -squarewell_search_range; j <= squarewell_search_range; ++ j)
    // for (int k = -squarewell_search_range; k <= squarewell_search_range; ++ k){
    //     const Cell& c1 = cell_grid[(particle_squarewell1->c[0] + i + n) % n]
    //                             [(particle_squarewell1->c[1] + j + n) % n]
    //                             [(particle_squarewell1->c[2] + k + n) % n];

    //     //New possible collision of collide particle 1 
    //     for (const int& l : c1.particle_list){
    //         PtrParticleEDMD& p1 = particle_squarewell1;
    //         PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
    //         if (p2 != particle_collide1 && p2 != particle_collide2)
    //             InitializeCollideNode(p1, p2);    
    //     }
        
    //     const Cell& c2 = cell_grid[(particle_collide2->c[0] + i + n) % n]
    //                             [(particle_collide2->c[1] + j + n) % n]
    //                             [(particle_collide2->c[2] + k + n) % n];
        
    //     //New possible collision of collide particle 2 
    //     for (const int& l : c2.particle_list){
    //         PtrParticleEDMD& p1 = particle_collide2;
    //         PtrParticleEDMD& p2 = particleEDMD_pool.at(l);
    //         if (p2 != particle_collide1 && p2 != particle_collide2)
    //             InitializeCollideNode(p1, p2);
    //     }    
    // }

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
}

//Execute reset event
void EDMD::ExecuteResetNode(const Node& node){
    const shared_ptr<ResetEvent>& reset_event = std::dynamic_pointer_cast<ResetEvent>(node->second);
    assert(reset_event != shared_ptr<ResetEvent>());

    // auto t_prime = [this](const double& t){ 
    //     if (gamma != 0)
    //         return m / gamma * (1 - std::exp(- gamma / m * t));
    //     else
    //         return t;
    // };

    // auto v_decay = [this](const double& t){
    //     if (gamma != 0)
    //         return std::exp(- gamma / m * t);
    //     else
    //         return 1.0;
    // };

    //Modify the particle
    for (auto& pi : particleEDMD_pool){
        const double& Delta_t = reset_event->t - pi->t;
        // const double& velocity_decay = std::exp(- gamma / m * Delta_t);

        // Update accumulated time stored in particle cross
        pi->t = 0; 

        for (int i = 0; i < 3; ++ i){
            // Calculate the new particle coordinate
            pi->r[i] += pi->v[i] * Delta_t;
            // Calculate the new particle velocity
            // pi->v[i] *= velocity_decay;
        }            
    }
    
    //Reset the time to decrease numeric error
    if (reset_event->epoch < ResetEvent::max_epoch){
        //Modify the event tree
        for (auto i = event_tree.begin(); i != event_tree.end(); ++ i){
            auto node = event_tree.extract(i);
            node.key() -= ResetEvent::TIME_PER_EPOCH_;
            node.mapped()->t -= ResetEvent::TIME_PER_EPOCH_;
            event_tree.insert(std::move(node));
        }

        EraseResetNode(node);
        ++ reset_event->epoch;
        InitializeResetNode();
    }
    else{
        //Clear the event tree
        event_tree.clear();
        //Close sample data output file
        // SampleEvent::kinetic_termperature_output.close();
    }
}

//Execute sample event
void EDMD::ExecuteSampleNode(const Node& node){
    const shared_ptr<SampleEvent>& sample_event = std::dynamic_pointer_cast<SampleEvent>(node->second);
    assert(sample_event != shared_ptr<SampleEvent>());

    /*-------------------------Sampling Desired Data-------------------------*/
    ++ SampleEvent::sampling_time;

    //Kinetic temperature
    // sample_event->kinetic_temperature = 0;
    // for (auto& pi : particleEDMD_pool){
    //     assert(sample_event->t >= pi->t);
    //     const double& velocity_decay = std::exp(- gamma / m * (sample_event->t - pi->t));
    //     const array<double, 3> sample_velocity(pi->v * velocity_decay);
    //     sample_event->kinetic_temperature += std::accumulate(sample_velocity.begin(), sample_velocity.end(), 0.0, 
    //     [](const double& sum_of_square, const double& i){ return sum_of_square + i * i;});
    // }
    // sample_event->kinetic_temperature /= static_cast<double>(3 * Ntot());

    // sample_event->PrintKineticTemperature();

    //Dump
    assert(std::filesystem::exists(SampleEvent::dump_directory));
    const int& epoch = ResetEvent::GetEpoch(node);
    std::stringstream dump_filename;
    dump_filename << "dump_" << SampleEvent::sampling_time << ".dat";
    // if (sample_event->kinetic_temperature != 0)
        Dump(SampleEvent::dump_directory, dump_filename.str(), epoch, sample_event->t);

    /*-----------------------------------------------------------------------*/
    const double& last_sampling_time = sample_event->t;
    EraseSampleNode(node);
    //Arithmetic sequences sampling time
    // InitializeSampleNode(last_sampling_time + SampleEvent::sampling_interval);
    //Geometric sequences sampling time
    InitializeSampleNode(last_sampling_time + SampleEvent::sampling_interval);
    SampleEvent::sampling_interval *= std::pow(10.0, 0.01);
}


//Execute event tree
void EDMD::ExecuteEventTree(){
    //Check if the tree is empty
    if (event_tree.empty())
        return;
    const Node& node_execute = event_tree.begin();
    const PtrEvent& event_execute = node_execute->second;

    cout << event_execute->t << ' ' << event_execute->Type() << endl;

    if (event_execute->Type() == "cross")
        ExecuteCrossNode(node_execute);
    // else if (event_execute->Type() == "rest")
    //     ExecuteRestNode(node_execute);
    else if (event_execute->Type() == "collide")
        ExecuteCollideNode(node_execute);
    else if (event_execute->Type() == "reset")
        ExecuteResetNode(node_execute);
    else if (event_execute->Type() == "sample")
        ExecuteSampleNode(node_execute);
    else if (event_execute->Type() == "squarewell")
        ExecuteSquareWellNode(node_execute);
    else
        assert(false);
}

void EDMD::ExecuteEDMD(const double& termination_time){
    
    //Initialize event tree
    event_tree.clear();
    //Reset event
    ResetEvent::SetMaxEpoch(static_cast<int>(std::ceil(termination_time / ResetEvent::TIME_PER_EPOCH_)));
    InitializeResetNode();
    //Sample event
    InitializeSampleNode(10);

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

            // for (const int& l : c.particle_list)
            // if (l > index(&pi, particleEDMD_pool)){
            // // if (pi->mobility || particleEDMD_pool.at(l)->mobility){
            //     InitializeCollideNode(pi, particleEDMD_pool.at(l));
            //     // cout << omp_get_thread_num() << endl;
            // }
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
