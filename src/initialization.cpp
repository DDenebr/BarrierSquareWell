#include <fstream>
#include <iterator>
#include <sstream>

#include "initialization.h"

using std::ifstream;
using std::stringstream;

System InitializeRandomHardSphereECMC(const Parameters& env){
    //Create an system
    System system(env);

    //ECMC parameters
    const double& event_chain_length = system.eta() > 0.5 ? system.sigma * 0.5 : system.sigma;
    const double& l = system.l;
    const double& L = system.L;
    const int& n = system.n;
    auto& gen = system.gen;
    std::uniform_int_distribution<int> random_int_generator(0, n - 1);
    std::uniform_real_distribution<double> random_double_generator(0, l);
    std::normal_distribution<double> gaussian_generator;
    auto& cell_grid = system.cell_grid;
    auto& particle_pool = system.particle_pool;
    const double& diameter_max = system.dmax();
    const int& overlap_search_range = static_cast<int>(std::ceil(diameter_max));

    //Insert particles according to particle data in parameters
    for (const auto& [type, pdata] : system.particle_data){
        const string& ptype = pdata.type;
        assert(type == ptype);
        const double& diameter = pdata.diameter;
        const double& mass = pdata.mass;
        const unsigned & number = pdata.number;

        for (int i = 0; i < number; ++ i){
            System::PtrParticle particle(make_shared<System::Particle>(type, diameter, mass, 0));

            vector<System::PtrParticle*> overlap_particle_list;
            overlap_particle_list.reserve(12);

            do{
                overlap_particle_list.clear();
                //A random position
                particle->location(array<int, 3>{random_int_generator(gen), random_int_generator(gen), random_int_generator(gen)},
                                array<double, 3>{random_double_generator(gen), random_double_generator(gen), random_double_generator(gen)});
                //Check overlap condition
                for (int i =- overlap_search_range; i <= overlap_search_range; ++ i)
                for (int j =- overlap_search_range; j <= overlap_search_range; ++ j)
                for (int k =- overlap_search_range; k <= overlap_search_range; ++ k){
                    //If overlapping, return
                    const System::Cell& c = cell_grid[(particle->c[0] + i + n) % n][(particle->c[1] + j + n) % n][(particle->c[2] + k + n) % n];
                    for (const int& l : c.particle_list){
                        const System::Particle& p1 = *particle;
                        const System::Particle& p2 = *particle_pool.at(l);
                        assert(p1.t == 0 && p2.t == 0);
                        if (particle != particle_pool.at(l))
                        if (system.distance(p1, p2, 0) < 0.5 * (p1.d + p2.d))
                            overlap_particle_list.emplace_back(&particle_pool.at(l));
                    }   
                }
            } while (overlap_particle_list.size() > 4);

            //Insert particle into particle pool
            particle_pool.emplace_back(particle);   

            //Add particle into cell
            system.AddParticleInCell(particle_pool.back());
            
            //Eliminate overlap condition
            int ecmc_times = 0;
            int relay_times_sum = 0;
            while (!overlap_particle_list.empty()){
 /*-------------------------------------------ECMC method-------------------------------------------*/                

                const array<double, 3>& move_direction(system.distancevector(**overlap_particle_list.front(), *particle, 0));
                const array<double, 3>& unit_direction(move_direction / norm2(move_direction));

                //Event chain
                map<double, System::PtrParticle*> event_chain;
                map<double, System::PtrParticle*> overlap_event_chain;
                string last_event_type("");
                double eliminate_displacement = -1;
                double current_chain_displacement = 0;
                //Chain termination
                const auto& terminate_chain = std::make_pair(event_chain_length, nullptr);

                //Event chain execute 
                int relay_times = 0;
                const int& max_relay_times = 200;
                System::PtrParticle* pivot = overlap_particle_list.front();
                //Use velocity to represent direction
                (*pivot)->velocity(unit_direction);
                do{
                /*----------------------Sort all possible events----------------------*/
                    //Terminal chain
                    event_chain.emplace(terminate_chain);
                /////////////////////////////////Cross/////////////////////////////////
                    array<int, 2> cross_direction{};              //crs_dir[0], xyz; crs_dir[1], +-1
                    // double t_crs = 0;
                    double cross_displacement = std::numeric_limits<double>::max();
                    //i=0, cross the bottom; i=1, cross the top
                    for (int i = 0; i < 2; ++i)
                    //x,y,z direction
                    for (int j = 0; j < 3; ++j){
                        const double& t = (i * l - (*pivot)->r[j]) / (*pivot)->v[j];
                        if (t > 0 && t < cross_displacement){
                            cross_displacement = t;
                            cross_direction[1] = -1 + 2 * i; 
                            cross_direction[0] = j;
                        }
                    }

                    if (last_event_type != "overlap")
                        event_chain.emplace(cross_displacement + current_chain_displacement, pivot);

                ////////////////////////////////Contact////////////////////////////////
                    if (last_event_type != "overlap"){
                        array<array<int, 3>, 2> search_range;
                        search_range[0].fill(-overlap_search_range);
                        search_range[1].fill(overlap_search_range);

                        //Contact search optimization
                        for (int i = 0; i < 2; ++i)
                        for (int j = 0; j < 3; ++j){
                            //Not search the cell behind particles
                            if (static_cast<double>(sgn((*pivot)->v[j])) * ((*pivot)->r[j] - 0.5 * l) + l * (static_cast<double>(overlap_search_range) - 0.5) > diameter_max)
                                search_range[i][j] = i + (sgn((*pivot)->v[j]) + 1) / 2 - overlap_search_range;
                            //If previous event is cross event, no need to search all of the neighbors
                            if ((*pivot)->r[j] == 0 || (*pivot)->r[j] == l)
                                search_range[i][j] = (1 - 2 * static_cast<int>((*pivot)->r[j] / l)) * overlap_search_range;
                        }
                        
                        for (int i = search_range[0][0]; i <= search_range[1][0]; ++i)
                        for (int j = search_range[0][1]; j <= search_range[1][1]; ++j)
                        for (int k = search_range[0][2]; k <= search_range[1][2]; ++k){
                            // c1 = &c[(pivot.c[0] + i + m) % m][(pivot.c[1] + j + m) % m][(pivot.c[2] + k + m) % m];
                            const System::Cell& c = cell_grid[((*pivot)->c[0] + i + n) % n][((*pivot)->c[1] + j + n) % n][((*pivot)->c[2] + k + n) % n];
                            for (const int &l : c.particle_list)
                            if (particle_pool.at(l) != *pivot){
                                
                                assert((*pivot)->t == 0 && particle_pool.at(l)->t == 0);
                                if (particle_pool.at(l)->v != array<double, 3>{})
                                    assert(false);

                                const array<double, 2>& contact_displacement(system.contacttime(**pivot, *particle_pool.at(l), ((*pivot)->d + particle_pool.at(l)->d) * 0.5));
                                

                                //Not contact
                                if (contact_displacement == array<double, 2>{-1, -1})
                                    continue;

                                //Just contact
                                if (contact_displacement[0] < 1e5 * std::numeric_limits<double>::epsilon() && 
                                    contact_displacement[0] > -1e5 * std::numeric_limits<double>::epsilon())           
                                    continue;
                                
                                //Eliminate Overlap
                                if ((contact_displacement[0] > 0 && contact_displacement[1] < 0) || (contact_displacement[0] < 0 && contact_displacement[1] > 0)){
                                    //Overlap eliminate event
                                    eliminate_displacement = contact_displacement[0] + current_chain_displacement;
                                    event_chain.emplace(eliminate_displacement, pivot);
                                    overlap_event_chain.emplace(eliminate_displacement, pivot != &particle_pool.back() ? pivot : &particle_pool.at(l));
                                    continue;
                                }

                                //Not contact
                                if (contact_displacement[0] < 0 && contact_displacement[1] < 0)
                                    continue;
                                
                                event_chain.emplace(contact_displacement[1] + current_chain_displacement, &particle_pool.at(l));
                            }
                        }    
                    }
                    
                /*----------------------Execute the nearest event----------------------*/

                    auto next = event_chain.begin();
                    const double& next_chain_displacement = next->first;
                    System::PtrParticle* next_pivot = next->second;
                    // Particle* p_next = next->second;
                    array<int, 3> new_c((*pivot)->c);
                    array<double, 3> new_r((*pivot)->r);
                    array<double, 3> new_v((*pivot)->v);

                    last_event_type = "contact";
                    
                    //Calculate new location
                    for (int i = 0; i < 3; ++ i)
                        new_r[i] += (next_chain_displacement - current_chain_displacement) * (*pivot)->v[i];

                    //Cross event, adjust cell
                    if (next_pivot == pivot)
                    if (overlap_event_chain.find(next_chain_displacement) == overlap_event_chain.end()){
                    // if (next_chain_displacement != eliminate_displacement){
                        for (int i = 0; i < 3; ++ i)
                        if (cross_direction[0] == i){
                            new_r[i] = static_cast<double>((1 - cross_direction[1]) / 2) * l;
                            new_c[i] = (new_c[i] + cross_direction[1] + n) % n;
                        }
                        last_event_type = "cross";
                    }
                    //Overlap elimination
                    else{
                        auto& overlap_particle = overlap_event_chain.find(next_chain_displacement)->second;
                        auto itr = std::find(overlap_particle_list.begin(), overlap_particle_list.end(), overlap_particle);
                        assert(itr != overlap_particle_list.end());
                        overlap_particle_list.erase(itr);
                        overlap_event_chain.erase(next_chain_displacement);
                        last_event_type = "overlap";
                    }
                    
                    //Move the particle
                    system.MoveParticle(*pivot, new_c, new_r);

                    //Update chain length
                    current_chain_displacement = next_chain_displacement;

                    //Update the event chain
                    if (next_pivot == pivot)
                        event_chain.erase(next);
                    else{
                        event_chain.clear();
                        overlap_event_chain.clear();
                    }

                    //Nullify velocity
                    (*pivot)->velocity(array<double, 3>{});
                    if (next_pivot != pivot && next_pivot != nullptr){
                        const array<double, 3>& new_direction(system.distancevector(**next_pivot, **pivot, 0));
                        new_v = new_direction / norm2(new_direction);
                    }
                     
                    //Relay to the next particle
                    if (next_pivot != nullptr)
                        pivot = next_pivot;
                    else
                        break;
                    
                    (*pivot)->velocity(new_v);
                    ++ relay_times;
                    assert(relay_times <= max_relay_times);
                } while(true);

                ++ ecmc_times;  
                relay_times_sum += relay_times;
            }
        }   
    }
    return system;
}

System InitializeFromDump(const string& dump_filename, const SquareWellCoefficient& coef){

    ifstream dump_in(dump_filename, std::ios::in);
    assert(dump_in.is_open());

    string line;
    //Read time
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));
    std::getline(dump_in, line);
    // const double& t = std::stod(line);
    std::stringstream epoch_time(line);
    int epoch; string time_per_epoch; string time;
    epoch_time >> epoch >> time_per_epoch >> time;
    // epoch_time >> epoch >> time;
    const double& t = std::stod(time);

    //Read particle number
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));
    std::getline(dump_in, line);
    const int& N = std::stoi(line);

    //Read box size
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen(' '));
    std::getline(dump_in, line);
    const double& L = std::stod(line);
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));

    //Read particle information
    vector<array<string, 9>> particle_inforamtion;
    map<string, pair<double, int>> particle_data;

    particle_inforamtion.reserve(N);
    dump_in.ignore(std::numeric_limits<std::streamsize>::max(), dump_in.widen('\n'));
    while (std::getline(dump_in, line)){
        stringstream ss(line);
        array<string, 9> particle;
        std::copy(std::istream_iterator<string>(ss), std::istream_iterator<string>(), particle.begin());
        particle_inforamtion.emplace_back(particle);

        //Update particle data
        const string& type = particle[1];
        const double& diameter = std::stod(particle[2]);
        if (particle_data.find(type) == particle_data.end())
            particle_data.emplace(type, pair<double, int>{diameter, 1});
        else{
            assert(particle_data.at(type).first == diameter);
            ++ particle_data.at(type).second;
        }
    }

    // const double& diameter_max = std::max_element(particle_data.begin(), particle_data.end(), 
    //     [](const pair<string, pair<double, int>>& largest, const pair<string, pair<double, int>>& first)
    //     {return largest.second.first < first.second.first;})->second.first;

    const int& n = static_cast<int>(std::floor(L / coef.width()));

    System system(Parameters(PeriodicalSquareBox(L, n), coef));
    for (const auto& [type, info] : particle_data)
        system.CreateAtom(type, 1, info.first, info.second);
    
    auto& cell_grid = system.cell_grid;
    auto& particle_pool = system.particle_pool;

    for (const auto& particle_data : particle_inforamtion){
        
        const string& type = particle_data[1];
        const double& d = std::stod(particle_data[2]);
        
        //Create particle
        System::PtrParticle particle(make_shared<System::Particle>(type, d, 1.0, 0.0));

        //Particle location and velocity
        array<double, 3> location{std::stod(particle_data[3]), std::stod(particle_data[4]), std::stod(particle_data[5])};
        const array<double, 3>& velocity{std::stod(particle_data[6]), std::stod(particle_data[7]), std::stod(particle_data[8])};
        array<int, 3> c{};
        array<double, 3> r{};
        for (int i = 0; i < 3; ++ i){
            if (location[i] < 0)
                location[i] = std::fmod(location[i], L) + L;
            if (location[i] > L)
                location[i] = std::fmod(location[i], L);
            assert(location[i] > 0 && location[i] < L);
            r[i] = std::fmod(location[i], system.l);
            c[i] = static_cast<int>(std::round((location[i] - r[i]) / system.l));
        }

        particle->location(c, r);
        particle->velocity(velocity);

        //Insert particle into particle pool
        particle_pool.emplace_back(particle);  

        //Add particle into cell
        system.AddParticleInCell(particle_pool.back());
    }

    return system;
}

