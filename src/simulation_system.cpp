#include <fstream>

#include "simulation_system.h"

using std::ofstream;

//Calculate the distance vector bewteen two particles
array<double, 3> System::distancevector(const Particle& p1, const Particle& p2, const double& t){
    //Coordinate of p1 and p2
    array<double, 3> r1{};
    array<double, 3> r2{};
    array<double, 3> dr{};
    double d = 0;
    for (int i = 0; i < 3; ++ i){

        r1[i] = static_cast<double>(p1.c[i]) * l + p1.r[i] + p1.v[i] * (t - p1.t);
        r2[i] = static_cast<double>(p2.c[i]) * l + p2.r[i] + p2.v[i] * (t - p2.t);

        dr[i] = r1[i] - r2[i];

        if (2 * dr[i] > L)
            dr[i] = std::fmod(dr[i], L) - L;
        if (2 * dr[i] < -L)
            dr[i] = std::fmod(dr[i], L) + L;
    }
    return dr;
}

//Calulate distance of two particles
double System::distance(const Particle& p1, const Particle& p2, const double& t){
    //square root of sum of squares
    return norm2(distancevector(p1, p2, t));
}

array<double, 2> System::contacttime(const Particle& p1, const Particle& p2, const double& len){
    
    if (p1.v == array<double, 3>{} && p2.v == array<double, 3>{})
        return array<double, 2>{-1, -1};

    const double& dt = p1.t - p2.t;
    array<double, 3> dr{};
    array<double, 3> dv{};
    double drdv = 0;
    double drdr = 0;
    double dvdv = 0;

    for (int i = 0; i < 3; ++ i){
        dr[i] = static_cast<double>(p1.c[i] - p2.c[i]) * l + p1.r[i] - p2.r[i] - p2.v[i] * dt;
        dv[i] = p1.v[i] - p2.v[i];
        if (2 * dr[i] > L)
            dr[i] = std::fmod(dr[i], L) - L;
        if (2 * dr[i] < -L)
            dr[i] = std::fmod(dr[i], L) + L;
        drdv += dr[i] * dv[i];
        drdr += dr[i] * dr[i];
        dvdv += dv[i] * dv[i];
    }

    assert(drdv < 0);
    //Discriminant (Delta)
    double Delta = drdv * drdv - dvdv * (drdr - len * len);

    if (Delta < 0)
        return array<double, 2>{-1, -1};
    
    array<double, 2> roots_Delta_t{
        (- drdv + std::sqrt(Delta)) / dvdv,
        (- drdv - std::sqrt(Delta)) / dvdv
    };

    array<double, 2> t{};

    for (int i = 0; i < 2; ++ i)
        t[i] = roots_Delta_t[i] + p1.t;
    
    return t;
}

void System::SwapParticle(PtrParticle& p1, PtrParticle& p2){
    //Check input validity
    const int& i1 = index(&p1, particle_pool);
    const int& i2 = index(&p2, particle_pool);

    //Check if p1, p2 are the same particle
    if (i1 == i2)
        return;
    
    //Adjust particle list in cell grid
    std::swap(*p1->ptr_cell_particle_list, *p2->ptr_cell_particle_list);
    //Swap particles
    // p1.swap(p2);
    std::swap(p1, p2);
}

bool System::AddParticle(const PtrParticle& particle_add){
    //Check input validity
    assert(particle_data.find(particle_add->type) != particle_data.end());
    
    //Check packing fraction
    const double& MAX_ETA = 0.4;
    if(eta() >= MAX_ETA)
        cout << "WARNNING: PACKING FRACTION IS TOO HIGH" << eta() <<  "(> 0.4)" << endl;

    //Check overlap

    //Largest particle diameter
    const double& d_max = dmax();
    //Size of cell list
    //cell_list_size better to be 1
    const int& cell_list_size = static_cast<int>(ceil(d_max / l));
    for (int i =- cell_list_size; i <= cell_list_size; ++ i)
    for (int j =- cell_list_size; j <= cell_list_size; ++ j)
    for (int k =- cell_list_size; k <= cell_list_size; ++ k){
        //If overlapping, return
        const Cell& c = cell_grid[(particle_add->c[0] + i + n) % n][(particle_add->c[1] + j + n) % n][(particle_add->c[2] + k + n) % n];
        for (const int& l : c.particle_list){
            const Particle& p1 = *particle_add;
            const Particle& p2 = *particle_pool.at(l);
            if (distance(p1, p2, p1.t) < 0.5 * (p1.d + p2.d))
                return false;
        }   
    }

    //Update particle number
    ++ particle_data[particle_add->type].number;

    //Insert particle into particle pool
    particle_pool.emplace_back(particle_add);   

    //Add particle into cell
    AddParticleInCell(particle_pool.back());

    //Swap particle to a random place
    std::uniform_int_distribution<int> Ui_dice(0, particle_pool.size() - 1);
    SwapParticle(particle_pool.back(), particle_pool.at(Ui_dice(gen)));
    
    return true;
}

void System::MoveParticle(PtrParticle& particle_move, const array<int, 3>& new_cell, const array<double, 3>& new_coordinate){

    //Check input validity
    const int& particle_index = index(&particle_move, particle_pool);
    for (int i = 0; i < 3; ++ i)
    assert(new_cell[i] >= 0 && new_cell[i] <= n && new_coordinate[i] >= 0 && new_coordinate[i] <= l);

    //If different cells
    if (particle_move->c != new_cell){
        //Remove particle from the old cell
        RemoveParticleFromCell(particle_move);
        //Update location
        particle_move->location(new_cell, new_coordinate);
        //Add particle to the new cell
        AddParticleInCell(particle_move);
    }
    else
        //Update location
        particle_move->location(new_cell, new_coordinate);
}

void System::RemoveParticle(PtrParticle& particle_remove)
{
    //Check input validity
    const int& particle_index = index(&particle_remove, particle_pool);

    //Update particle number
    assert(particle_data[particle_remove->type].number > 0);
    -- particle_data[particle_remove->type].number;

    //Swap the particle to be removed with the last particle
    SwapParticle(particle_remove, particle_pool.back());

    //Remove particle from cell
    RemoveParticleFromCell(particle_pool.back());

    //Erase the particle
    particle_pool.pop_back();
}

void System::InitializeParticleVelocity(const double& kinetic_temperature){
    const double& T = kinetic_temperature;
    const double& velocity_standard_deviation = std::sqrt(kB * T / mass);
    std::normal_distribution<double> gaussian_generator(0, velocity_standard_deviation);
    for (PtrParticle& pi : particle_pool)
    pi->velocity(array<double, 3>{gaussian_generator(gen), gaussian_generator(gen), gaussian_generator(gen)});
};

void System::Dump(const path& dump_directory, const string& dump_filename, const unsigned& epoch, const double& current_time){

    std::filesystem::create_directories(dump_directory);
    assert(std::filesystem::exists(dump_directory));

    ofstream dump_out((dump_directory / dump_filename).string(), std::ios::trunc);
    assert(dump_out.is_open());

    dump_out << "ITEM: TIMESTEP" << endl;
    dump_out << epoch << ' ' << std::hexfloat << current_time << endl;
    // dump_out << real_time << endl;

    dump_out << "ITEM: NUMBER OF ATOMS" << endl;
    dump_out << Ntot() << endl;

    dump_out << "ITEM: BOX BOUNDS pp pp pp" << endl;
    for(int i = 0; i < 3; ++ i)
    dump_out << "0 " << std::hexfloat << L << endl;
    // dump_out << "0 " << L << endl;

    dump_out << "ITEM: ATOMS id type diameters x y z vx vy vz" << endl;
    for (const auto& particle : particle_pool){
        dump_out << index(&particle, particle_pool) << ' ';
        dump_out << particle->type << ' ';
        dump_out << particle->d << ' ';
        for (int i = 0; i < 3; ++ i)
        dump_out << std::hexfloat << static_cast<double>(particle->b[i]) * L + static_cast<double>(particle->c[i]) * l + particle->r[i] + particle->v[i] * (current_time - particle->t) << ' ';
        // dump_out << (static_cast<double>(particle->c[i]) * l + particle->r[i] + particle->v[i] * t_prime(current_time - particle->t)) << ' ';
        for (int i = 0; i < 3; ++ i)
        dump_out << std::hexfloat << particle->v[i] * (current_time - particle->t) << ' ';
        // dump_out << particle->v[i] * v_decay(current_time - particle->t) << ' ';
        dump_out << endl;
    }

    dump_out.close();
}

bool System::CheckOverlap(){
    for (const PtrParticle& p : particle_pool)
    for (int i = -1; i <= 1; ++ i)
    for (int j = -1; j <= 1; ++ j)
    for (int k = -1; k <= 1; ++ k){
        const Cell& c = cell_grid[(p->c[0] + i + n) % n][(p->c[1] + j + n) % n][(p->c[2] + k + n) % n];
        for (const int& l : c.particle_list)
        if (l > index(&p, particle_pool))
        if (distance(*p, *particle_pool.at(l), p->t > particle_pool.at(l)->t ? p->t : particle_pool.at(l)->t) < 0.5 * (p->d +  particle_pool.at(l)->d)){
            cerr << index(&p, particle_pool) << ' ' << l << endl;
            return false;
        }       
    }
    return true;  
};


