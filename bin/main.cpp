#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <omp.h>

#include "event_driven_molecular_dynamics.h"
#include "initialization.h"
#include "forward_flux_sampling.h"

using std::cout;
using std::endl;
using std::filesystem::path;
using std::stringstream;

int main (int argc, char* argv[]){

    //General parameters
    const double& temperature = std::stod(argv[1]);
    const double& density = std::stod(argv[2]);
    const double& particle_mass = 1;
    const double& particle_diameter = 1;
    const unsigned& particle_number = std::stoul(argv[3]);
    
    //Square well parameters
    const double& squarewell_width = std::stod(argv[4]);
    const double& squarewell_depth = -1;
    const double& barrier_height = std::stod(argv[5]);

    const int& simulation_no = std::stoi(argv[6]);
    const unsigned& input_layer_size = std::stoul(argv[7]);

    //Simulation box
    const double& simulation_box_length = std::cbrt(particle_number / density);
    const unsigned& box_cell_number = std::floor(simulation_box_length / squarewell_width);

    //Output path
    path working_folder = "/home/raop0002/BarrierSquareWell";
    path log_folder = working_folder / "log" / "test";
    path dump_folder = working_folder / "dump" / "test";

    assert(std::filesystem::exists(working_folder));
    std::filesystem::create_directories(log_folder);
    std::filesystem::create_directories(dump_folder);

    assert(std::filesystem::exists(log_folder));
    assert(std::filesystem::exists(dump_folder));

    //Job name
    stringstream job_name_stream;
    job_name_stream << std::setprecision(3) 
                    << "T" << temperature 
                    << "_rho" << density 
                    << "_lambda" << squarewell_width 
                    << "_Eb" << barrier_height 
                    << "_N" << particle_number
                    << "_" << simulation_no;
    const string& job_name = job_name_stream.str();
    
    //Redirect the log output
    // auto log_output = std::freopen((log_folder / (job_name + ".o")).c_str(), "w", stdout);
    // auto log_error = std::freopen((log_folder / (job_name + ".e")).c_str(), "w", stderr);
    // auto log_output = std::freopen((log_folder / ("ffs.o")).c_str(), "w", stdout);
    // auto log_error = std::freopen((log_folder / ("ffs.e")).c_str(), "w", stderr);
    //Dump path
    path job_dump_folder = dump_folder / job_name;
    std::filesystem::create_directories(job_dump_folder);
    
    Unit unit;
    PeriodicalSquareBox box(simulation_box_length, box_cell_number);
    SquareWellCoefficient coef(squarewell_width, squarewell_depth, barrier_height);
    ParticleData pdata("0", particle_mass, particle_diameter, particle_number);

    Parameters env(box, coef);
    env.CreateAtom(pdata);
    
//////////////////////////////////////////////////////////////////////////////////////////

    FFS ffs(working_folder / "FFS", 20, 130, vector<double>{30, 50, 70, 90, 110}, env);
    // ffs.InitalizeFirstInterface(simulation_no, 1, 100, 0.1);
    ffs.ExecuteInterface(2, simulation_no, input_layer_size, 1, 1000, 0.1);

//////////////////////////////////////////////////////////////////////////////////////////

    // // System sys(InitializeRandomHardSphereECMC(env));
    // System sys(InitializeFromDump("/home/raop0002/BarrierSquareWell/dump/T1_rho0.6_lambda1.5_Eb0_N20000/dump_393.dat", coef));
    // // assert(sys.CheckOverlap());
    // cout << sys.eta() << endl;

    // EDMD edmd(sys);

    // const double& initial_sample_time = 0.01;
    // const double& initial_sample_interval = 0.01 * (std::pow(10, 0.01) - 1);
    // const vector<string>& sample_properties{"temperature", "pressure", "time"};
    // const double& bath_temperature = temperature;
    // const double& thermostat_frequency = 5 * particle_number;
    // const double& simulation_time = 0.2;

    // edmd.SetSamplingParameters("Geometric", initial_sample_interval, sample_properties, job_dump_folder);
    // edmd.SetResetParameters(0.1);
    // // edmd.InitializeParticleVelocity(temperature * 0.8);
    // edmd.ExecuteEDMD(simulation_time, bath_temperature, thermostat_frequency, initial_sample_time);

    return 0;
}