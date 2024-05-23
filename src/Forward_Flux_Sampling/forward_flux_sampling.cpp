#include <iterator>
#include <list>

#include <omp.h>

#include "forward_flux_sampling.h"

using std::list;

void FFS::InitalizeFirstInterface(const unsigned& experiment_id, const double& bath_temperature, const double& simulation_time, const double& sample_interval){

    path initial_configuration_filename = FFS_folder / "init.dat";

    path initial_interface_dump_folder(FFS_folder / "Lambda0");
    std::filesystem::create_directories(initial_interface_dump_folder);

    path experiment_dump_folder(initial_interface_dump_folder / ("exp" + std::to_string(experiment_id)));
    std::filesystem::create_directories(experiment_dump_folder);

    assert(std::filesystem::exists(initial_configuration_filename));
    assert(std::filesystem::exists(initial_interface_dump_folder));
    assert(std::filesystem::exists(experiment_dump_folder));

    EDMD initial_edmd(InitializeFromDump(initial_configuration_filename.string(), parameters));

    const double& initial_sample_time = sample_interval;
    const double& initial_sample_interval = sample_interval;
    const vector<string>& properties{"bubble"};
    const double& thermostat_frequency = 5 * initial_edmd.Ntot();

    const string& master_log_filename = (initial_interface_dump_folder / "master_log.log").string();
    std::ofstream master_log_out(master_log_filename.c_str(), std::ios::app);
    const string& experiment_log_filename = (experiment_dump_folder / ("exp" + std::to_string(experiment_id) + "_log.log")).string();
    std::ofstream experiment_log_out(experiment_log_filename.c_str(), std::ios::trunc);
    list<double> max_bubble_list;
    int dump_index = 0;

    assert(master_log_out.is_open());
    assert(experiment_log_out.is_open());

    const double& lambda_0 = lambda_vec.front();

    function<void(const double&, EDMD&)> sample_operation = [&](const double& t, EDMD& edmd){

        const double& max_bubble_volume = edmd.MaximumBubbleVolume(t);

        max_bubble_list.emplace_back(max_bubble_volume);

        if (max_bubble_list.back() > lambda_0 && *(++ max_bubble_list.rbegin()) < lambda_0){
            const int& epoch = edmd.GetEpoch();
            const double& time_per_epoch = edmd.GetTimePerEpoch();
            std::stringstream dump_filename;
            dump_filename << "dump_" << dump_index ++ << ".dat";
            // if (sample_event->kinetic_temperature != 0)
            edmd.Dump(experiment_dump_folder, dump_filename.str(), epoch, time_per_epoch, t);
        }
        
        experiment_log_out << t << ' ' << max_bubble_list.back() << ' ' << max_bubble_list.size() << endl;
    };

    function<void(const double&, EDMD&)> terminate_operation = [&](const double& t, EDMD& edmd){
        master_log_out << experiment_id << ' ' << t + edmd.GetEpoch() * edmd.GetTimePerEpoch() << ' ' << dump_index << endl;
        experiment_log_out.close();
        master_log_out.close();
    };

    initial_edmd.SetSamplingParameters("Arithmetic", initial_sample_interval, properties, initial_interface_dump_folder);
    initial_edmd.SetSamplingOperation(sample_operation);
    initial_edmd.SetTerminateOperation(terminate_operation);
    // edmd.InitializeParticleVelocity(temperature * 0.8);
    initial_edmd.ExecuteEDMD(simulation_time, bath_temperature, thermostat_frequency, initial_sample_time);
}

void FFS::ExecuteInterface(const unsigned& layer_id, const unsigned& experiment_id, const unsigned& dump_pool_size, const double& bath_temperature, const double& simulation_time, const double& sample_interval){
    
    // const unsigned& layer_id = 1;
    path input_configuration_folder = FFS_folder / ("Lambda" + std::to_string(layer_id - 1)) / ("exp" + std::to_string(experiment_id));
    assert(std::filesystem::exists(input_configuration_folder));

    path output_configuration_folder = FFS_folder / ("Lambda" + std::to_string(layer_id));
    std::filesystem::create_directories(output_configuration_folder);
    assert(std::filesystem::exists(output_configuration_folder));

    path experiment_dump_folder(output_configuration_folder / ("exp" + std::to_string(experiment_id)));
    std::filesystem::create_directories(experiment_dump_folder);
    assert(std::filesystem::exists(experiment_dump_folder));

    //read input weight
    vector<double> weight_vec(dump_pool_size, 1);
    
    if (layer_id > 1){
        std::ifstream weight_in(input_configuration_folder / ("exp" + std::to_string(experiment_id) + "_log.log"));
        assert(weight_in.is_open());  
        string line;
        while (std::getline(weight_in, line)){
            std::stringstream ss(line);
            unsigned input_index, last_index;
            double w, bv;
            ss >> input_index >> last_index >> w >> bv;
            if (input_index >= 0)
                weight_vec.at(input_index) = w;
        }
        weight_in.close();
    }

    const unsigned& try_times = 100;
    unsigned success_times = 0;
    double success_times_weighted = 0;

    const double& terminate_probability = 0.5;

    const string& master_log_filename = (output_configuration_folder / "master_log.log").string();
    std::ofstream master_log_out(master_log_filename.c_str(), std::ios::app);
    const string& experiment_log_filename = (experiment_dump_folder / ("exp" + std::to_string(experiment_id) + "_log.log")).string();
    std::ofstream experiment_log_out(experiment_log_filename.c_str(), std::ios::trunc);

    assert(master_log_out.is_open());
    assert(experiment_log_out.is_open());

    #pragma omp parallel for num_threads(64) schedule(dynamic)
    for (int i = 0; i < try_times; ++ i){
        
        std::mt19937 gen(std::chrono::high_resolution_clock::now().time_since_epoch().count() + omp_get_thread_num());
        // std::uniform_int_distribution<unsigned> Ui_dice(0, dump_pool_size - 1);
        std::discrete_distribution<int> Uid_dice(weight_vec.begin(), weight_vec.end());
        std::uniform_real_distribution<double> Ur_dice(0,1);

        const unsigned& input_index = Uid_dice(gen);

        path input_configuration_filename = input_configuration_folder / ("dump_" + std::to_string(input_index) + ".dat");
        assert(std::filesystem::exists(input_configuration_filename));

        EDMD layer_edmd(InitializeFromDump(input_configuration_filename.string(), parameters));

        const double& initial_sample_time = sample_interval;
        const double& initial_sample_interval = sample_interval;
        const vector<string>& properties{"bubble"};
        const double& thermostat_frequency = 5 * layer_edmd.Ntot();

        const double& lambda_i = lambda_vec.at(layer_id);
        double weight = 1;

        function<void(const double&, EDMD&)> sample_operation = [&](const double& t, EDMD& edmd){

            const double& max_bubble_volume = edmd.MaximumBubbleVolume(t);

            //if reach current interface, terminate
            if (max_bubble_volume > lambda_i){
                const int& epoch = edmd.GetEpoch();
                const double& time_per_epoch = edmd.GetTimePerEpoch();
                std::stringstream dump_filename;
                dump_filename << "dump_" << success_times << ".dat";
                
                #pragma omp critical
                {
                    experiment_log_out << success_times << ' ' << i << ' ' << weight << ' ' << max_bubble_volume << endl;
                    ++ success_times;
                    success_times_weighted += weight;
                }

                edmd.Dump(experiment_dump_folder, dump_filename.str(), epoch, time_per_epoch, t);
                
                //terminate the run
                edmd.InitializeTerminateNode(std::nextafter(t, std::numeric_limits<double>::max()));
            }

            //Terminate
            if (max_bubble_volume < lambda_A){

                #pragma omp critical
                experiment_log_out << -1 << ' ' << i << ' ' << weight << ' ' << max_bubble_volume << endl;

                edmd.InitializeTerminateNode(std::nextafter(t, std::numeric_limits<double>::max()));
            }
            //Pruning
            for (int j = layer_id - 2; j >= 0 ; -- j)
            if (max_bubble_volume < lambda_vec.at(j))
            //continue
            if (Ur_dice(gen) > terminate_probability)
                weight /= (1 - terminate_probability);
            //terminate
            else{

                #pragma omp critical
                experiment_log_out << -1 << ' ' << i << ' ' << weight << ' ' << max_bubble_volume << endl;

                edmd.InitializeTerminateNode(std::nextafter(t, std::numeric_limits<double>::max()));
            }
        };

        function<void(const double&, EDMD&)> terminate_operation = [&](const double& t, EDMD& edmd){
        };

        layer_edmd.SetSamplingParameters("Arithmetic", initial_sample_interval, properties, output_configuration_folder);
        layer_edmd.SetSamplingOperation(sample_operation);
        layer_edmd.SetTerminateOperation(terminate_operation);
        // edmd.InitializeParticleVelocity(temperature * 0.8);
        layer_edmd.ExecuteEDMD(simulation_time, bath_temperature, thermostat_frequency, initial_sample_time);
    }

    master_log_out << experiment_id << ' ' << try_times << ' ' << success_times_weighted << endl;

    experiment_log_out.close();
    master_log_out.close();
}