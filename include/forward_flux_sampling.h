#ifndef FORWARD_FLUX_SAMPLING_INCLUDED
#define FORWARD_FLUX_SAMPLING_INCLUDED

#include "event_driven_molecular_dynamics.h"
#include "event.h"
#include "initialization.h"

class FFS
{
    public:

    FFS(const path& ffs_folder, const double& lambda_A, const double& lambda_B, const vector<double>& lambdas, const Parameters& para) :
        FFS_folder(ffs_folder), lambda_A(lambda_A), lambda_B(lambda_B), lambda_vec(lambdas), parameters(para) {
        std::sort(lambda_vec.begin(), lambda_vec.end());
        assert(lambda_vec.front() > lambda_A && lambda_vec.back() < lambda_B);
    };

    void InitalizeFirstInterface(const unsigned& experiment_id, const double& bath_temperature, const double& simulation_time, const double& sample_interval);
    void ExecuteInterface(const unsigned& layer_id, const unsigned& experiment_id, const unsigned& dump_pool_size, const double& bath_temperature, const double& simulation_time, const double& sample_interval);

    private:
    Parameters parameters;
    path FFS_folder;
    double lambda_A;
    double lambda_B;
    vector<double> lambda_vec;
    
};

#endif