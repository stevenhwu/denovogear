//
// Created by steven on 2/3/16.
//

#ifndef DENOVOGEAR_MUTATION_STATS_H
#define DENOVOGEAR_MUTATION_STATS_H

#include <string>
#include <vector>

#include <dng/matrix.h>

class MutationStats {


public:
    float get_mup() const {
        return mup;
    }

//private://TODO: surely it's not public
    float mup;
    float lld;
    float llh;
    float mux;

    bool has_single_mut;
    float mu1p;

    std::string dnt;
    std::string dnl;
    int32_t dnq;
    int32_t dnc;

    dng::IndividualVector posterior_probabilities;
    dng::IndividualVector genotype_likelihoods;
    std::vector<float> node_mup;
    std::vector<float> node_mu1p;

//internal?!
private:
    double min_prob;
    double logdata;
    double logdata_nomut;


public:
    double store_mup(const double logdata_nomut, const double logdata);


    MutationStats(double min_prob);

    bool check_mutation_prob_lt_threshold();

    void store_scaled_log_likelihood(double scale);

    void store_genotype_likelihood(const std::vector<dng::depth_t> &depths, const dng::peel::workspace_t &workspace);

    void store_posterior_probabilities(const dng::peel::workspace_t &workspace);
};


#endif //DENOVOGEAR_MUTATION_STATS_H
