//
// Created by steven on 2/3/16.
//
#include <iostream>

#include <dng/peeling.h>

#include "mutation_stats.h"

MutationStats::MutationStats(double min_prob) : min_prob(min_prob) {

}



double MutationStats::store_mup(const double logdata_nomut, const double logdata){
    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]

    this->logdata_nomut = logdata_nomut;
    this->logdata = logdata;
    mup = -std::expm1(logdata_nomut - logdata) ;//+ 10000;
    return mup;
}

bool MutationStats::check_mutation_prob_lt_threshold() {//overkill?
    return mup < min_prob;
}

void MutationStats::store_scaled_log_likelihood(double scale) {
    lld = (logdata + scale) / M_LN10;
    llh = logdata / M_LN10;
}

void MutationStats::store_genotype_likelihood(const std::vector<dng::depth_t> &depths,
                                              const dng::peel::workspace_t &workspace) {

    genotype_likelihoods.resize(workspace.num_nodes);
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::size_t pos = workspace.library_nodes.first + u;
        genotype_likelihoods[pos] = workspace.lower[pos].log() / M_LN10;
    }
}

void MutationStats::store_posterior_probabilities(const dng::peel::workspace_t &workspace) {


    posterior_probabilities.resize(workspace.num_nodes);
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
//        posterior_probabilities[i] = workspace.upper[i] * workspace.lower[i];
        posterior_probabilities[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        posterior_probabilities[i] /= posterior_probabilities[i].sum();
        std::cout << posterior_probabilities[i].sum() << "\t";
    }
    std::cout << ""<< std::endl;

}
