//
// Created by steven on 2/3/16.
//
#include <iostream>

#include <dng/peeling.h>

#include <dng/mutation_stats.h>
#include <dng/hts/bcf.h>

MutationStats::MutationStats(double min_prob) : min_prob(min_prob) {

}



void MutationStats::set_mup(const double logdata_nomut, const double logdata){
    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    this->logdata_nomut = logdata_nomut;
    this->logdata = logdata;
    mup = -std::expm1(logdata_nomut - logdata) ;

}

bool MutationStats::check_mutation_prob_lt_threshold() {//overkill?
    return mup < min_prob;
}

void MutationStats::set_scaled_log_likelihood(double scale) {
    lld = (logdata + scale) / M_LN10;
    llh = logdata / M_LN10;
}

void MutationStats::set_genotype_likelihood(const dng::peel::workspace_t &workspace, const int depth_size) {

    genotype_likelihoods.resize(workspace.num_nodes);
    for(std::size_t u = 0; u < depth_size; ++u) {
        std::size_t pos = workspace.library_nodes.first + u;
        genotype_likelihoods[pos] = workspace.lower[pos].log() / M_LN10;
    }
}

void MutationStats::set_posterior_probabilities(const dng::peel::workspace_t &workspace) {


    posterior_probabilities.resize(workspace.num_nodes);
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
//        posterior_probabilities[i] = workspace.upper[i] * workspace.lower[i];
        posterior_probabilities[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        posterior_probabilities[i] /= posterior_probabilities[i].sum();
//        std::cout << posterior_probabilities[i].sum() << "\t";
    }
//    std::cout << ""<< std::endl;

}

void MutationStats::set_node_mup(const std::vector<double> &event, int start_index) {
    node_mup.resize(event.size(), hts::bcf::float_missing);
    for(size_t i = start_index; i < event.size(); ++i) {
        node_mup[i] = static_cast<float>(event[i]);
    }
}

