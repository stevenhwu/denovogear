//
// Created by steven on 2/3/16.
//
#include <iostream>

#include <dng/peeling.h>

#include <dng/mutation_stats.h>


MutationStats::MutationStats(double min_prob) : min_prob(min_prob) { }


bool MutationStats::set_mutation_prob(const double logdata_nomut, const double logdata){
    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    this->logdata_nomut = logdata_nomut;
    this->logdata = logdata;
    mup = -std::expm1(logdata_nomut - logdata) ;
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
        //TODO: Eigen 3.3 might have log10()
    }
}

void MutationStats::set_posterior_probabilities(const dng::peel::workspace_t &workspace) {

    posterior_probabilities.resize(workspace.num_nodes);
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
        posterior_probabilities[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        posterior_probabilities[i] /= posterior_probabilities[i].sum();
    }

}

void MutationStats::set_exactly_one_mutation(double total){
    mu1p = total * (1.0 - mup);
    has_single_mut = (mu1p / mup) >= min_prob;
}

void MutationStats::set_node_mup(const std::vector<double> &event,
                                 std::size_t first_nonfounder_index) {
    set_node_core(node_mup, event, first_nonfounder_index);
}

void MutationStats::set_node_mu1p(std::vector<double> &event, int total,
                                  std::size_t first_nonfounder_index) {
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        event[i] = event[i] / total;
    }
    set_node_core(node_mu1p, event, first_nonfounder_index);
}

void MutationStats::set_node_core(std::vector<float> &stats, const std::vector<double> &event,
                                  std::size_t first_nonfounder_index) {
    stats.resize(event.size(), hts::bcf::float_missing);
    for(size_t i = first_nonfounder_index; i < event.size(); ++i) {
        stats[i] = static_cast<float>(event[i]);
    }
    //HOWTO: use std::copy with cast??    std::copy(event.begin()+first_nonfounder_index, event.end(), stats.begin() );
}


float MutationStats::get_mutation_prob() const {
    return mup;
}

bool MutationStats::get_has_single_mut() const {
    return has_single_mut;
}

const dng::GenotypeArray &MutationStats::inspect_posterior_at(int index) const {
    return posterior_probabilities[index];
}

const dng::GenotypeArray &MutationStats::inspect_genotype_at(int index) const {
    return genotype_likelihoods[index];
}


void MutationStats::output_to_vcf(hts::bcf::Variant &record){
    //TODO: to implement, hack hack

//
//record.info("MUP", stats.mup);
//record.info("LLD", stats.lld);
//record.info("LLH", stats.llh);
//record.info("MUX", stats.mux);
//record.info("MU1P", stats.mu1p);





};