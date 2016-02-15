//
// Created by steven on 2/3/16.
//

#ifndef DENOVOGEAR_MUTATION_STATS_H
#define DENOVOGEAR_MUTATION_STATS_H

#include <string>
#include <vector>

#include <dng/matrix.h>
#include <dng/peeling.h>
#include <dng/hts/bcf.h>
class MutationStats {


public:
    float get_mup() const {
        return mup;
    }

    const dng::GenotypeArray &inspect_posterior_at(int index) const {
        return posterior_probabilities[index];
    }

    const dng::GenotypeArray &inspect_genotype_at(int index) const {
        return genotype_likelihoods[index];
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


    MutationStats(double min_prob);

    bool set_mup(const double logdata_nomut, const double logdata);

    bool check_mutation_prob_lt_threshold();

    void set_scaled_log_likelihood(double scale);

    void set_genotype_likelihood(const dng::peel::workspace_t &workspace, const int depth_size);

    void set_posterior_probabilities(const dng::peel::workspace_t &workspace);

    void set_node_mup(const std::vector<double> &event, int first_nonfounder_index);

    //TODO: Not implemented yet
    void add_stats();

    void calculate_stats();

    void output_to_vcf(hts::bcf::Variant &record);


};
#endif //DENOVOGEAR_MUTATION_STATS_H
