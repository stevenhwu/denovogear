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

    MutationStats(double min_prob);

    float get_mutation_prob() const;

    bool get_has_single_mut() const;

    const dng::GenotypeArray &inspect_posterior_at(int index) const;

    const dng::GenotypeArray &inspect_genotype_at(int index) const;


    bool set_mutation_prob(const double logdata_nomut, const double logdata);

    void set_scaled_log_likelihood(double scale);

    void set_genotype_likelihood(const dng::peel::workspace_t &workspace, const int depth_size);

    void set_posterior_probabilities(const dng::peel::workspace_t &workspace);

    void set_exactly_one_mutation(double total);

    void set_node_mup(const std::vector<double> &event, std::size_t first_nonfounder_index);
    void set_node_mu1p(std::vector<double> &event, int total, std::size_t first_nonfounder_index);


    //TODO: Not implemented yet
    void output_to_vcf(hts::bcf::Variant &record);


//private://TODO: surely these should not be public. Refactor these while working with call.cc
    //TODO: float vs double?
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

    enum class output_levels {BASIC, COMPLETE, EVERYTHING};
    //TODO: Record different sets of variables, [basic, complete, everytihng], something like the following
    //basic [ mup, lld, llh, genotype_likelihood]?
    //complete [basci, mux, posterior, node_mup]
    //everything [complete, has_single_mut, mu1p, dnt, dnl, dnq, dnc, node_mu1p]


private:

    double min_prob;
    double logdata;//TODO: Maybe hack these away, using lld, llh
    double logdata_nomut;

    void set_node_core(std::vector<float> &stats, const std::vector<double> &event,
                       std::size_t first_nonfounder_index);

};
#endif //DENOVOGEAR_MUTATION_STATS_H
