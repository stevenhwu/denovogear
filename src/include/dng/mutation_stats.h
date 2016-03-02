//
// Created by steven on 2/3/16.
//

#ifndef DENOVOGEAR_MUTATION_STATS_H
#define DENOVOGEAR_MUTATION_STATS_H

#include <string>
#include <vector>

#include <dng/hts/bcf.h>

#include <dng/matrix.h>
#include <dng/peeling.h>
#include <dng/pedigree.h>
#include <dng/mutation.h>


class MutationStats {

public:

    MutationStats(double min_prob);


    bool CalculateMutationProb(const dng::peel::workspace_t &work_nomut,
                               const dng::peel::workspace_t &work_full);

    void SetScaledLogLikelihood(double scale);

    void SetGenotypeLikelihoods(const dng::peel::workspace_t &workspace,
                                const int depth_size);

    void SetPosteriorProbabilities(const dng::peel::workspace_t &workspace);


    void CalculateExpectedMutation(dng::peel::workspace_t &work_full,
                                   dng::TransitionVector &mean_matrices);

    void CalculateNodeMutation(dng::peel::workspace_t &work_full,
                               dng::TransitionVector &posmut_transition_matrices);


    void CalculateDenovoMutation(dng::peel::workspace_t &work_nomut,
                                 dng::TransitionVector &onemut_transition_matrices,
                                 const dng::Pedigree &pedigree);



    void SetGenotypeRelatedStats(const int (&acgt_to_refalt_allele)[5],
                                 const int (&refalt_to_acgt_allele)[5],
                                 const uint32_t n_alleles,
                                 const std::size_t ref_index,
                                 const std::size_t num_nodes,
                                 const std::size_t library_start);


    void RecordBasicStats(hts::bcf::Variant &record);

    void RecordGenotypeStats(hts::bcf::Variant &record);

    void RecordSingleMutationStats(hts::bcf::Variant &record);


    enum class OutputLevel {
        Basic, Complete, Singcle
    };
    //TODO: Record different sets of variables,
    // [basic, complete, everytihng/single], something like the following
    //basic [ mup, lld, llh, genotype_likelihood]?
    //complete [basic, mux, posterior, node_mup]
    //single [complete, has_single_mut, mu1p, dnt, dnl, dnq, dnc, node_mu1p]


    void SetNodeMup(const std::vector<double> &event,
                    std::size_t first_nonfounder_index);

    void SetNodeMu1p(std::vector<double> &event, double total,
                     std::size_t first_nonfounder_index);




private://TODO: surely these should not be public. Refactor these while working with call.cc
public:
    float mup_;
    float lld_;
    float llh_;
    float mux_;

    bool has_single_mut_;
    float mu1p_;

    std::string dnt_;
    std::string dnl_;
    int32_t dnq_;
    [[deprecated]] int32_t dnc_;

    dng::IndividualVector posterior_probabilities_;
    dng::IndividualVector genotype_likelihoods_;
    std::vector<float> node_mup_;
    std::vector<float> node_mu1p_;

    std::vector<int32_t> best_genotypes_; //.resize(2 * num_nodes);
    std::vector<int32_t> genotype_qualities_;//(num_nodes);
    std::vector<float> gp_scores_; //(gt_count * num_nodes);
    std::vector<float> gl_scores_; //(gt_count * num_nodes, hts::bcf::float_missing);

private:

    double min_prob_;
    //TODO: Maybe hack these away, using lld, llh
    double logdata_;
    double logdata_nomut_;

    void SetExactlyOneMutation(double total);


    void SetNodeCore(std::vector<float> &stats,
                     const std::vector<double> &event,
                     std::size_t first_nonfounder_index);

};

#endif //DENOVOGEAR_MUTATION_STATS_H
