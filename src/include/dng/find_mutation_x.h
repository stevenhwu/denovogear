

#pragma once
#ifndef DENOVOGEAR_FIND_MUTATION_X_H
#define DENOVOGEAR_FIND_MUTATION_X_H

#include <iostream>

#include <dng/hts/bam.h>
#include <dng/hts/bcf.h>
#include <dng/pedigree.h>
#include <dng/likelihood.h>

#include <dng/relationship_graph.h>
#include <dng/mutation_other_patterns.h>

#include "mutation_stats.h"

//using namespace dng::task;
using namespace dng;

class AbstractFindMutations {//TODO: Abstract or Base? decide later

public:
    virtual void SomeFunction() = 0;
    //Hack some common functions, Mostly from the original FindMutations

    struct FindMutationParams {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    bool CalculateMutation(const std::vector<depth_t> &depths, std::size_t ref_index,
                           MutationStats &mutation_stats);


protected:

    const dng::RelationshipGraph &ship_graph;
    FindMutationParams params_;

    dng::peel::workspace_t work_full_;
    dng::peel::workspace_t work_nomut_;

    dng::TransitionVector full_transition_matrices_;
    dng::TransitionVector nomut_transition_matrices_;
    dng::TransitionVector posmut_transition_matrices_;
    dng::TransitionVector onemut_transition_matrices_;
    dng::TransitionVector mean_matrices_;

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;
    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)


    bool CalculateMutationProb(MutationStats &stats);

    void CalculateDenovoMutation(MutationStats &mutation_stats);

    //http://www.urbandictionary.com/define.php?term=shipping
    AbstractFindMutations(const dng::RelationshipGraph &ship_graph,
                          FindMutationParams &params) :
            ship_graph{ship_graph}, params_(params),
            genotype_likelihood_{params.params_a, params.params_b},
            work_nomut_(ship_graph.CreateWorkspace()) {



        // Use a parent-independent mutation model, which produces a
        // beta-binomial
        genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq,
                                              {params_.ref_weight, 0, 0, 0});
        genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq,
                                              {0, params_.ref_weight, 0, 0});
        genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq,
                                              {0, 0, params_.ref_weight, 0});
        genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq,
                                              {0, 0, 0, params_.ref_weight});
        genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

        // Calculate mutation expectation matrices
        full_transition_matrices_.assign(work_nomut_.num_nodes, {});
        nomut_transition_matrices_.assign(work_nomut_.num_nodes, {});
        posmut_transition_matrices_.assign(work_nomut_.num_nodes, {});
        onemut_transition_matrices_.assign(work_nomut_.num_nodes, {});
        mean_matrices_.assign(work_nomut_.num_nodes, {});

    }
};



class FindMutationsXLinked : public AbstractFindMutations {

public:

    void SomeFunction() { };

    FindMutationsXLinked(const RelationshipGraph &ship_graph,
                         FindMutationParams &params);

    void init();

    void run(const std::vector<depth_t> &depths,
             int ref_index, MutationStats &mutation_stats);


};

#endif //DENOVOGEAR_FIND_MUTATION_H
