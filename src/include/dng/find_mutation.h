/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *
 * This file is part of DeNovoGear.
 *
 * DeNovoGear is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once
#ifndef DENOVOGEAR_FIND_MUTATION_H
#define DENOVOGEAR_FIND_MUTATION_H

#include <iostream>

#include <dng/hts/bam.h>
#include <dng/hts/bcf.h>
#include <dng/pedigree.h>
#include <dng/likelihood.h>
#include <dng/mutation.h>
#include <dng/mutation_stats.h>
#include <dng/detail/unit_test.h>



//using namespace dng::task;
using namespace dng;

std::vector<std::pair<std::string, uint32_t>> parse_contigs(const bam_hdr_t *hdr);

std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr);


class FindMutations {
public:

    DNG_UNIT_TEST(test_constructor);
    DNG_UNIT_TEST(test_prior);
    DNG_UNIT_TEST(test_full_transition);
    DNG_UNIT_TEST(test_operator);


    struct FindMutationParams {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    FindMutations(double min_prob, const Pedigree &pedigree, FindMutationParams params);


    //TODO(SW): either use this function, or replace operator() with this
    bool CalculateMutation(const std::vector<depth_t> &depths, std::size_t ref_index,
                           MutationStats &mutation_stats);

    //TODO(SW): use mutation_stats.cc
    struct [[deprecated]] stats_t {
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

        IndividualVector posterior_probabilities;
        IndividualVector genotype_likelihoods;
        std::vector<float> node_mup;
        std::vector<float> node_mu1p;
    };

    [[deprecated]] bool old_operator(const std::vector<depth_t> &depths,
                                     int ref_index, stats_t *stats);


protected:

    const dng::Pedigree &pedigree_;
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

    [[deprecated]] double min_prob_;
    [[deprecated]] dng::peel::workspace_t work_;

    [[deprecated]] std::array<double, 5> max_entropies_;
    [[deprecated]] std::vector<double> event_;


    bool CalculateMutationProb(MutationStats &stats);

    void CalculateDenovoMutation(MutationStats &mutation_stats);


};

class AbstractFindMutations{//TODO: Abstract or Base? decide later

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

//    FindMutations(double min_prob, const Pedigree &pedigree, params_t params);

//    bool operator()(const std::vector<depth_t> &depths, int ref_index,
//                    stats_t *stats);

    //TODO: either place with this function, or replace operator() with this
    bool calculate_mutation(const std::vector<depth_t> &depths, int ref_index, MutationStats &mutation_stats);

//protected://TODO: HACK:

    const dng::Pedigree &pedigree_;
    FindMutationParams params_;
//    double min_prob_;
    dng::peel::workspace_t work_;

    dng::TransitionVector full_transition_matrices_;
    dng::TransitionVector nomut_transition_matrices_;
    dng::TransitionVector posmut_transition_matrices_;
    dng::TransitionVector onemut_transition_matrices_;
    dng::TransitionVector mean_matrices_;

    // Model genotype likelihoods as a mixture of two dirichlet multinomials
    // TODO: control these with parameters
    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;

    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
//    std::array<double, 5> max_entropies_;
    std::vector<double> event_;


    AbstractFindMutations(const Pedigree &pedigree,
                          FindMutationParams params) :
            pedigree_{pedigree}, params_(params),
            genotype_likelihood_{params.params_a, params.params_b},
            work_(pedigree.CreateWorkspace()) {



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
        full_transition_matrices_.assign(work_.num_nodes, {});
        nomut_transition_matrices_.assign(work_.num_nodes, {});
        posmut_transition_matrices_.assign(work_.num_nodes, {});
        onemut_transition_matrices_.assign(work_.num_nodes, {});
        mean_matrices_.assign(work_.num_nodes, {});

    }
};

#endif //DENOVOGEAR_FIND_MUTATION_H
