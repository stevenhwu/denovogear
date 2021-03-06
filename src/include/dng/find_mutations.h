/*
 * Copyright (c) 2016 Steven H. Wu
 * Copyright (c) 2016 Reed A. Cartwright
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
 *           Reed A. Cartwright <reed@cartwrig.ht>
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
#ifndef DNG_FIND_MUTATIONS_H_
#define DNG_FIND_MUTATIONS_H_


#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include <dng/detail/unit_test.h>
#include <dng/relationship_graph.h>
#include <dng/likelihood.h>
#include <dng/matrix.h>
#include <dng/mutation_stats.h>
#include <dng/peeling.h>



namespace dng {

class FindMutations {
public:
    struct params_t {
        double theta;
        std::array<double, 4> nuc_freq;
        double ref_weight;

        dng::genotype::DirichletMultinomialMixture::params_t params_a;
        dng::genotype::DirichletMultinomialMixture::params_t params_b;
    };

    struct stats_t {
        float mup;
        float lld;
        [[deprecated]]float llh;
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

    FindMutations(double min_prob, const RelationshipGraph &graph,
            params_t params);

    ~FindMutations();

    bool operator()(const std::vector<depth_t> &depths, int ref_index,
            stats_t *stats);


protected:
    const dng::RelationshipGraph relationship_graph_;

    params_t params_;

    [[deprecated]] double min_prob_;

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

    [[deprecated]] std::array<double, 5> max_entropies_;
    [[deprecated]] std::vector<double> event_;


    bool CalculateMutationProb(MutationStats &stats);

    void CalculateDenovoMutation(MutationStats &mutation_stats);

    DNG_UNIT_TEST(test_constructor);
    DNG_UNIT_TEST(test_prior);
    DNG_UNIT_TEST(test_full_transition);
    DNG_UNIT_TEST(test_operator);

};

} // namespace dng

#endif /* DNG_FIND_MUTATIONS_H_ */
