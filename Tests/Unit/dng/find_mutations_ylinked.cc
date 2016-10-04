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


#define BOOST_TEST_MODULE dng::find_mutation

#include <iostream>
#include <fstream>

#include <dng/task/call.h>
#include <dng/find_mutations_ylinked.h>

#include "../boost_test_helper.h"
#include "fixture_trio_workspace.h"
#include "fixture_read_test_from_file.h"

void ReduceReadDepth(std::vector<depth_t> &read_depths,
        const std::vector<int> &keep_library_index) {
    std::vector<depth_t> read_depths_keep, temp;
    for (auto i : keep_library_index) {
        read_depths_keep.push_back(read_depths[i]);
    }
    temp = read_depths;
    read_depths = read_depths_keep;
    read_depths_keep = temp;

}

struct FixtureFindMutationTrio : public  ReadTrioFromFile {
    std::string fixture;

    dng::RelationshipGraph r_graph;

    FixtureFindMutationTrio(std::string s = "FixtureFindMutationTrio")
            : ReadTrioFromFile(), fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        r_graph.Construct(io_pedigree, rgs, InheritancePattern::Y_LINKED, arg.mu,
                          arg.mu_somatic, arg.mu_library);
        rgs.KeepTheseOnly(r_graph.keep_library_index());
        ReduceReadDepth(read_depths, r_graph.keep_library_index());
    }

    ~FixtureFindMutationTrio() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }
};

struct FixtureFindMutationM12 : public  ReadM12FromFile {

//    double min_prob;
//
    dng::RelationshipGraph r_graph;
////    dng::peel::workspace_t workspace;
//
//    int ref_index = 2;
//    std::vector<depth_t> read_depths{1};
//
//    FindMutations::params_t test_param_1 {0, {{0,0,0,0}}, 0,
//                                          std::string{"0,0,0,0"},
//                                          std::string{"0,0,0,0"} };

    FixtureFindMutationM12(std::string s = "FixtureFindMutationM12") : ReadM12FromFile(), fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        r_graph.Construct(io_pedigree, rgs, InheritancePattern::Y_LINKED, arg.mu,
                          arg.mu_somatic, arg.mu_library);
        rgs.KeepTheseOnly(r_graph.keep_library_index());
        ReduceReadDepth(read_depths, r_graph.keep_library_index());

    }
    ~FixtureFindMutationM12() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }


};
namespace dng {

BOOST_FIXTURE_TEST_CASE(test_constructor, FixtureFindMutationTrio) {


    FindMutationsYLinked find_mutation{min_prob, r_graph, test_param_1};

    BOOST_CHECK_EQUAL(find_mutation.min_prob_, 0.01);
    BOOST_CHECK_EQUAL(find_mutation.params_.theta, 0.001);
    BOOST_CHECK_EQUAL(find_mutation.params_.ref_weight, 1);

    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};
    auto freqs1 = find_mutation.params_.nuc_freq;
    boost_check_close_vector(expect_freqs, freqs1);

    std::vector<std::array<double, 4>> expect_gamma{{0.98, 0.0005, 0.0005, 1.04},
                                                    {0.02, 0.075, 0.005, 1.18}};
    auto gamma_a = find_mutation.params_.params_a;
    BOOST_CHECK_EQUAL(gamma_a.pi, expect_gamma[0][0]);
    BOOST_CHECK_EQUAL(gamma_a.phi, expect_gamma[0][1]);
    BOOST_CHECK_EQUAL(gamma_a.epsilon, expect_gamma[0][2]);
    BOOST_CHECK_EQUAL(gamma_a.omega, expect_gamma[0][3]);

    auto gamma_b = find_mutation.params_.params_b;
    BOOST_CHECK_EQUAL(gamma_b.pi, expect_gamma[1][0]);
    BOOST_CHECK_EQUAL(gamma_b.phi, expect_gamma[1][1]);
    BOOST_CHECK_EQUAL(gamma_b.epsilon, expect_gamma[1][2]);
    BOOST_CHECK_EQUAL(gamma_b.omega, expect_gamma[1][3]);

//#if CALCULATE_ENTROPY == 1
//    auto event = find_mutation.event_;
//    BOOST_CHECK_EQUAL(event.size(), 5);
//    for (int k = 0; k < 5; ++k) {
//        BOOST_CHECK_EQUAL(event[k], 0);
//    }
//#endif
}



BOOST_FIXTURE_TEST_CASE(test_prior, FixtureFindMutationTrio) {

    std::vector<std::array<double, 4>> expected_prior {
        {0.999300699300699, 0.0001998001998002, 0.0001998001998002, 0.0002997002997003},
        {0.0002997002997003, 0.999200799200799, 0.0001998001998002, 0.0002997002997003},
        {0.0002997002997003, 0.0001998001998002, 0.999200799200799, 0.0002997002997003},
        {0.0002997002997003, 0.0001998001998002, 0.0001998001998002, 0.999300699300699},
        {0.3, 0.2, 0.2, 0.3}
    };

    FindMutationsYLinked find_mutation {min_prob, r_graph, test_param_1};
    auto *pArray = find_mutation.genotype_prior_;

    for (int i = 0; i < 5; ++i) {
        auto prior_array = pArray[i];
        boost_check_close_vector(expected_prior[i], prior_array);
    }
}


BOOST_FIXTURE_TEST_CASE(test_full_transition, FixtureFindMutationTrio) {

    std::array<double, 4> freqs{0.3, 0.2, 0.2, 0.3};
    double mu = 1e-8;

    auto orig = f81::matrix(2*mu, freqs);
    auto exp_somatic_full = mitosis_haploid_matrix(orig);
    auto exp_somatic_nomut = mitosis_haploid_matrix(orig, 0);
    auto exp_somatic_posmut = exp_somatic_full - exp_somatic_nomut;
    auto exp_somatic_onemut = mitosis_haploid_matrix(orig, 1);
    auto exp_somatic_mean = mitosis_haploid_mean_matrix(orig);

    std::vector<TransitionMatrix> exp_full {{}, exp_somatic_full};
    std::vector<TransitionMatrix> exp_nomut {{}, exp_somatic_nomut};
    std::vector<TransitionMatrix> exp_posmut {{}, exp_somatic_posmut};
    std::vector<TransitionMatrix> exp_onemut {{}, exp_somatic_onemut};
    std::vector<TransitionMatrix> exp_mean {{}, exp_somatic_mean};

    FindMutationsYLinked find_mutation {min_prob, r_graph, test_param_1};
    auto full_matrices = find_mutation.full_transition_matrices_;
    auto nomut_matrices = find_mutation.nomut_transition_matrices_;
    auto posmut_matrices = find_mutation.posmut_transition_matrices_;
    auto onemut_matrices = find_mutation.onemut_transition_matrices_;
    auto mean_matrices = find_mutation.mean_matrices_;
    for (int i = 0; i < 2; ++i) {
        boost_check_matrix(exp_full[i], full_matrices[i]);
        boost_check_matrix(exp_nomut[i], nomut_matrices[i]);
        boost_check_matrix(exp_posmut[i], posmut_matrices[i]);
        boost_check_matrix(exp_onemut[i], onemut_matrices[i]);
        boost_check_matrix(exp_mean[i], mean_matrices[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_genotype, FixtureFindMutationTrio) {

    std::vector<std::array<double, 4>> expected_genotype {
        {1,1,1,1},
        {2.950305540749871563719e-18, 2.950305540749871563719e-18, 1.000000000000000000000e+00, 2.950305540749871563719e-18}
    };

    FindMutationsYLinked::stats_t mutation_stats;
//    MutationStats stats(min_prob);
    FindMutationsYLinked find_mutation{min_prob, r_graph, test_param_1};
    find_mutation(read_depths, ref_index, &mutation_stats);

    IndividualVector lower = find_mutation.work_full_.lower;

    for (int i = 1; i < 2; ++i) {
        GenotypeArray genotype_array = lower[i];
        boost_check_close_vector(expected_genotype[i], lower[i]);
    }
}


BOOST_FIXTURE_TEST_CASE(test_operator, FixtureFindMutationTrio) {

    std::vector<std::array<double, 4>> expected_nomut_lower {
        {2.95030548493328e-18, 2.95030547695948e-18, 0.999999978378379, 2.95030548493328e-18},
        {2.950305540749871563719e-18, 2.950305540749871563719e-18, 1.000000000000000000000e+00, 2.950305540749871563719e-18}
    };
    std::vector<std::array<double, 4>> expected_full_lower {
        {5.40540532629463e-09, 5.40540532629463e-09, 0.999999978378379, 5.40540532629463e-09},
        {2.950305540749871563719e-18, 2.950305540749871563719e-18, 1.000000000000000000000e+00, 2.950305540749871563719e-18}
    };

    std::vector<std::array<double, 4>> expected_upper {
        {0.0002997002997003, 0.0001998001998002, 0.999200799200799, 0.0002997002997003},
    };

    FindMutationsYLinked::stats_t mutation_stats;
//    MutationStats stats(min_prob);
    FindMutationsYLinked find_mutation{min_prob, r_graph, test_param_1};
    find_mutation(read_depths, ref_index, &mutation_stats);


    IndividualVector full_lower = find_mutation.work_full_.lower;
    IndividualVector nomut_lower = find_mutation.work_nomut_.lower;
    IndividualVector upper = find_mutation.work_full_.upper;

    for (int i = 0; i < 2; ++i) {
        boost_check_close_vector(expected_full_lower[i], full_lower[i]);
        boost_check_close_vector(expected_nomut_lower[i], nomut_lower[i]);
    }
    boost_check_close_vector(expected_upper[0], upper[0]);
    double expected_nomut_ret = -0.000799541952039;
    double expected_full_ret = -0.000799541947716;

    double expected_mup = 4.32344416343e-12;

//    logdata_nomut_ = work_nomut.forward_result;
//    logdata_ = work_full.forward_result;

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
//    const double expected_mup = -std::expm1(result_nomut - result_full);
//    const double expected_lld = (result_full +scale) / M_LN10;
//    const double expected_llh = expected_full_ret / M_LN10;
//
    BOOST_CHECK_CLOSE(expected_mup, mutation_stats.mup, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
//    BOOST_CHECK_CLOSE(expected_lld, mutation_stats.lld, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
//    BOOST_CHECK_CLOSE(expected_llh, mutation_stats.llh, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
//
//    //Test posterior
//    full_matrices = find_mutation.full_transition_matrices_;
//    workspace.CleanupFast();
//    dng::peel::up(workspace, family[0], full_matrices );
//    dng::peel::to_father(workspace, family[1], full_matrices );
//    dng::peel::up(workspace, family[2], full_matrices );
//
//    dng::peel::up_reverse(workspace, family[2], full_matrices );
//    dng::peel::to_father_reverse(workspace, family[1], full_matrices );
//    dng::peel::up_reverse(workspace, family[0], full_matrices );
//
//    dng::GenotypeArray expected_posterior;
//    for (std::size_t i = 0; i < workspace.num_nodes; ++i) {
//        expected_posterior = workspace.lower[i] * workspace.upper[i];
//        expected_posterior /= expected_posterior.sum();
//
//        boost_check_close_vector(expected_posterior,
//        		mutation_stats.posterior_probabilities[i]);
//    }



}
//
//
//BOOST_FIXTURE_TEST_CASE(test_calculate_mutation_expected, TrioWorkspace) {
//
//    FindMutationsYLinked find_mutation{min_prob, r_graph, test_param_1};
//
//    FindMutationsYLinked::stats_t mutation_stats;
////    MutationStats mutation_stats(min_prob);
//    find_mutation(read_depths, ref_index, &mutation_stats);
//
//    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mup, 1e-6);
//    BOOST_CHECK_SMALL(-30.5967 - mutation_stats.lld, 1e-4);
////    BOOST_CHECK_SMALL(-6.68014 - mutation_stats.llh, 1e-5);
//    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mux, 1e-6);
//    BOOST_CHECK_SMALL(0.116189 - mutation_stats.mu1p, 1e-6);
//
//    BOOST_CHECK_EQUAL("GGxGG>GT", mutation_stats.dnt);
//    BOOST_CHECK_EQUAL("LB-NA12878:Solexa-135852", mutation_stats.dnl);
//    BOOST_CHECK_EQUAL(42, mutation_stats.dnq);
//
//#if CALCULATE_ENTROPY == 1
//    BOOST_CHECK_EQUAL(100, mutation_stats.dnc);
//#endif
//
//}


BOOST_FIXTURE_TEST_CASE(test_operator2, FixtureFindMutationM12) {


    FindMutationsYLinked::stats_t mutation_stats;
//    MutationStats stats(min_prob);
    FindMutationsYLinked find_mutation{min_prob, r_graph, test_param_1};
//    find_mutation(read_depths, ref_index, &mutation_stats);


//    IndividualVector full_lower = find_mutation.work_full_.lower;
//    IndividualVector nomut_lower = find_mutation.work_nomut_.lower;
//    IndividualVector upper = find_mutation.work_full_.upper;

    BOOST_CHECK_EQUAL(-10, mutation_stats.mup);


}
} // namespace dng
