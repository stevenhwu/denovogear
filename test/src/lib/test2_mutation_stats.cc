//
// Created by steven on 2/9/16.
//


#define BOOST_TEST_MODULE dng::lib::mutation_stats

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <dng/mutation_stats.h>

#include "boost_test_helper.h"
#include <fixture/fixture_random_family.h>

using namespace dng;
namespace utf = boost::unit_test;

const int NUM_TEST = 100;


struct Fixture {

    std::string fixture;
    std::uniform_real_distribution<double> rand_unif_log;
    std::uniform_real_distribution<double> rand_unif;

    double min_prob = 0.01;
    Fixture(std::string s = "") : fixture(s) {

        rand_unif_log = std::uniform_real_distribution<double>(std::log(1e-20), 0);
        rand_unif = std::uniform_real_distribution<double>(0,1);

        BOOST_TEST_MESSAGE("set up fixture: Fixture " << s);
    }

    ~Fixture() {
        BOOST_TEST_MESSAGE("tear down fixture: Fixture " << fixture);
    }

};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }



BOOST_FIXTURE_TEST_SUITE(test_mutation_stats_suite, Fixture )

BOOST_AUTO_TEST_CASE(test_set_mup, *utf::fixture(&setup, &teardown)) {

    MutationStats stats (min_prob);

    for (int t = 0; t <NUM_TEST; ++t) {

        double ln_mut = rand_unif_log(random_gen_mt);
        double ln_no_mut = rand_unif_log(random_gen_mt);

        double prob_mu = std::exp(ln_mut);
        double prob_no_mut = std::exp(ln_no_mut);

        double expected = 1 - (prob_no_mut / prob_mu);
        double expected_log =  - std::expm1(ln_no_mut - ln_mut);

        stats.set_mutation_prob(ln_no_mut, ln_mut);
        BOOST_CHECK_CLOSE(expected, expected_log, 0.01);
        BOOST_CHECK_CLOSE(expected_log, stats.get_mutation_prob(), BOOST_CLOSE_THRESHOLD);
    }

    for (int t = 0; t <NUM_TEST; ++t) {

        double prob_mu = rand_unif(random_gen_mt);
        double prob_no_mut = rand_unif(random_gen_mt);

        double ln_mu = std::log(prob_mu);
        double ln_no_mut = std::log(prob_no_mut);

        double expected = 1 - (prob_no_mut / prob_mu);

        stats.set_mutation_prob(ln_no_mut, ln_mu);
        BOOST_CHECK_CLOSE(expected, stats.get_mutation_prob(), BOOST_CLOSE_THRESHOLD);
    }

}





//void MutationStats::set_posterior_probabilities(const dng::peel::workspace_t &workspace) {
//
//
//    posterior_probabilities.resize(workspace.num_nodes);
//    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
////        posterior_probabilities[i] = workspace.upper[i] * workspace.lower[i];
//        posterior_probabilities[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
//        posterior_probabilities[i] /= posterior_probabilities[i].sum();
////        std::cout << posterior_probabilities[i].sum() << "\t";
//    }
////    std::cout << ""<< std::endl;
//
//}

BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(test_mutation_stats_suite2, RandomWorkspace )

BOOST_AUTO_TEST_CASE(test_set_posterior_probabilities, *utf::fixture(&setup, &teardown)) {

    init_workspace();

    dng::IndividualVector expected_probs;
    expected_probs.reserve(workspace.num_nodes);

    for (int i = 0; i < workspace.num_nodes; ++i) {
        double sum = 0;
        for (int j = 0; j < 10; ++j) {
            expected_probs[i][j] = workspace.upper[i][j] * workspace.lower[i][j];
            sum += expected_probs[i][j];
        }
        expected_probs[i] /= sum;
    }

    MutationStats stats (0.1);
    stats.set_posterior_probabilities(workspace);
    for (int i = 0; i < workspace.num_nodes; ++i) {
        boost_check_close_vector(expected_probs[i], stats.inspect_posterior_at(i));

    }

}

BOOST_AUTO_TEST_SUITE_END()
