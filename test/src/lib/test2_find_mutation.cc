//
// Created by steven on 2/4/16.
//


#define BOOST_TEST_MODULE dng::lib::find_mutation

#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>

#include <dng/task/call.h>

#include <fstream>
#include <src/helpers/find_mutation_helper.h>

#include "../test_call.h"
#include "boost_test_helper.h"

#include "../../../src/vt/boost_utils.h"

using namespace dng;
namespace utf = boost::unit_test;

const int NUM_TEST = 100;


// TODO: Example of BOOST_DATA_TEST_CASE and BOOST_PARAM_TEST_CASE.
// TODO: Should be able to replace the for loop with these.
// TODO: Might not be able to use fixture.
//
//
//std::vector<int> test_types;//
//
//BOOST_AUTO_TEST_SUITE(suite1,
//  * utf::fixture<Fx>(std::string("FX"))
//  * utf::fixture<Fx>(std::string("FX2")))
//
//  BOOST_AUTO_TEST_CASE(test1, * utf::fixture(&setup, &teardown))
//  {
//    BOOST_TEST_MESSAGE("running test1");
//    BOOST_TEST(true);
//  }
//
//  BOOST_AUTO_TEST_CASE(test2)
//  {
//    BOOST_TEST_MESSAGE("running test2");
//    BOOST_TEST(true);
//  }
//
//BOOST_DATA_TEST_CASE( test_case_arity1, data::xrange(5), my_var )
//{
//    BOOST_TEST_MESSAGE("running data: ");
//    BOOST_TEST((my_var <= 4 && my_var >= 0));
//}
//BOOST_PARAM_TEST_CASE(test_function, params_begin, params_end);
//BOOST_AUTO_TEST_SUITE_END()


//BOOST_AUTO_TEST_SUITE(test_peeling_suite,  * utf::fixture<Fx>(std::string("FX")) )


struct RandomFamily {

    std::string fixture;


    double min_prob;
    dng::Pedigree pedigree;
    dng::peel::workspace_t workspace;

    FindMutations::params_t test_param_1 {0, {0,0,0,0}, 0,
                                          std::string("0,0,0,0"),
                                          std::string("0,0,0,0") };

    typedef dng::task::Call task_type;

    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    RandomFamily(std::string s = "") : fixture(s) {

        po::options_description ext_desc_, int_desc_;
        po::positional_options_description pos_desc_;
        po::variables_map vm_;

        int argc=4;
        char *argv[argc+1];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";
        argv[2] = (char*) "testdata/sample_5_3/ceu.ped"; //"pedFile";
        argv[3] = (char*) "testdata/sample_5_3/test1.vcf"; //test1.bam
//        argv[2] = (char*) "testDataSW/ceu3.ped"; //"pedFile";
//        argv[3] = (char*) "testDataSW/test3.vcf"; //test1.bam

        add_app_args(ext_desc_, static_cast<typename task_type::argument_type &>(arg));
        int_desc_.add_options()
                ("input", po::value< std::vector<std::string> >(&arg.input), "input files")
                ;
        int_desc_.add(ext_desc_);
        pos_desc_.add("input", -1);
        po::store(po::command_line_parser(argc, argv)
                          .options(int_desc_).positional(pos_desc_).run(), vm_);
        po::notify(vm_);

        // Parse pedigree from file
        dng::io::Pedigree ped;
        std::ifstream ped_file(arg.ped);
        ped.Parse(utils::istreambuf_range(ped_file));

        dng::ReadGroups rgs;
        std::vector<hts::File> indata;
        std::vector<hts::bcf::File> bcfdata;
        for (auto &&str : arg.input) {
            indata.emplace_back(str.c_str(), "r");
            if (indata.back().is_open()) {
                continue;
            }
            throw std::runtime_error("unable to open input file '" + str + "'.");
        }
        bcfdata.emplace_back(std::move(indata[0]));
        rgs.ParseSamples(bcfdata[0]);

        pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

        std::array<double, 4> freqs;
        auto f = util::parse_double_list(arg.nuc_freqs, ',', 4);
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);

        test_param_1 = FindMutations::params_t {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]};


        int min_qual = arg.min_basequal;
        min_prob = arg.min_prob;


    }

    double setup_workspace(int ref_index, std::vector<depth_t> &read_depths ){

        workspace.Resize(5);
        workspace.founder_nodes = std::make_pair(0, 2);
        workspace.germline_nodes = std::make_pair(0, 2);
        workspace.somatic_nodes = std::make_pair(2, 2);
        workspace.library_nodes = std::make_pair(2, 5);

        std::array<double, 4> prior {};
        prior.fill(0);
        prior[ref_index] = test_param_1.ref_weight;
        auto genotype_prior_prior = population_prior(test_param_1.theta, test_param_1.nuc_freq, prior);
        workspace.SetFounders(genotype_prior_prior);

        std::vector<std::string> expect_gamma{"0.98, 0.0005, 0.0005, 1.04",
                                              "0.02, 0.075,  0.005,  1.18"};
        dng::genotype::DirichletMultinomialMixture genotype_likelihood_{
                dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[0]},
                dng::genotype::DirichletMultinomialMixture::params_t {expect_gamma[1]}  };
        double scale = 0.0, stemp;
        for (std::size_t u = 0; u < read_depths.size(); ++u) {
            std::tie(workspace.lower[2 + u], stemp) =
                    genotype_likelihood_(read_depths[u], ref_index);
            scale += stemp;
        }
        return scale;
    }

    ~RandomFamily() {
        BOOST_TEST_MESSAGE("tear down fixture " << fixture);
    }


};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }



BOOST_FIXTURE_TEST_SUITE(test_find_mutation_suite, Fx)

BOOST_AUTO_TEST_CASE(test_constructor, *utf::fixture(&setup, &teardown)) {

    FindMutationsGetter find_mutation {min_prob, pedigree, test_param_1};

    BOOST_CHECK_EQUAL(find_mutation.getMin_prob_() , 0.01);

    //getParams_()
    BOOST_CHECK_EQUAL(find_mutation.getParams_().theta, 0.001);//test_param_1.theta);
    BOOST_CHECK_EQUAL(find_mutation.getParams_().ref_weight, 1);//test_param_1.ref_weight);

    std::array<double, 4> expect_freqs{0.3, 0.2, 0.2, 0.3};
    auto freqs1 = find_mutation.getParams_().nuc_freq;
    boost_check_vector(expect_freqs, freqs1, 4);
//    for (int j = 0; j < 4; ++j) {
//        BOOST_CHECK_EQUAL(freqs1[j], expect_freqs[j]);
//    }
    std::vector<std::array<double, 4>> expect_gamma{{0.98, 0.0005, 0.0005, 1.04},
                                                    {0.02, 0.075,  0.005,  1.18}};
    auto gamma_a = find_mutation.getParams_().params_a;
    BOOST_CHECK_EQUAL(gamma_a.pi, expect_gamma[0][0]);
    BOOST_CHECK_EQUAL(gamma_a.phi, expect_gamma[0][1]);
    BOOST_CHECK_EQUAL(gamma_a.epsilon, expect_gamma[0][2]);
    BOOST_CHECK_EQUAL(gamma_a.omega, expect_gamma[0][3]);

    auto gamma_b = find_mutation.getParams_().params_b;
    BOOST_CHECK_EQUAL(gamma_b.pi, expect_gamma[1][0]);
    BOOST_CHECK_EQUAL(gamma_b.phi, expect_gamma[1][1]);
    BOOST_CHECK_EQUAL(gamma_b.epsilon, expect_gamma[1][2]);
    BOOST_CHECK_EQUAL(gamma_b.omega, expect_gamma[1][3]);

    auto event = find_mutation.getEvent_();
    BOOST_CHECK_EQUAL(event.size(), 5);
    for (int k = 0; k < 5; ++k) {
        BOOST_CHECK_EQUAL(event[k], 0);
    }
}



BOOST_AUTO_TEST_CASE(test_prior, *utf::fixture(&setup, &teardown)) {

/* R Code
prior0<- c(0,0,0,0)
freq<- c(0.3,0.2,0.2,0.3)
theta <- 0.001
tempComb<- combn(4,2)
allComb<-cbind(c(1,1), tempComb[,1:3], c(2,2), tempComb[,4:5], c(3,3), tempComb[,6], c(4,4))
result<- vector(mode="list", length=5)
for(i in 1:5){
    prior<- prior0
    if(i<5){
        prior[i]<- prior[i]+1
    }
    alpha <- freq*theta + prior
    alpha_sum <- sum(alpha)
    scale<- alpha_sum * (1+alpha_sum)
    result[[i]]<- apply(allComb, 2, function(x){
        if(x[1]==x[2]){
            return(alpha[x[1]]*(1+alpha[x[1]])/scale)
        }
        else{
            return(2*alpha[x[1]]*alpha[x[2]]/scale)
        }
    })
}

s<- sapply(result, function(x){
    paste(x, sep="", collapse=", ")
})
cat("{", paste(s, collapse="}, \n{"), "}\n" )
*/

    std::vector<std::array<double, 10>> expected_prior {
        {0.998951118846171, 0.000199760259730275, 0.000199760259730275, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 5.99101048876161e-08, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
        {0.000149820194797706, 0.000299610434542968, 5.99101048876161e-08, 8.98651573314242e-08, 0.998801318621409, 0.000199740289695312, 0.000299610434542968, 9.98701448476561e-05, 5.99101048876161e-08, 0.000149820194797706},
        {0.000149820194797706, 5.99101048876161e-08, 0.000299610434542968, 8.98651573314242e-08, 9.98701448476561e-05, 0.000199740289695312, 5.99101048876161e-08, 0.998801318621409, 0.000299610434542968, 0.000149820194797706},
        {0.000149820194797706, 5.99101048876161e-08, 5.99101048876161e-08, 0.000299640389595412, 9.98701448476561e-05, 3.99400699250774e-08, 0.000199760259730275, 9.98701448476561e-05, 0.000199760259730275, 0.998951118846171},
        {0.29979020979021, 0.00011988011988012, 0.00011988011988012, 0.00017982017982018, 0.19984015984016, 7.99200799200799e-05, 0.00011988011988012, 0.19984015984016, 0.00011988011988012, 0.29979020979021 }
    };

    FindMutationsGetter find_mutation {min_prob, pedigree, test_param_1};
    auto *pArray = find_mutation.getGenotype_prior_();

    for (int i = 0; i < 5; ++i) {
        auto prior_array = pArray[i];
        boost_check_vector(expected_prior[i], prior_array);
    }


}

BOOST_AUTO_TEST_CASE(test_full_transition, *utf::fixture(&setup, &teardown)) {

    std::array<double, 4> freqs{0.3, 0.2, 0.2, 0.3};
    double mu = 1e-8;
    auto dad = f81::matrix(3*mu, freqs);
    auto mom = f81::matrix(3*mu, freqs);

    auto exp_germline_full = meiosis_diploid_matrix(dad, mom);
    auto exp_germline_nomut = meiosis_diploid_matrix(dad, mom, 0);
    auto exp_germline_posmut = exp_germline_full - exp_germline_nomut;
    auto exp_germline_onemut = meiosis_diploid_matrix(dad, mom, 1);
    auto exp_germline_mean = meiosis_diploid_mean_matrix(dad, mom);


    auto orig = f81::matrix(2*mu, freqs);
    auto exp_somatic_full = mitosis_diploid_matrix(orig);
    auto exp_somatic_nomut = mitosis_diploid_matrix(orig, 0);
    auto exp_somatic_posmut = exp_somatic_full - exp_somatic_nomut;
    auto exp_somatic_onemut = mitosis_diploid_matrix(orig, 1);
    auto exp_somatic_mean = mitosis_diploid_mean_matrix(orig);


    FindMutationsGetter find_mutation {min_prob, pedigree, test_param_1};
    auto full_matrices = find_mutation.getFull_transition_matrices_();
    auto nomut_matrices = find_mutation.getNomut_transition_matrices_();
    auto posmut_matrices = find_mutation.getPosmut_transition_matrices_();
    auto onemut_matrices = find_mutation.getOnemut_transition_matrices_();
    auto mean_matrices = find_mutation.getMean_matrices_();

    std::vector<TransitionMatrix> exp_full {{},{},exp_germline_full, exp_somatic_full, exp_somatic_full};
    std::vector<TransitionMatrix> exp_nomut {{},{},exp_germline_nomut, exp_somatic_nomut, exp_somatic_nomut};
    std::vector<TransitionMatrix> exp_posmut {{},{},exp_germline_posmut, exp_somatic_posmut, exp_somatic_posmut};
    std::vector<TransitionMatrix> exp_onemut {{},{},exp_germline_onemut, exp_somatic_onemut, exp_somatic_onemut};
    std::vector<TransitionMatrix> exp_mean {{},{},exp_germline_mean, exp_somatic_mean, exp_somatic_mean};

    for (int i = 0; i < 5; ++i) {
        boost_check_matrix(exp_full[i],   full_matrices[i]);
        boost_check_matrix(exp_nomut[i],  nomut_matrices[i]);
        boost_check_matrix(exp_posmut[i], posmut_matrices[i]);
        boost_check_matrix(exp_onemut[i], onemut_matrices[i]);
        boost_check_matrix(exp_mean[i],   mean_matrices[i]);
    }

}


BOOST_AUTO_TEST_CASE(test_operator, *utf::fixture(&setup, &teardown)) {

    int ref_index = 3;
    std::vector<depth_t> read_depths(3);
    for (int j = 0; j < read_depths.size(); ++j) {
        read_depths[j].counts[0] = 20 + j;
        read_depths[j].counts[1] = j*j;
        read_depths[j].counts[2] = 1;
        read_depths[j].counts[3] = 0;
        read_depths[j].counts[j%3] += 30;
    }

    std::vector<peel::family_members_t> family {
            {1,4},
            {0,1,2},
            {0,3}
    };
    min_prob = 0;
    FindMutations::stats_t stats;
    FindMutationsGetter find_mutation{min_prob, pedigree, test_param_1};
    find_mutation(read_depths, ref_index, &stats);


    double scale = setup_workspace(ref_index, read_depths);

    //Test basic stats
    auto numut_matrices = find_mutation.getNomut_transition_matrices_();
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], numut_matrices );
    dng::peel::to_father_fast(workspace, family[1], numut_matrices );
    dng::peel::up(workspace, family[2], numut_matrices );
    double result_nomut = log((workspace.lower[0] * workspace.upper[0]).sum());


    auto full_matrices = find_mutation.getFull_transition_matrices_();
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], full_matrices );
    dng::peel::to_father(workspace, family[1], full_matrices );
    dng::peel::up(workspace, family[2], full_matrices );
    double result_full = log((workspace.lower[0] * workspace.upper[0]).sum());
//    std::cout << "Expected: " << result_nomut << std::endl;
//    std::cout << "Expected full: " << result_full << std::endl;


    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    const double expected_mup = -std::expm1(result_nomut - result_full);
    const double expected_lld = (result_full +scale) / M_LN10;
    const double expected_llh = result_full / M_LN10;

    BOOST_CHECK_CLOSE(expected_mup, stats.mup, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(expected_lld, stats.lld, BOOST_CLOSE_THRESHOLD);
    BOOST_CHECK_CLOSE(expected_llh, stats.llh, BOOST_CLOSE_THRESHOLD);


    std::cout << stats.mup << "\t" << stats.llh << std::endl;
    std::cout << "EXP: " << expected_mup << "\t" << expected_llh << std::endl;

    //Test posterior
    full_matrices = find_mutation.getFull_transition_matrices_();
    workspace.CleanupFast();
    dng::peel::up(workspace, family[0], full_matrices );
    dng::peel::to_father(workspace, family[1], full_matrices );
    dng::peel::up(workspace, family[2], full_matrices );

    dng::peel::up_reverse(workspace, family[2], full_matrices );
    dng::peel::to_father_reverse(workspace, family[1], full_matrices );
    dng::peel::up_reverse(workspace, family[0], full_matrices );

    dng::GenotypeArray expected_posterior;
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
        expected_posterior = WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, i);
        expected_posterior /= expected_posterior.sum();

        boost_check_vector(expected_posterior, stats.posterior_probabilities[i]);
    }



}


BOOST_AUTO_TEST_CASE(test_calculate_mutation, *utf::fixture(&setup, &teardown)) {

    int ref_index = 3;
    std::vector<depth_t> read_depths(3);
    for (int j = 0; j < read_depths.size(); ++j) {
        read_depths[j].counts[0] = 20 + j;
        read_depths[j].counts[1] = j * j;
        read_depths[j].counts[2] = 1;
        read_depths[j].counts[3] = 0;
    }

    min_prob = 0; //Pass everything
    FindMutations::stats_t stats;
    FindMutationsGetter find_mutation{min_prob, pedigree, test_param_1};
    find_mutation(read_depths, ref_index, &stats);

    MutationStats mutation_stats(min_prob);
    find_mutation.calculate_mutation(read_depths, ref_index, mutation_stats);

    BOOST_CHECK_EQUAL(stats.mup, mutation_stats.mup);
    BOOST_CHECK_EQUAL(stats.llh, mutation_stats.llh);
    BOOST_CHECK_EQUAL(stats.lld, mutation_stats.lld);
    BOOST_CHECK_EQUAL(stats.mux, mutation_stats.mux);
    BOOST_CHECK_EQUAL(stats.has_single_mut, mutation_stats.has_single_mut);
    BOOST_CHECK_EQUAL(stats.mu1p, mutation_stats.mu1p);
    BOOST_CHECK_EQUAL(stats.dnt, mutation_stats.dnt);
    BOOST_CHECK_EQUAL(stats.dnl, mutation_stats.dnl);
    BOOST_CHECK_EQUAL(stats.dnc, mutation_stats.dnc);

    for (int i = 0; i < stats.posterior_probabilities.size(); ++i) {
        boost_check_vector(stats.posterior_probabilities[i], mutation_stats.posterior_probabilities[i]);
    }
    for (int i = 2; i < stats.genotype_likelihoods.size(); ++i) {
        boost_check_vector(stats.genotype_likelihoods[i], mutation_stats.genotype_likelihoods[i]);
//        std::cout << i << "\t" << stats.genotype_likelihoods[i][0] << "\t" << mutation_stats.genotype_likelihoods[i][0] << std::endl;
//        std::cout << stats.genotype_likelihoods[i][1] << "\t" << mutation_stats.genotype_likelihoods[i][1] << std::endl;
    }
    for (int i = 2; i < stats.node_mup.size(); ++i) {
        BOOST_CHECK_EQUAL(stats.node_mup[i], mutation_stats.node_mup[i]);
    }

//    std::cout << stats.dnt << "\t" << stats.dnl  << "\t" << stats.dnq<< "\t" << stats.dnc << std::endl;
//    for (auto p : stats.posterior_probabilities) {
//        std::cout << p.prod() << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : stats.genotype_likelihoods) {
//        std::cout << p.prod() << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : stats.node_mup) {
//        std::cout << p << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : stats.node_mu1p) {
//        std::cout << p << "\t";
//    }
//    std::cout << "\n========================" << std::endl;

//    IndividualVector posterior_probabilities;
//    IndividualVector genotype_likelihoods;
//    std::vector<float> node_mup;
//    std::vector<float> node_mu1p;


//    MutationStats mutation_stats(find_mutation.min_prob_);
//    find_mutation.calculate_mutation(read_depths, ref_index, mutation_stats);
//    std::cout << "==== call calculate_mutation() ======" << std::endl;
//    std::cout << mutation_stats.mup << "\t" << mutation_stats.llh << "\t" << mutation_stats.lld <<
//    "\t" << mutation_stats.mux << "\t?:" << mutation_stats.has_single_mut << "\t" <<
//    mutation_stats.mu1p << std::endl;
//    std::cout << mutation_stats.dnt << "\t" << mutation_stats.dnl << "\t" <<
//    mutation_stats.dnq << "\t" << mutation_stats.dnc << std::endl;
//    for (auto p : mutation_stats.posterior_probabilities) {
//        std::cout << p.prod() << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : mutation_stats.genotype_likelihoods) {
//        std::cout << p.prod() << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : mutation_stats.node_mup) {
//        std::cout << p << "\t";
//    }
//    std::cout << std::endl;
//    for (auto p : mutation_stats.node_mu1p) {
//        std::cout << p << "\t";
//    }
//    std::cout << "\n========================" << std::endl;

}




BOOST_AUTO_TEST_SUITE_END()

