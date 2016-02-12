//
// Created by steven on 2/10/16.
//


#define BOOST_TEST_MODULE lib::mutation

#include <boost/test/unit_test.hpp>
#include <dng/mutation.h>

#include <boost_test_helper.h>
#include <iostream>
#include <gdbm.h>

using namespace dng;
namespace utf = boost::unit_test;

struct CreateMutationMatrix {

    std::string fixture;

    MutationMatrix equal_mutation_matrix;
    MutationMatrix unequal_mutation_matrix;
    std::array<double, 4> equal_freq {0.25, 0.25, 0.25, 0.25};
    std::array<double, 4> unequal_freq {0.1, 0.2, 0.3, 0.4};


    CreateMutationMatrix(std::string s = "") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: CreateMutationMatrix " << fixture);

        equal_mutation_matrix = f81::matrix(1e-6, equal_freq);
        unequal_mutation_matrix = f81::matrix(1e-6, unequal_freq);

    }

    ~CreateMutationMatrix() {
        BOOST_TEST_MESSAGE("tear down fixture: CreateMutationMatrix " << fixture);
    }


};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }


BOOST_FIXTURE_TEST_SUITE(test_mutation_suite, CreateMutationMatrix )


BOOST_AUTO_TEST_CASE(test_f81_equal, *utf::fixture(&setup, &teardown)) {

    auto diag = equal_mutation_matrix.diagonal();
    std::vector<double> expected_diag = {0.999999, 0.999999, 0.999999, 0.999999};
    std::vector<double> expected_off_diag = {3.333331111e-07, 3.333331111e-07, 3.333331111e-07, 3.333331111e-07};

    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = equal_mutation_matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}

BOOST_AUTO_TEST_CASE(test_f81, *utf::fixture(&setup, &teardown)) {

    auto diag = unequal_mutation_matrix.diagonal();
    std::vector<double> expected_diag = {0.9999987143, 0.9999988571, 0.9999990000, 0.9999991429};
    std::vector<double> expected_off_diag = {1.428570408e-07, 2.857140816e-07, 4.285711225e-07, 5.714281633e-07};


    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = unequal_mutation_matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}


BOOST_AUTO_TEST_CASE(test_f81_clazy, *utf::fixture(&setup, &teardown)) {

    auto matrix = f81::matrix(1e-3, {0.01, 0.1, 0.19, 0.7});
    auto diag = matrix.diagonal();
    std::vector<double> expected_diag = {0.9978677587, 0.9980615989, 0.9982554390, 0.9993538663};
    std::vector<double> expected_off_diag = {0.0000215377905, 0.0002153779050, 0.0004092180195, 0.0015076453352};


    boost_check_close_vector(expected_diag, diag);
    for (int i = 0; i < 4; ++i) {
        auto actual = matrix.col(i);
        Eigen::Array4d expected = Eigen::Array4d::Constant(expected_off_diag[i]);
        expected(i) = expected_diag[i];
        boost_check_close_vector(expected, actual);
    }

}



    /* R Code

options(digits=10)

mu6 <- 1e-6
mu3 <- 1e-3
equal_freq <- rep(0.25,4)
unequal_freq <- seq(0.1, 0.4, by=0.1)
crazy_freq <- c(0.01, 0.1, 0.19, 0.7)

f81<- function(mu, freq){
    beta <- 1/(1- sum(freq^2) )
    exp_beta <- exp(-beta*mu)
    diag<- exp_beta+ freq*(1- exp_beta )
    off_diag<- freq*(1-exp_beta)
    return( list(diag=diag, off_diag=off_diag))
}

f81( mu6, equal_freq)
f81( mu6, unequal_freq)
f81( mu3, crazy_freq)
### Result
> f81( mu6, equal_freq)
$diag
[1] 0.999999 0.999999 0.999999 0.999999

$off_diag
[1] 3.333331111e-07 3.333331111e-07 3.333331111e-07 3.333331111e-07

> f81( mu6, unequal_freq)
$diag
[1] 0.9999987143 0.9999988571 0.9999990000 0.9999991429

$off_diag
[1] 1.428570408e-07 2.857140816e-07 4.285711225e-07 5.714281633e-07

> f81( mu3, crazy_freq)
$diag
[1] 0.9978677587 0.9980615989 0.9982554390 0.9993538663

$off_diag
[1] 0.0000215377905 0.0002153779050 0.0004092180195 0.0015076453352


     */



BOOST_AUTO_TEST_SUITE_END()



BOOST_FIXTURE_TEST_SUITE(test_mutation_suite, CreateMutationMatrix )

BOOST_AUTO_TEST_CASE(test_mitosis_haploid_matrix, *utf::fixture(&setup, &teardown)) {


    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            expected_transition_matrix_negative(i, j) = unequal_mutation_matrix(i, j);
        }
    }
    TransitionMatrix expected_copy_unequal_mutation_matrix = unequal_mutation_matrix;


    TransitionMatrix expected_transition_matrix_zero = TransitionMatrix::Zero(4, 4);
    expected_transition_matrix_zero.diagonal() = unequal_mutation_matrix.diagonal();


    TransitionMatrix expected_transition_matrix_one = unequal_mutation_matrix;
    expected_transition_matrix_one.diagonal() = Eigen::Vector4d::Zero(4);



    auto actual_negative = mitosis_haploid_matrix(unequal_mutation_matrix, -1);
    boost_check_matrix(expected_transition_matrix_negative, actual_negative);
    boost_check_matrix(expected_copy_unequal_mutation_matrix, actual_negative);

    auto actual_zero = mitosis_haploid_matrix(unequal_mutation_matrix, 0);
    boost_check_matrix(expected_transition_matrix_zero, actual_zero);

    auto actual_one = mitosis_haploid_matrix(unequal_mutation_matrix, 1);
    boost_check_matrix(expected_transition_matrix_one, actual_one);

}

BOOST_AUTO_TEST_CASE(test_mitosis_diploid_matrix, *utf::fixture(&setup, &teardown)) {

    TransitionMatrix ret = TransitionMatrix::Zero(10, 10);
//    for(int i = 0; i < 10; ++i) {
//        int h = unfolded_diploid_genotypes_upper[i];
//        for(int j = 0; j < 16; ++j) {
//            int k = folded_diploid_genotypes[j];
//            ret(i, k) += (mutype < 0 || mutype == mitotic_diploid_mutation_counts[h][j]) ?
//                         kronecker_product(m, m, h, j) : 0.0;
//        }
//    }
//
    auto kp = kroneckerProduct(unequal_mutation_matrix, unequal_mutation_matrix);

    //TODO: REDO THIS PART!! doesn't test anything here. Calculate these in R
    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(10, 10);
    for(int m = 0; m < 10; ++m) {
        int i = unfolded_diploid_genotypes_upper[m];

    for (int j = 0; j < 16; ++j) {
            int n = folded_diploid_genotypes[j];
            expected_transition_matrix_negative(m, n) += kp.coeff(i, j);
        }
    }

//    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_zero = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_one = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_two = TransitionMatrix::Zero(10, 10);


    auto actual_negative = mitosis_diploid_matrix(unequal_mutation_matrix, -1);
    boost_check_matrix(expected_transition_matrix_negative, actual_negative);

    auto actual_zero = mitosis_diploid_matrix(unequal_mutation_matrix, 0);
//    boost_check_matrix(expected_transition_matrix_zero, actual_zero);

    auto actual_one = mitosis_diploid_matrix(unequal_mutation_matrix, 1);
//    boost_check_matrix(expected_transition_matrix_one, actual_one);

    auto actual_two = mitosis_diploid_matrix(unequal_mutation_matrix, 2);
//    boost_check_matrix(expected_transition_matrix_two, actual_two);


}

BOOST_AUTO_TEST_CASE(test_meiosis_haploid_matrix, *utf::fixture(&setup, &teardown)) {


//    inline TransitionMatrix meiosis_haploid_matrix(const MutationMatrix &m,
//                                                   int mutype = -1) {
//    TransitionMatrix ret{10, 4};
//    auto mm = mitosis_haploid_matrix(m, mutype);
//    for (int i = 0; i < 10; ++i) {
//        int a = folded_diploid_nucleotides[i][0];
//        int b = folded_diploid_nucleotides[i][1];
//        for (int j = 0; j < 4; ++j) {
//            ret(i, j) = 0.5 * (mm(a, j) + mm(b, j));
//        }
//    }


    TransitionMatrix expected_transition_matrix_negative = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_zero = TransitionMatrix::Zero(10, 10);

    TransitionMatrix expected_transition_matrix_one = TransitionMatrix::Zero(10, 10);


    auto actual_negative = meiosis_haploid_matrix(unequal_mutation_matrix, -1);
//    boost_check_matrix(expected_transition_matrix_negative, actual_negative);

    auto actual_zero = meiosis_haploid_matrix(unequal_mutation_matrix, 0);
//    boost_check_matrix(expected_transition_matrix_zero, actual_zero);

    auto actual_one = meiosis_haploid_matrix(unequal_mutation_matrix, 1);
//    boost_check_matrix(expected_transition_matrix_one, actual_one);

}




//        inline TransitionMatrix mitosis_haploid_matrix(const MutationMatrix &m,
//                                               int mutype = -1) {
//}
//
//inline TransitionMatrix mitosis_diploid_matrix(const MutationMatrix &m,
//                                               int mutype = -1) {
//}
//
//inline TransitionMatrix meiosis_haploid_matrix(const MutationMatrix &m,
//                                               int mutype = -1) {
//}
//
//inline TransitionMatrix meiosis_diploid_matrix(const MutationMatrix &mdad,
//                                               const MutationMatrix &mmom, int mutype = -1) {
//}
//
//inline TransitionMatrix mitosis_haploid_mean_matrix(const MutationMatrix &m) {
//}
//
//inline TransitionMatrix mitosis_diploid_mean_matrix(const MutationMatrix &m) {
//}
//
//inline TransitionMatrix meiosis_haploid_mean_matrix(const MutationMatrix &m) {
//}
//
//inline TransitionMatrix meiosis_diploid_mean_matrix(const MutationMatrix &mdad,
//                                                    const MutationMatrix &mmom) {
//}
//
//namespace f81 {
//
//// Construct an F81 mutation matrix
//inline MutationMatrix matrix(double mu, const std::array<double, 4> &nuc_freq) {
//    double beta = 1.0;
//    for(auto d : nuc_freq) {
//        beta -= d * d;
//    }
//    double p = -expm1(-mu / beta);
//
//    MutationMatrix ret;
//    for(int i = 0; i < 4; ++i) {
//        for(int j = 0; j < 4; ++j) {
//            ret(i, j) = nuc_freq[j] * p;
//        }
//        ret(i, i) += 1.0 - p;
//    }
//    return ret;
//}
//
//// Calculate the maximum likelihood estimate of mu given q fraction of observed differences
//inline double estimate(double q, const std::array<double, 4> &nuc_freq) {
//    double beta = 1.0;
//    for(auto d : nuc_freq) {
//        beta -= d * d;
//    }
//    return -beta * log(1.0 - q / beta);
//}
//
//} // namespace f81
//
//
//inline GenotypeArray population_prior(double theta,
//                                           const std::array<double, 4> &nuc_freq,
//                                           const std::array<double, 4> &prior) {
//    double alpha[4] = {
//            theta *nuc_freq[0] + prior[0], theta *nuc_freq[1] + prior[1],
//            theta *nuc_freq[2] + prior[2], theta *nuc_freq[3] + prior[3]
//    };
//    double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];
//    GenotypeArray ret{10};
//    ret <<
//    alpha[0]*(1.0 + alpha[0]) / alpha_sum / (1.0 + alpha_sum), // AA
//            2.0 * alpha[0]*(alpha[1]) / alpha_sum / (1.0 + alpha_sum), // AC
//            2.0 * alpha[0]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum), // AG
//            2.0 * alpha[0]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // AT
//            alpha[1]*(1.0 + alpha[1]) / alpha_sum / (1.0 + alpha_sum), // CC
//            2.0 * alpha[1]*(alpha[2]) / alpha_sum / (1.0 + alpha_sum), // CG
//            2.0 * alpha[1]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // CT
//            alpha[2]*(1.0 + alpha[2]) / alpha_sum / (1.0 + alpha_sum), // GG
//            2.0 * alpha[2]*(alpha[3]) / alpha_sum / (1.0 + alpha_sum), // GT
//            alpha[3]*(1.0 + alpha[3]) / alpha_sum / (1.0 + alpha_sum); // GG
//    return ret;
//}


BOOST_AUTO_TEST_SUITE_END()

