//
// Created by steven on 1/15/16.
//
// Replace some with with BOOST later
//

#ifndef DENOVOGEAR_BOOST_TEST_HELPER_H
#define DENOVOGEAR_BOOST_TEST_HELPER_H

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

//TODO: Where should this file go?
//FIXME: too many global
const double ABS_TEST_THRESHOLD = 1e-10;
const double BOOST_CLOSE_THRESHOLD = 1e-5;//
//TODO: What is a good cut? 1e-8?// sqrt(std::numeric_limits<double>::epsilon());



template<typename V, typename V2>
void boost_check_vector(V &expected, V2 &result) {

    BOOST_CHECK_EQUAL(expected.size(), result.size());
    for (int i = 0; i < expected.size(); ++i) {
        BOOST_CHECK_CLOSE(expected[i], result[i], BOOST_CLOSE_THRESHOLD);
    }
}

template<typename V, typename V2>
void boost_check_vector(V &expected, V2 &result, int expected_size) {
    
    BOOST_CHECK_EQUAL(expected_size, expected.size());
    boost_check_vector(expected, result);
}



template<typename M>
void boost_check_matrix(M &expected, M &result) {
    BOOST_CHECK_EQUAL(expected.rows(), result.rows());
    BOOST_CHECK_EQUAL(expected.cols(), result.cols());

    for (int i = 0; i < expected.rows(); ++i) {
        for (int j = 0; j < expected.cols(); ++j) {
            BOOST_CHECK_CLOSE(expected(i, j), result(i, j), BOOST_CLOSE_THRESHOLD);
        }
    }
}


template<typename M>
void boost_check_matrix(M &expected, M &result, int expected_rows, int expected_cols) {
    BOOST_CHECK_EQUAL(expected_rows, expected.rows());
    BOOST_CHECK_EQUAL(expected_cols, expected.cols());

    boost_check_matrix(expected, result);
}




#endif //DENOVOGEAR_BOOST_TEST_HELPER_H
