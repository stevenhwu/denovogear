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
#ifndef DENOVOGEAR_ASSERT_UTILS_H
#define DENOVOGEAR_ASSERT_UTILS_H

#include <cassert>


const double ASSERT_CLOSE_THRESHOLD = 1e-6;

//TODO(SW): HACK: FIXME:
template<typename A, typename B>
void AssertEqual(A expected, B actual) {
    assert(expected == actual);
};

template<typename A>
void AssertEqual(A expected, A actual) {
    assert(expected == actual);
};

template<typename A>
void AssertNear(A expected, A actual) {
    if (!(expected == 0 && actual == 0)) {
        assert(((expected - actual) / expected) < ASSERT_CLOSE_THRESHOLD);
    }

};

//template<typename A>
//void AssertVectorNear(A expected, A actual){
//    AssertEqual(expected.size(), actual.size());
//    for (int j = 0; j < expected.size(); ++j) {
//        AssertNear(expected[j], actual[j]);
//    }
//};

template<typename A, typename B>
void AssertVectorEqual(A expected, B actual) {
    AssertEqual(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        AssertEqual(expected[j], actual[j]);
    }
};


template<typename A, typename B>
void AssertVectorNear(A expected, B actual) {
    AssertEqual(expected.size(), actual.size());
    for (int j = 0; j < expected.size(); ++j) {
        AssertNear(expected[j], actual[j]);
    }
};


template<typename A>
void AssertEigenMatrixNear(A expected, A actual) {
    AssertEqual(expected.rows(), actual.rows());
    AssertEqual(expected.cols(), actual.cols());
    for (int j = 0; j < expected.rows(); ++j) {
        for (int k = 0; k < expected.cols(); ++k) {
            AssertNear(expected(j, k), actual(j, k));
        }
    }
}


#endif //DENOVOGEAR_ASSERT_UTILS_H
