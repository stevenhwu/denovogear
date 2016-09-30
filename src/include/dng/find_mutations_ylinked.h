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
#ifndef DNG_FIND_MUTATIONS_YLINKED_H_
#define DNG_FIND_MUTATIONS_YLINKED_H_


#include <cstdlib>
#include <vector>

//#include <iostream>
//#include <iomanip>
//#include <ctime>
//#include <chrono>
#include <sstream>
#include <string>

//#include <boost/range/algorithm/replace.hpp>
//#include <boost/range/algorithm/max_element.hpp>

//#include <boost/algorithm/string.hpp>
#include <dng/matrix.h>
#include <dng/task/call.h>
#include <dng/relationship_graph.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>
#include <dng/io/utility.h>
#include <dng/mutation_stats.h>
#include <dng/find_mutations_abstract.h>

namespace dng {


class FindMutationsYLinked : public FindMutationsAbstract {


public:


    FindMutationsYLinked(double min_prob, const RelationshipGraph &graph,
            params_t params);

    ~FindMutationsYLinked();

    bool operator()(const std::vector<depth_t> &depths, int ref_index,
                    stats_t *stats);

//    bool InheritanceY(const std::vector<depth_t> &depths, int ref_index,
//                    stats_t *stats);


protected:
    void SetupTransitionMatrix();

//    const dng::RelationshipGraph relationship_graph_;
//
//    params_t params_;
//
//    double min_prob_;
//
//    dng::peel::workspace_t work_full_;
//    dng::peel::workspace_t work_nomut_;
//
//
//    dng::TransitionVector full_transition_matrices_;
//    dng::TransitionVector nomut_transition_matrices_;
//    dng::TransitionVector posmut_transition_matrices_;
//    dng::TransitionVector onemut_transition_matrices_;
//    dng::TransitionVector mean_matrices_;
//
//    // Model genotype likelihoods as a mixture of two dirichlet multinomials
//    // TODO: control these with parameters
//    dng::genotype::DirichletMultinomialMixture genotype_likelihood_;
//
//    dng::GenotypeArray genotype_prior_[5]; // Holds P(G | theta)
//
//    std::vector<int> keep_library_index_;
//
//    bool CalculateMutationProb(MutationStats &mutation_stats);
//
//    void CalculateDenovoMutation(MutationStats &mutation_stats);
//
//
//    std::array<double, 5> max_entropies_;
//    std::vector<double> event_;
//
//    DNG_UNIT_TEST(test_constructor);
//    DNG_UNIT_TEST(test_prior);
//    DNG_UNIT_TEST(test_full_transition);
//    DNG_UNIT_TEST(test_operator);


};
} // namespace dng

#endif /* DNG_FIND_MUTATIONS_YLINKED_H_ */