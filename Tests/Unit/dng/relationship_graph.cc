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


#define BOOST_TEST_MODULE dng::relationship_graph

#include <iostream>

#include <dng/relationship_graph.h>

#include "../boost_test_helper.h"
#include <boost/test/unit_test.hpp>
#include <test/unit_test.hpp>
#include "fixture_read_trio_from_file.h"
#include <boost/test/unit_test_suite.hpp>
namespace utf = boost::unit_test;


struct FixturePedigree : public ReadTrioFromFile{


    dng::RelationshipGraph pedigree;
//    dng::PedigreeV2 pedigree_v2;

    FixturePedigree(std::string s = "FixturePedigree") : ReadTrioFromFile(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
//        pedigree_v2.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

    }

    ~FixturePedigree() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};

void setup() { BOOST_TEST_MESSAGE("set up fun"); }

void teardown() { BOOST_TEST_MESSAGE("tear down fun"); }

//BOOST_FIXTURE_TEST_SUITE(test_pedigree_suite, FixturePedigree )
BOOST_FIXTURE_TEST_CASE_NO_DECOR(aoeu, FixturePedigree){

}


/*

0     1
|     |
---|---
|  |  |
3  2  4

*/
BOOST_FIXTURE_TEST_CASE_NO_DECOR(test_constructor, FixturePedigree ) {

    BOOST_CHECK_EQUAL(5, pedigree.num_nodes() );

    auto workspace = pedigree.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(2, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(2, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(5, workspace.library_nodes.second);


    auto labels = pedigree.labels();

    const std::vector<std::string> expected_labels = {
        "GL-1", // founder 1
        "GL-2", // founder 2
        "LB-NA12878:Solexa-135852",  // lib 1
        "LB-NA12891:Solexa-135851",  // lib 2
        "LB-NA12892:Solexa-135853"   // lib 3

//#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878:Solexa-135852   NA12891:Solexa-135851   NA12892:Solexa-135853
    };
    for (int j = 0; j < 5; ++j) {
        BOOST_CHECK_EQUAL(expected_labels[j], labels[j]);
    }
/* Code relate to labels
//#define DNG_GL_PREFIX "GL-"
//#define DNG_SM_PREFIX "SM-" // define also in newick.cc
//#define DNG_LB_PREFIX "LB-"
//    // Add the labels for the germline nodes
//    labels[0] = DNG_GL_PREFIX "unknown";
//    for (size_t i = 1; i < num_members; ++i) {
//        labels[i] = DNG_GL_PREFIX + pedigree.name(i);
//    }
//    for(auto && a : rgs.libraries()) {
//        vertex_t v = add_vertex(pedigree_graph);
//        labels[v] = DNG_LB_PREFIX + a;
//    }
//
//    if(!labels[u].empty()) {
//        labels_.push_back(labels[u]);
//    } else {
//        labels_.push_back(DNG_SM_PREFIX "unnamed_node_" + util::to_pretty(vid));
//    }
*/

    auto transitions = pedigree.transitions();
    auto size_t_negative_one = static_cast<size_t>(-1);
    std::vector<RelationshipGraph::transition_t> expected_transitions = {
            { RelationshipGraph::TransitionType::Founder, size_t_negative_one,
					size_t_negative_one, 0, 0 },
            { RelationshipGraph::TransitionType::Founder, size_t_negative_one,
					size_t_negative_one, 0, 0 },
            { RelationshipGraph::TransitionType::Germline, 0, 1, arg.mu
					+ arg.mu_somatic + arg.mu_library, arg.mu
					+ arg.mu_somatic + arg.mu_library },
            { RelationshipGraph::TransitionType::Somatic, 0, size_t_negative_one,
					arg.mu_somatic + arg.mu_library, 0 },
			{ RelationshipGraph::TransitionType::Somatic, 1, size_t_negative_one,
					arg.mu_somatic + arg.mu_library, 0 }
    };
    //Transition related code
//    arg.mu, arg.mu_somatic, arg.mu_library);
//    if(edge_types[*ei] == EdgeType::Meiotic) {
//        lengths[*ei] *= mu;
//    } else if(edge_types[*ei] == EdgeType::Mitotic) {
//        lengths[*ei] *= mu_somatic;
//    } else if(edge_types[*ei] == EdgeType::Library) {
//        lengths[*ei] *= mu_library;
//    }
//
//    for(std::size_t i = first_founder_; i < first_nonfounder_; ++i) {
//        transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 0};
//    }
//    transitions_[child] = {tt, parent, static_cast<size_t>(-1), lengths[*pos], 0};
//    transitions_[child] = { TransitionType::Germline, dad, mom, lengths[*pos], lengths[*(pos + 1)]  };

    for (int k = 0; k < 5; ++k) {
        auto expected = expected_transitions[k];
        auto actual = transitions[k];
        BOOST_CHECK(expected.type == actual.type);
        BOOST_CHECK_EQUAL(expected.parent1, actual.parent1);
        BOOST_CHECK_EQUAL(expected.parent2, actual.parent2);
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

    }


}

//BOOST_AUTO_TEST_CASE(test_pedigree_equal) {
//
//    bool is_equal = pedigree.equal(pedigree_v2);
//    BOOST_CHECK(is_equal);
//}


//BOOST_AUTO_TEST_SUITE_END()

namespace dng {
BOOST_FIXTURE_TEST_CASE_NO_DECOR(test_pedigree_inspect, FixturePedigree) {

    BOOST_CHECK_EQUAL(5, pedigree.num_nodes());

    std::vector<peel::family_members_t> family = pedigree.family_members_;
    std::vector<peel::family_members_t> expected_family = {
            {1, 4},
            {0, 1, 2},
            {0, 3}
    };

    for (int f = 0; f < expected_family.size(); ++f) {
        boost_check_equal_vector(expected_family, family);
    }

	std::vector<decltype(peel::op::NUM)> ops = pedigree.peeling_ops_;
	std::vector<decltype(peel::op::NUM)> expected_ops = { peel::op::UP,
			peel::op::TOFATHER, peel::op::UP };
	boost_check_equal_vector(expected_ops, ops);

	std::vector<decltype(peel::op::NUM)> functions_ops =
			pedigree.peeling_functions_ops_;
	std::vector<decltype(peel::op::NUM)> expected_functions_ops = {
			peel::op::UPFAST, peel::op::TOFATHERFAST, peel::op::UP };

    boost_check_equal_vector(expected_functions_ops, functions_ops);

}
typedef std::pair<int, int> PairIndex;
void boost_check_equal_pair_index(PairIndex expected_index, PairIndex result_index){

    BOOST_CHECK_EQUAL(expected_index.first, result_index.first);
    BOOST_CHECK_EQUAL(expected_index.second, result_index.second);

}

BOOST_FIXTURE_TEST_CASE_NO_DECOR(test_parse_io_pedigree, ReadTrioFromFile) {

    RelationshipGraph relationship_graph;

    relationship_graph.SetupFirstNodeIndex(ped);
    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(4, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(3, relationship_graph.first_nonfounder_);

    // Construct a graph of the pedigree and somatic information
    Graph pedigree_graph(relationship_graph.first_somatic_);
//
//        PrintDebugEdges("========== VERISON 2 =========\ninit pedigree", pedigree_graph);
//
//        auto edge_types = get(boost::edge_type, pedigree_graph);
//        auto lengths = get(boost::edge_length, pedigree_graph);
//        auto labels = get(boost::vertex_label, pedigree_graph);
//    //    auto groups = get(vertex_group, pedigree_graph);
//    //    auto families = get(edge_family, pedigree_graph);
//
//        // Add the labels for the germline nodes
//        labels[0] = DNG_GL_PREFIX "unknown";
//        for (size_t i = 1; i < first_somatic_; ++i) {
//            labels[i] = DNG_GL_PREFIX + pedigree.name(i);
//        }
//
//
//        // Go through rows and construct the pedigree part.
//
    relationship_graph.ParseIoPedigree(pedigree_graph, ped);
    BOOST_CHECK_EQUAL(7, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(6, num_edges(pedigree_graph));

    RelationshipGraph::PropVertexIndex index = get(boost::vertex_index, pedigree_graph);


//    typedef std::pair<RelationshipGraph::PropVertexIndex,
//            RelationshipGraph::PropVertexIndex> PairIndex;
//    typedef std::pair<std::size_t, std::size_t> PairIndex;
    typedef std::pair<int, int> PairIndex;
    std::vector<PairIndex> index_vector;

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {

        PairIndex pi = std::make_pair(index[source(*ei, pedigree_graph)],
                index[target(*ei, pedigree_graph)]);
        index_vector.push_back(pi);
    }

    sort(index_vector.begin(), index_vector.end());

    std::vector<PairIndex> expected_vector {
            std::make_pair(1,2),
            std::make_pair(1,3),
            std::make_pair(1,4),
            std::make_pair(2,3),
            std::make_pair(2,5),
            std::make_pair(3,6) };
    for (int i = 0; i < expected_vector.size(); ++i) {
        auto result = std::find(index_vector.begin(), index_vector.end(), expected_vector[i]);
//        BOOST_CHECK_EQUAL(i, result-index_vector.begin());
//        std::cout << (expected_vector[i] == index_vector[i]) << std::endl;
//        BOOST_CHECK(expected_vector[i] == index_vector[i]);
        boost_check_equal_pair_index(expected_vector[i], index_vector[i]);
    }

//    a = std::search(index_vector.begin(), index_vector.end(), std::make_pair(2,3));
//    std::cout << (a==index_vector.end()) << std::endl;

//        std::cout << " ==END==" << std::endl;
//        std::cout << "Founder, Non_F, Lib, Somatic: " << first_founder_ << "\t"
//                << first_nonfounder_ << "\t" << first_library_ << "\t"
//                << first_somatic_ << std::endl;
//        std::cout << std::endl;
//}

}


} // namespace dng

