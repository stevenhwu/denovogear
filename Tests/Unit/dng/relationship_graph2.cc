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


#define BOOST_TEST_MODULE dng::relationship_graph2

#include <dng/relationship_graph.h>

#include <iostream>

#include "../boost_test_helper.h"
#include "fixture_read_trio_from_file.h"

namespace utf = boost::unit_test;

#define DNG_GL_PREFIX "GL-" //HACK: Refactor this later

struct FixturePedigreeMid {


    dng::RelationshipGraph pedigree;
//    dng::PedigreeV2 pedigree_v2;

    std::string fixture;

    dng::io::Pedigree io_pedigree;
    dng::ReadGroups rgs;

    typedef dng::task::Call task_type;
    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    FixturePedigreeMid(std::string s = "FixturePedigreeMid") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        po::options_description ext_desc, int_desc;
        po::positional_options_description pos_desc;
        po::variables_map vm;

        int argc=4;
        char *argv[argc];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";

        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/relationship_graph.ped");
        argv[2] = (char*) ped_filename.data();

        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/relationship_graph.vcf");
        argv[3] = (char*) vcf_filename.data();

        add_app_args(ext_desc,
                     static_cast<typename task_type::argument_type &>(arg));
        int_desc.add_options()("input",
                               po::value<std::vector<std::string> >(&arg.input),
                               "input files");
        int_desc.add(ext_desc);
        pos_desc.add("input", -1);
        po::store(
                po::command_line_parser(argc, argv).options(int_desc).positional(
                        pos_desc).run(), vm);
        po::notify(vm);

        // Parse pedigree from file

        std::ifstream ped_file(arg.ped);
        io_pedigree.Parse(istreambuf_range(ped_file));


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
        arg.mu= 0.05;
        arg.mu_somatic = 0.07;
        arg.mu_library = 0.11;

    }

    ~FixturePedigreeMid() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }

};


[[deprecated]] typedef std::pair<int, int> PairIndex;
//TODO(SW): tuple vs struct?
typedef std::tuple<int, int, graph::EdgeType, float> EdgeInfo;

struct EdgeInfo2{
    int source_vertex;
    int target_vertex;
    graph::EdgeType type;
    float edge_length;
};

std::vector<EdgeInfo> extract_edge_info(Graph &pedigree_graph) {

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto edge_length = get(boost::edge_length, pedigree_graph);
    auto node_index = get(boost::vertex_index, pedigree_graph);
    std::vector<EdgeInfo> edge_info_vector;
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        EdgeInfo pi = std::make_tuple(node_index[source(*ei, pedigree_graph)],
                                      node_index[target(*ei, pedigree_graph)],
                                      edge_types[*ei], edge_length[*ei]);
//        std::cout <<  static_cast<std::underlying_type<graph::EdgeType>::type>(edge_types[*ei]) << std::cout;
//        std::cout <<  static_cast<std::underlying_type<graph::EdgeType>::type> (graph::EdgeType::Library) << std::endl;
//        std::cout <<  static_cast<int> (graph::EdgeType::Library) << std::endl;
//        std::cout <<  static_cast<int> (edge_types[*ei]) << "\t" << edge_length[*ei] << std::endl;
        edge_info_vector.push_back(pi);
    }
    sort(edge_info_vector.begin(), edge_info_vector.end());
    return edge_info_vector;
}

//void boost_check_equal_pair_index(PairIndex expected_index, PairIndex result_index){
//    BOOST_CHECK_EQUAL(expected_index.first, result_index.first);
//    BOOST_CHECK_EQUAL(expected_index.second, result_index.second);
//}


void boost_check_equal_edge(EdgeInfo expected, EdgeInfo actual){

    BOOST_CHECK_EQUAL(std::get<0>(expected), std::get<0>(actual));
    BOOST_CHECK_EQUAL(std::get<1>(expected), std::get<1>(actual));
    BOOST_CHECK(std::get<2>(expected) == std::get<2>(actual));
    BOOST_CHECK_EQUAL(std::get<3>(expected), std::get<3>(actual));
    BOOST_TEST_MESSAGE("Edge:" << std::get<0>(expected) << "-" << std::get<1>(expected) <<
            " Actual: " << std::get<0>(actual) << "-" << std::get<1>(actual) ) ;
}
/*

0     1
|     |
---|---
|  |  |
3  2  4

*/
namespace dng {
BOOST_FIXTURE_TEST_CASE(test_constructor, FixturePedigreeMid ) {

    pedigree.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);

    BOOST_CHECK_EQUAL(22, pedigree.num_nodes() );

    auto workspace = pedigree.CreateWorkspace();
    BOOST_CHECK_EQUAL(0, workspace.founder_nodes.first);
    BOOST_CHECK_EQUAL(6, workspace.founder_nodes.second);
    BOOST_CHECK_EQUAL(0, workspace.germline_nodes.first);
    BOOST_CHECK_EQUAL(10, workspace.germline_nodes.second);
    BOOST_CHECK_EQUAL(10, workspace.somatic_nodes.first);
    BOOST_CHECK_EQUAL(10, workspace.somatic_nodes.second);
    BOOST_CHECK_EQUAL(10, workspace.library_nodes.first);
    BOOST_CHECK_EQUAL(22, workspace.library_nodes.second);

    auto labels = pedigree.labels();

    const std::vector<std::string> expected_labels = {
            "GL-1", "GL-2",
            "GL-4", "GL-5",
            "GL-9", "GL-10",
            "GL-3", "GL-6",
            "GL-11", "GL-8",
            "LB-NA12001:Solexa-001", "LB-NA12002:Solexa-002",
            "LB-NA12003:Solexa-003", "LB-NA12004:Solexa-004",
            "LB-NA12005:Solexa-005", "LB-NA12006:Solexa-006",
            "LB-NA12007:Solexa-007", "LB-NA12008:Solexa-008",
            "LB-NA12009:Solexa-009", "LB-NA12010:Solexa-010",
            "LB-NA12011:Solexa-011", "LB-NA12012:Solexa-012"
    };
    boost_check_equal_vector(expected_labels, labels);


    auto transitions = pedigree.transitions();
    auto s_max = static_cast<size_t>(-1);
    double mu_somatic_library = arg.mu + arg.mu_somatic + arg.mu_library;
    double somatic_library = arg.mu_somatic + arg.mu_library;
    std::vector<RelationshipGraph::transition_t> expected_transitions = {
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
            {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},

            {RelationshipGraph::TransitionType::Germline, 0, 1, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 2, 3, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 5, 4, arg.mu, arg.mu},
            {RelationshipGraph::TransitionType::Germline, 6, 7, arg.mu, arg.mu},

            {RelationshipGraph::TransitionType::Somatic, 0, s_max,
                    somatic_library, 0},
			{RelationshipGraph::TransitionType::Somatic, 1, s_max,
			        somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 6, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 2, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 3, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 7, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Germline, 6, 7,
                    mu_somatic_library, mu_somatic_library},
            {RelationshipGraph::TransitionType::Somatic, 9, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Somatic, 4, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 5, s_max,
                    somatic_library, 0},
            {RelationshipGraph::TransitionType::Somatic, 8, s_max,
                    somatic_library, 0},

            {RelationshipGraph::TransitionType::Germline, 8, 9,
                    mu_somatic_library, mu_somatic_library}
    };

    BOOST_CHECK_EQUAL(expected_transitions.size(), transitions.size());
    for (int k = 0; k < expected_transitions.size(); ++k) {
        std::cout << k << std::endl;
        auto expected = expected_transitions[k];
        auto actual = transitions[k];
        BOOST_CHECK(expected.type == actual.type);
        BOOST_CHECK_EQUAL(expected.parent1, actual.parent1);
        BOOST_CHECK_EQUAL(expected.parent2, actual.parent2);
        BOOST_CHECK_CLOSE(expected.length1, actual.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(expected.length2, actual.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);

    }


}

//1-2    3-4
// |      |
// 7------8   5-6
//   |  |      |
//   10 11-----9
//          |
//          12



BOOST_FIXTURE_TEST_CASE(test_pedigree_inspect, FixturePedigreeMid) {

    pedigree.Construct(io_pedigree, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
    BOOST_CHECK_EQUAL(22, pedigree.num_nodes());

    std::vector<peel::family_members_t> family = pedigree.family_members_;
    std::vector<peel::family_members_t> expected_family = {
            {2, 13},
            {3, 14},
            {2, 3, 7},
            {5, 19},
            {4, 18},
            {5, 4, 8},
            {8, 20},
            {8, 9, 21},
            {9, 17},
            {7, 15},
            {6, 7, 9, 16},
            {6, 12},
            {1, 11},
            {0, 1, 6},
            {0, 10}
    };

    for (int f = 0; f < expected_family.size(); ++f) {
        boost_check_equal_vector(expected_family[f], family[f]);
    }

	std::vector<decltype(peel::op::NUM)> ops = pedigree.peeling_ops_;
    std::vector<decltype(peel::op::NUM)> expected_ops = {
            peel::op::UP,
            peel::op::UP,
			peel::op::TOCHILD,
			peel::op::UP,
			peel::op::UP,
			peel::op::TOCHILD,
			peel::op::UP,
            peel::op::TOMOTHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP
    };
	boost_check_equal_vector(expected_ops, ops);

	std::vector<decltype(peel::op::NUM)> functions_ops =
			pedigree.peeling_functions_ops_;
    std::vector<decltype(peel::op::NUM)> expected_functions_ops = {
            peel::op::UPFAST,
            peel::op::UPFAST,
            peel::op::TOCHILDFAST,
            peel::op::UPFAST,
            peel::op::UPFAST,
            peel::op::TOCHILDFAST,
            peel::op::UPFAST,
            peel::op::TOMOTHERFAST,
            peel::op::UP,
            peel::op::UPFAST,
            peel::op::TOFATHERFAST,
            peel::op::UP,
            peel::op::UPFAST,
            peel::op::TOFATHERFAST,
            peel::op::UP
    };
    boost_check_equal_vector(expected_functions_ops, functions_ops);

}


BOOST_FIXTURE_TEST_CASE(test_parse_io_pedigree, FixturePedigreeMid){

    RelationshipGraph relationship_graph;

    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(7, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(13, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(0, relationship_graph.first_library_);

    Graph pedigree_graph(relationship_graph.first_somatic_);

    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    BOOST_CHECK_EQUAL(25, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(29, num_edges(pedigree_graph));

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    std::vector<EdgeInfo> expected_vector {
            std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
            std::make_tuple(1, 7, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(1, 13, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(2, 7, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(2, 14, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(3, 4, graph::EdgeType::Spousal, 0),
            std::make_tuple(3, 8, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(3, 15, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(4, 8, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(4, 16, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(5, 9, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(5, 17, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(6, 5, graph::EdgeType::Spousal, 0),
            std::make_tuple(6, 9, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(6, 18, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(7, 8, graph::EdgeType::Spousal, 0),
            std::make_tuple(7, 10, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(7, 11, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(7, 19, graph::EdgeType::Mitotic, 1.0),
            std::make_tuple(8, 10, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(8, 11, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(8, 20, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(9, 11, graph::EdgeType::Spousal, 0),
            std::make_tuple(9, 12, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(9, 21, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(10, 22, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(11, 12, graph::EdgeType::Meiotic, 1.0),
            std::make_tuple(11, 23, graph::EdgeType::Mitotic, 1.0),

            std::make_tuple(12, 24, graph::EdgeType::Mitotic, 1.0),

    };
    BOOST_CHECK_EQUAL(expected_vector.size(), edge_vector.size());
    for (int i = 0; i < expected_vector.size(); ++i) {
        boost_check_equal_edge(expected_vector[i], edge_vector[i]);
    }


}


BOOST_FIXTURE_TEST_CASE(test_add_lib_from_rgs, FixturePedigreeMid) {


    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(7, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(13, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(25, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(37, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(37, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(41, num_edges(pedigree_graph));

    std::vector<EdgeInfo> expected_edges{
        std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
        std::make_tuple(1, 7, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(1, 13, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(2, 7, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(2, 14, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(3, 4, graph::EdgeType::Spousal, 0),
        std::make_tuple(3, 8, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(3, 15, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(4, 8, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(4, 16, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(5, 9, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(5, 17, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(6, 5, graph::EdgeType::Spousal, 0),
        std::make_tuple(6, 9, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(6, 18, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(7, 8, graph::EdgeType::Spousal, 0),
        std::make_tuple(7, 10, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(7, 11, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(7, 19, graph::EdgeType::Mitotic, 1.0),
        std::make_tuple(8, 10, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(8, 11, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(8, 20, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(9, 11, graph::EdgeType::Spousal, 0),
        std::make_tuple(9, 12, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(9, 21, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(10, 22, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(11, 12, graph::EdgeType::Meiotic, 1.0),
        std::make_tuple(11, 23, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(12, 24, graph::EdgeType::Mitotic, 1.0),

        std::make_tuple(13, 25, graph::EdgeType::Library, 1.0),
        std::make_tuple(14, 26, graph::EdgeType::Library, 1.0),
        std::make_tuple(15, 28, graph::EdgeType::Library, 1.0),
        std::make_tuple(16, 29, graph::EdgeType::Library, 1.0),
        std::make_tuple(17, 33, graph::EdgeType::Library, 1.0),
        std::make_tuple(18, 34, graph::EdgeType::Library, 1.0),
        std::make_tuple(19, 27, graph::EdgeType::Library, 1.0),
        std::make_tuple(20, 30, graph::EdgeType::Library, 1.0),
        std::make_tuple(21, 35, graph::EdgeType::Library, 1.0),
        std::make_tuple(22, 31, graph::EdgeType::Library, 1.0),
        std::make_tuple(23, 32, graph::EdgeType::Library, 1.0),
        std::make_tuple(24, 36, graph::EdgeType::Library, 1.0),

    };

    std::vector<std::string> expected_vertex{
        "GL-unknown",
        "GL-1", "GL-2",
        "GL-4", "GL-5",
        "GL-9", "GL-10",
        "GL-3", "GL-6", "GL-11",
        "GL-7", "GL-8",
        "GL-12",

        "SM-NA12001", "SM-NA12002",
        "SM-NA12004", "SM-NA12005",
        "SM-NA12009", "SM-NA12010",
        "SM-NA12003", "SM-NA12006", "SM-NA12011",
        "SM-NA12007", "SM-NA12008",
        "SM-NA12012",

        "LB-NA12001:Solexa-001", "LB-NA12002:Solexa-002",
        "LB-NA12003:Solexa-003", "LB-NA12004:Solexa-004",
        "LB-NA12005:Solexa-005", "LB-NA12006:Solexa-006",
        "LB-NA12007:Solexa-007", "LB-NA12008:Solexa-008",
        "LB-NA12009:Solexa-009", "LB-NA12010:Solexa-010",
        "LB-NA12011:Solexa-011", "LB-NA12012:Solexa-012"

    };

    for (int v = 0; v < relationship_graph.num_nodes_; ++v) {
        BOOST_CHECK_EQUAL(expected_vertex[v], labels[v]);
    }

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }

}


BOOST_FIXTURE_TEST_CASE(test_update_edge_lengths, FixturePedigreeMid) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    std::vector<EdgeInfo> expected_edges{
        std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
        std::make_tuple(1, 7, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(1, 13, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(2, 7, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(2, 14, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(3, 4, graph::EdgeType::Spousal, 0),
        std::make_tuple(3, 8, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(3, 15, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(4, 8, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(4, 16, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(5, 9, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(5, 17, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(6, 5, graph::EdgeType::Spousal, 0),
        std::make_tuple(6, 9, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(6, 18, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(7, 8, graph::EdgeType::Spousal, 0),
        std::make_tuple(7, 10, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(7, 11, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(7, 19, graph::EdgeType::Mitotic, 0.07),
        std::make_tuple(8, 10, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(8, 11, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(8, 20, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(9, 11, graph::EdgeType::Spousal, 0),
        std::make_tuple(9, 12, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(9, 21, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(10, 22, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(11, 12, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(11, 23, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(12, 24, graph::EdgeType::Mitotic, 0.07),

        std::make_tuple(13, 25, graph::EdgeType::Library, 0.11),
        std::make_tuple(14, 26, graph::EdgeType::Library, 0.11),
        std::make_tuple(15, 28, graph::EdgeType::Library, 0.11),
        std::make_tuple(16, 29, graph::EdgeType::Library, 0.11),
        std::make_tuple(17, 33, graph::EdgeType::Library, 0.11),
        std::make_tuple(18, 34, graph::EdgeType::Library, 0.11),
        std::make_tuple(19, 27, graph::EdgeType::Library, 0.11),
        std::make_tuple(20, 30, graph::EdgeType::Library, 0.11),
        std::make_tuple(21, 35, graph::EdgeType::Library, 0.11),
        std::make_tuple(22, 31, graph::EdgeType::Library, 0.11),
        std::make_tuple(23, 32, graph::EdgeType::Library, 0.11),
        std::make_tuple(24, 36, graph::EdgeType::Library, 0.11),

    };

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }

}


BOOST_FIXTURE_TEST_CASE(test_simplify_pedigree, FixturePedigreeMid) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(7, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(13, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(25, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(37, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(37, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(27, num_edges(pedigree_graph));

    std::vector<EdgeInfo> expected_edges{

        std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
        std::make_tuple(1, 7, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(1, 25, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(2, 7, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(2, 26, graph::EdgeType::Mitotic, 0.18),

        std::make_tuple(3, 4, graph::EdgeType::Spousal, 0),
        std::make_tuple(3, 8, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(3, 28, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(4, 8, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(4, 29, graph::EdgeType::Mitotic, 0.18),

        std::make_tuple(5, 9, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(5, 33, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(6, 5, graph::EdgeType::Spousal, 0),
        std::make_tuple(6, 9, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(6, 34, graph::EdgeType::Mitotic, 0.18),

        std::make_tuple(7, 8, graph::EdgeType::Spousal, 0),
        std::make_tuple(7, 11, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(7, 27, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(7, 31, graph::EdgeType::Meiotic, 0.23),

        std::make_tuple(8, 11, graph::EdgeType::Meiotic, 0.05),
        std::make_tuple(8, 30, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(8, 31, graph::EdgeType::Meiotic, 0.23),


        std::make_tuple(9, 11, graph::EdgeType::Spousal, 0),
        std::make_tuple(9, 35, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(9, 36, graph::EdgeType::Meiotic, 0.23),

        std::make_tuple(11, 32, graph::EdgeType::Mitotic, 0.18),
        std::make_tuple(11, 36, graph::EdgeType::Meiotic, 0.23),

    };

    BOOST_CHECK_EQUAL(expected_edges.size(), edge_vector.size());
    for (int i = 0; i < expected_edges.size(); ++i) {
        boost_check_equal_edge(expected_edges[i], edge_vector[i]);
    }

}


BOOST_FIXTURE_TEST_CASE(test_update_labels_node_ids, FixturePedigreeMid) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);
    std::vector<size_t> node_ids(relationship_graph.num_nodes_, -1);
    relationship_graph.UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    std::vector<EdgeInfo> edge_vector = extract_edge_info(pedigree_graph);
    auto labels = get(boost::vertex_label, pedigree_graph);

    BOOST_CHECK_EQUAL(0, relationship_graph.first_founder_);
    BOOST_CHECK_EQUAL(6, relationship_graph.first_nonfounder_);
    BOOST_CHECK_EQUAL(10, relationship_graph.first_somatic_);
    BOOST_CHECK_EQUAL(10, relationship_graph.first_library_);
    BOOST_CHECK_EQUAL(22, relationship_graph.num_nodes_);
    BOOST_CHECK_EQUAL(37, num_vertices(pedigree_graph));
    BOOST_CHECK_EQUAL(27, num_edges(pedigree_graph));

    //TODO(SW): which one, same "result" but different ideas
    auto s_max = static_cast<std::size_t>(-1);
//    auto size_t_max = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> expected_node_ids = {
            s_max, 0, 1, 2, 3,
            4, 5, 6, 7, 8,
            s_max, 9, s_max, s_max, s_max,
            s_max, s_max, s_max, s_max, s_max,
            s_max, s_max, s_max, s_max, s_max,

            10, 11, 12, 13, 14, 15,
            16, 17, 18, 19, 20, 21
            };

    std::vector<std::string> expected_labels {

        "GL-1", "GL-2",
        "GL-4", "GL-5",
        "GL-9", "GL-10",
        "GL-3", "GL-6", "GL-11",
//        "GL-7",
        "GL-8",
//        "GL-12",

        "LB-NA12001:Solexa-001", "LB-NA12002:Solexa-002",
        "LB-NA12003:Solexa-003", "LB-NA12004:Solexa-004",
        "LB-NA12005:Solexa-005", "LB-NA12006:Solexa-006",
        "LB-NA12007:Solexa-007", "LB-NA12008:Solexa-008",
        "LB-NA12009:Solexa-009", "LB-NA12010:Solexa-010",
        "LB-NA12011:Solexa-011", "LB-NA12012:Solexa-012"


    };

    boost_check_equal_vector(expected_labels, relationship_graph.labels_);
    boost_check_equal_vector(expected_node_ids, node_ids);
}


BOOST_FIXTURE_TEST_CASE(test_create_families_info, FixturePedigreeMid) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);
    std::vector<size_t> node_ids(relationship_graph.num_nodes_, -1);
    relationship_graph.UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    RelationshipGraph::family_labels_t family_labels;//(num_families);
    std::vector<vertex_t> pivots;//(num_families, dummy_index);
    relationship_graph.CreateFamiliesInfo(pedigree_graph, family_labels, pivots);


    std::vector<std::vector<EdgeInfo>> expected_family_labels {

        {std::make_tuple(3, 28, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(4, 29, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(3, 4, graph::EdgeType::Spousal, 0),
                std::make_tuple(4, 8, graph::EdgeType::Meiotic, 0.05),
                std::make_tuple(3, 8, graph::EdgeType::Meiotic, 0.05)
        },
        {std::make_tuple(6, 34, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(5, 33, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(6, 5, graph::EdgeType::Spousal, 0),
                std::make_tuple(5, 9, graph::EdgeType::Meiotic, 0.05),
                std::make_tuple(6, 9, graph::EdgeType::Meiotic, 0.05)
        },
        {std::make_tuple(9, 35, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(9, 11, graph::EdgeType::Spousal, 0),
                std::make_tuple(11, 36, graph::EdgeType::Meiotic, 0.23),
                std::make_tuple(9, 36, graph::EdgeType::Meiotic, 0.23)
        },
        {std::make_tuple(11, 32, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(8, 30, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(7, 8, graph::EdgeType::Spousal, 0),
                std::make_tuple(8, 11, graph::EdgeType::Meiotic, 0.05),
                std::make_tuple(7, 11, graph::EdgeType::Meiotic, 0.05),
                std::make_tuple(8, 31, graph::EdgeType::Meiotic, 0.23),
                std::make_tuple(7, 31, graph::EdgeType::Meiotic, 0.23)
        },
        {std::make_tuple(7, 27, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(2, 26, graph::EdgeType::Mitotic, 0.18)},
        {std::make_tuple(1, 2, graph::EdgeType::Spousal, 0),
                std::make_tuple(2, 7, graph::EdgeType::Meiotic, 0.05),
                std::make_tuple(1, 7, graph::EdgeType::Meiotic, 0.05)
        },
        {std::make_tuple(1, 25, graph::EdgeType::Mitotic, 0.18)}

    };
    std::vector<vertex_t> expected_pivots {3, 4, 8, 6, 5, 9, 9, 11, 11, 8, 7, 7, 2, 1, 0};

    boost_check_equal_vector(expected_pivots, pivots);
//    auto families = get(boost::edge_family, pedigree_graph);

    BOOST_CHECK_EQUAL(expected_family_labels.size(), family_labels.size());
    for (int i = 0; i < expected_family_labels.size(); ++i) {
        std::cout << "Families : " << i << "\tPivots: " << pivots[i] << "\t";
        boost::detail::edge_desc_impl<boost::undirected_tag, unsigned long> x;
        for (int j = 0; j < expected_family_labels[i].size(); ++j) {
            auto f = family_labels[i][j];

//                auto edge_types = get(boost::edge_type, pedigree_graph);
//                auto edge_length = get(boost::edge_length, pedigree_graph);
//                auto node_index = get(boost::vertex_index, pedigree_graph);
//                EdgeInfo pi = std::make_tuple(node_index[source(f, pedigree_graph)],
//                            node_index[target(j, pedigree_graph)],
//                            edge_types[j], edge_length[j]);
            EdgeInfo actual = std::make_tuple(f.m_source, f.m_target,
                    boost::get(boost::edge_type, pedigree_graph, f),
                    boost::get(boost::edge_length, pedigree_graph, f) );

            boost_check_equal_edge(expected_family_labels[i][j], actual);

        }
    }

}


BOOST_FIXTURE_TEST_CASE(test_create_peeling_ops, FixturePedigreeMid) {

    Graph pedigree_graph(io_pedigree.member_count());

    RelationshipGraph relationship_graph;
    relationship_graph.SetupFirstNodeIndex(io_pedigree);
    relationship_graph.ParseIoPedigree(pedigree_graph, io_pedigree);
    relationship_graph.AddLibrariesFromReadGroups(pedigree_graph, rgs);
    double expected_mu = 0.05;
    double expected_mu_somatic = 0.07;
    double expected_mu_library = 0.11;
    relationship_graph.UpdateEdgeLengths(pedigree_graph, expected_mu,
                                         expected_mu_somatic,
                                         expected_mu_library);
    relationship_graph.SimplifyPedigree(pedigree_graph);
    std::vector<size_t> node_ids(relationship_graph.num_nodes_, -1);
    relationship_graph.UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    RelationshipGraph::family_labels_t family_labels;//(num_families);
    std::vector<vertex_t> pivots;//(num_families, dummy_index);
    relationship_graph.CreateFamiliesInfo(pedigree_graph, family_labels, pivots);
    relationship_graph.CreatePeelingOps(pedigree_graph, node_ids, family_labels, pivots);


    std::vector<std::size_t> expected_roots {0};
    size_t s_max = static_cast<size_t>(-1);
    std::vector<RelationshipGraph::transition_t> expected_transitions {
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},
        {RelationshipGraph::TransitionType::Founder, s_max, s_max, 0, 0},

        {RelationshipGraph::TransitionType::Germline, 0, 1, 0.05, 0.05},
        {RelationshipGraph::TransitionType::Germline, 2, 3, 0.05, 0.05},
        {RelationshipGraph::TransitionType::Germline, 5, 4, 0.05, 0.05},
        {RelationshipGraph::TransitionType::Germline, 6, 7, 0.05, 0.05},

        {RelationshipGraph::TransitionType::Somatic, 0, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 1, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 6, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 2, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 3, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 7, s_max, 0.18, 0},

        {RelationshipGraph::TransitionType::Germline, 6, 7, 0.23, 0.23},
        {RelationshipGraph::TransitionType::Somatic, 9, s_max, 0.18, 0},

        {RelationshipGraph::TransitionType::Somatic, 4, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 5, s_max, 0.18, 0},
        {RelationshipGraph::TransitionType::Somatic, 8, s_max, 0.18, 0},

        {RelationshipGraph::TransitionType::Germline, 8, 9, 0.23, 0.23}

    };

    std::vector<decltype(peel::op::NUM)> expected_peeling_ops = {
            peel::op::UP,
            peel::op::UP,
            peel::op::TOCHILD,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOCHILD,
            peel::op::UP,
            peel::op::TOMOTHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP,
            peel::op::UP,
            peel::op::TOFATHER,
            peel::op::UP
    };


    std::vector<peel::family_members_t> expected_family_members {
        {2,13},
        {3,14},
        {2,3,7},

        {5,19},
        {4,18},
        {5,4,8},

        {8,20},
        {8,9,21},

        {9,17},
        {7,15},
        {6,7,9,16},

        {6,12},
        {1,11},
        {0,1,6},
        {0,10}
    };

    boost_check_equal_vector(expected_roots, relationship_graph.roots_);

    BOOST_CHECK_EQUAL(expected_transitions.size(), relationship_graph.transitions_.size());
    for (int i = 0; i < expected_transitions.size(); ++i) {
        auto exp = expected_transitions[i];
        auto acutal = relationship_graph.transitions_[i];

        BOOST_CHECK(exp.type == acutal.type);
//        BOOST_CHECK_EQUAL(static_cast<int>(exp.type), static_cast<int>(acutal.type));
        BOOST_CHECK_EQUAL(exp.parent1, acutal.parent1);
        BOOST_CHECK_EQUAL(exp.parent2, acutal.parent2);
        BOOST_CHECK_CLOSE(exp.length1, acutal.length1, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
        BOOST_CHECK_CLOSE(exp.length2, acutal.length2, BOOST_CLOSE_PERCENTAGE_THRESHOLD);
    }

    boost_check_equal_vector(expected_peeling_ops, relationship_graph.peeling_ops_);

    BOOST_CHECK_EQUAL(expected_family_members.size(), relationship_graph.family_members_.size());
    for (int i = 0; i < expected_family_members.size(); ++i) {
        auto exp = expected_family_members[i];
        auto actual = relationship_graph.family_members_[i];
        boost_check_equal_vector(exp, actual);
    }



}


} // namespace dng





