/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#include <dng/pedigree_v2.h>
#include <dng/graph.h>
#include <dng/mutation.h>

#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/timer/timer.hpp>
#include "../vt/assert_utils.h"

#define DNG_GL_PREFIX "GL-"
#define DNG_SM_PREFIX "SM-" // define also in newick.cc
#define DNG_LB_PREFIX "LB-"

/*
RULES FOR LINKING READ GROUPS TO PEOPLE.

0) Read all read groups in bam files. Map RG to Library to Sample.

1) If tissue information is present, build a tissue tree connecting samples
   to zygotic genotypes.  Build mutation matrices for every unique branch
   length.

2) If no tissue information is present in pedigree, check to see if there is a
   sample with the same label as the individual.  If so, connect that sample to
   the individual with a somatic branch of length 1.

3) If a sample has multiple libraries, connect the libraries to sample with a
   library branch.

4) If a library has multiple read-groups, concat the read-groups.

*/

bool dng::PedigreeV2::Construct(const io::Pedigree &pedigree,
                              dng::ReadGroups &rgs,
                              double mu, double mu_somatic, double mu_library) {
    using namespace boost;
    using namespace std;

    //TODO: try to remove these if possible
//    graph_traits<Graph>::edge_iterator ei, ei_end;
//    graph_traits<Graph>::vertex_iterator vi, vi_end;

//Maybe not here
//    typedef property_map<Graph, vertex_index_t>::type IndexMap;

    // Setup counts //const??
    std::size_t num_members = pedigree.member_count();
//    const vertex_t dummy_index = 0;

    // Count the founders
    first_founder_ = 0; //Can this be something other than 0??
    first_somatic_ = num_members;

    for (first_nonfounder_ = first_founder_; first_nonfounder_ < num_members;
         ++first_nonfounder_) {
        if (pedigree.table()[first_nonfounder_].dad != 0 &&
            pedigree.table()[first_nonfounder_].mom != 0) {
            break;
        }
    }
//    std::size_t num_libraries = rgs.libraries().size(); //??const??

    // Construct a graph of the pedigree and somatic information
    Graph pedigree_graph(num_members);

    std::cout << "222222222222222222\ninit pedigree: V E: " << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
//    std::cout << pedigree.table().size() << "\t" << pedigree.table().size() << "\t" << rgs.libraries().size() << std::endl;


    auto edge_types = get(edge_type, pedigree_graph);
    auto lengths = get(edge_length, pedigree_graph);
    auto labels = get(vertex_label, pedigree_graph);
//    auto groups = get(vertex_group, pedigree_graph);
//    auto families = get(edge_family, pedigree_graph);
    // Add the labels for the germline nodes


    labels[0] = DNG_GL_PREFIX "unknown";
    for (size_t i = 1; i < num_members; ++i) {
        labels[i] = DNG_GL_PREFIX + pedigree.name(i);
    }
    PrintEdges("Start", pedigree_graph);
    // Go through rows and construct the pedigree part.


    ParseIoPedigree(pedigree_graph, pedigree);
    AddLibrariesFromReadGroups(pedigree_graph, rgs, labels);


    vector<size_t> node_ids(num_nodes_, -1);

    UpdateEdgeLengths(pedigree_graph,  mu,  mu_somatic,  mu_library);
    SimplifyPedigree(pedigree_graph, edge_types, lengths);
    UpdateLabelsNodeIds(pedigree_graph, rgs, node_ids);

    PrintEdges("After update position()", pedigree_graph);

    // Reset Family Information
    roots_.clear();
    roots_.reserve(16);
    family_members_.clear();
    family_members_.reserve(128);
    peeling_ops_.clear();
    peeling_ops_.reserve(128);

    // Calculate the connected components.  This defines independent sections
    // of the graph.

//    std::size_t num_groups = connected_components(pedigree_graph, groups);

    // Calculate the biconnected components and articulation points.
    // This defines "nuclear" families and pivot individuals.
    // Nodes which have no edges will not be part of any family.
//    vector<vertex_t> articulation_vertices;
//    std::size_t num_families = biconnected_components(pedigree_graph, families,
//                                                      back_inserter(articulation_vertices)).first;

//    //LOGIC:: only need 2 variables??
    family_labels_t family_labels;//(num_families);
    vector<vertex_t> pivots;//(num_families, dummy_index);
//
    CreateFamiliesInfo(pedigree_graph, family_labels, pivots);


    CreatePeelingOps(pedigree_graph, node_ids, family_labels, pivots);

    ConstructPeelingMachine();

#ifdef DNG_DEVEL
    PrintMachine(cerr);
#else
    //PrintMachine(cerr);
#endif

    return true;
}


void dng::PedigreeV2::UpdateLabelsNodeIds(dng::Graph &pedigree_graph, dng::ReadGroups rgs, std::vector<size_t> &node_ids){

    auto labels = get(boost::vertex_label, pedigree_graph);

    std::cout << num_nodes_ << std::endl;
    labels_.clear();
    labels_.reserve(128);
    size_t vid = 0;
    for(std::size_t u = 0; u < num_nodes_ ; ++u) {
        if(out_degree(u, pedigree_graph) == 0) {
            continue;
        }
        if(!labels[u].empty()) {
            labels_.push_back(labels[u]);
        } else {
            labels_.push_back(DNG_SM_PREFIX "unnamed_node_" + util::to_pretty(vid));
        }
        node_ids[u] = vid++;
    }
    num_nodes_ = vid;
    std::cout << num_nodes_ << std::endl;



    std::cout << "After remove some nodes: new V labels_sizse: " << vid << "\t" << labels_.size() << std::endl;
    for (auto item : labels_) {
        std::cout << item << std::endl;
    }
    std::cout << "node_ids:"<< std::endl;
    for (int j = 0; j < node_ids.size(); ++j) {
        std::cout << j << " -> " << node_ids[j] << std::endl;
    }

    EraseRemovedLibraries(rgs, node_ids);//TODO: Do we really need this??


//    size_t vid = num_nodes_;
    auto update_position = [&node_ids, vid](size_t pos) -> size_t {
        for(; pos < node_ids.size() && node_ids[pos] == -1; ++pos)
            /*noop*/;
        return (pos < node_ids.size()) ? node_ids[pos] : vid;
    };

    first_founder_ = update_position(first_founder_);
    first_nonfounder_ = update_position(first_nonfounder_);
    first_somatic_ = update_position(first_somatic_);
    first_library_ = update_position(first_library_);

    PrintEdges("Test", pedigree_graph);



}


void dng::PedigreeV2::UpdateEdgeLengths(dng::Graph &pedigree_graph, double mu, double mu_somatic, double mu_library) {
    boost::graph_traits<dng::Graph>::edge_iterator ei, ei_end;
    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);

    for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        std::cout << *ei << "\t" << lengths[*ei] << " " ;
        if(edge_types[*ei] == dng::graph::EdgeType::Meiotic) {
            lengths[*ei] *= mu;
        } else if(edge_types[*ei] == dng::graph::EdgeType::Mitotic) {
            lengths[*ei] *= mu_somatic;
        } else if(edge_types[*ei] == dng::graph::EdgeType::Library) {
            lengths[*ei] *= mu_library;
        }
        std::cout << "->newRate: " << lengths[*ei] << std::endl;
    }
    std::cout << "Founder, Non_F, Lib, Somatic: " << this->first_founder_ << "\t" << this->first_nonfounder_ << "\t" <<
    this->first_library_ << "\t" << this->first_somatic_ << std::endl;
}

void dng::PedigreeV2::AddLibrariesFromReadGroups(Graph &pedigree_graph, const dng::ReadGroups &rgs,
                                            PropVertexLabel &labels) {

    dng::PedigreeV2::PrintEdges("After clear_vertex()", pedigree_graph);

    first_library_ = num_vertices(pedigree_graph);

    // Add library nodes to graph
    for(auto && a : rgs.libraries()) {
        dng::vertex_t v = add_vertex(pedigree_graph);
        labels[v] = DNG_LB_PREFIX + a;
        std::cout << "Add lib: " << labels[v] << std::endl;
    }
    dng::PedigreeV2::PrintEdges("After add libs", pedigree_graph);

    dng::PedigreeV2::ConnectSomaticToLibraries(pedigree_graph, rgs, labels);
    num_nodes_ = num_vertices(pedigree_graph);

    dng::PedigreeV2::PrintEdges("After connect somatic", pedigree_graph);
}
//
//void dng::Pedigree::ConstructPeelingMachine() {
//    using namespace dng::peel;
//    peeling_functions_.clear();
//    peeling_functions_.reserve(peeling_ops_.size());
//    peeling_functions_ops_.reserve(peeling_ops_.size());
//    peeling_reverse_functions_.reserve(peeling_ops_.size());
//    std::vector<std::size_t> lower_written(num_nodes_, -1);
//    for(std::size_t i = 0 ; i < peeling_ops_.size(); ++i) {
//        auto a = peeling_ops_[i];
//        const auto &fam = family_members_[i];
//        auto w = fam[info[a].writes_to];
//        int b = a;
//        switch(a) {
//        case op::DOWN:
//            // If the lower of the parent has never been written to, we can use the fast version
//            b = (lower_written[fam[0]] == -1) ? op::UPFAST + b : b;
//            break;
//        case op::TOCHILD:
//            // If the we only have one child, we can use the fast version
//            b = (fam.size() == 3) ? op::UPFAST + b : b;
//            break;
//        case op::TOMOTHER:
//        case op::TOFATHER:
//        case op::UP:
//            // If the lower of the destination has never been written to, we can use the fast version
//            b = (lower_written[w] == -1) ? op::UPFAST + b : b;
//            break;
//        default:
//            assert(false); // should never get here
//            break;
//        }
//        peeling_functions_ops_.push_back(static_cast<decltype(op::NUM)>(b));
//        peeling_functions_.push_back(functions[b]);
//        peeling_reverse_functions_.push_back(reverse_functions[b]);
//
//        // If the operation writes to a lower value, make note of it
//        if(info[a].writes_lower) {
//            lower_written[w] = i;
//        }
//    }
//}

//

//template <class Config>
//testV(typename Config::vertex_descriptor v){
//
//    std::cout <<"V: "<< v<< std::endl;
//};

//detail {

//=========================================================================
// Adjacency List Generator

//template <class Graph, class VertexListS, class OutEdgeListS,
//        class DirectedS, class VertexProperty, class EdgeProperty,
//        class GraphProperty, class EdgeListS>
//struct adj_list_gen
//using namespace boost;

//typedef boost::detail::edge_descriptor EE;
//template <class Config>
//void testE(EdgeProperties e, dng::Graph &pedigree_graph){
//
//
//     for (dng::vertex_t w = 12; w > 0; --w) {
//         dng::vertex_t v = w - 1;
//
//        auto rng = out_edges(v, pedigree_graph);
//         for (auto it = rng.first; it != rng.second; ++it) {
//             dng::edge_t ee = *it;
//             std::cout << "E: "<< e[ee]<< std::endl;
//         }
//    }
//
//};


void dng::PedigreeV2::ParseIoPedigree(dng::Graph &pedigree_graph, const dng::io::Pedigree &pedigree){

    for (auto &row : pedigree.table()) {

        vertex_t child = row.child;
        vertex_t dad = row.dad;
        vertex_t mom = row.mom;
        std::cout << "===C_DM: " << child  << "\t" << dad << "\t" << mom << std::endl;

        if (child == 0) {
            continue;
        }
        if (dad == mom && dad != 0) {
            // Selfing is not supported

            throw std::runtime_error("Unable to construct peeler for pedigree; selfing is not supported");
//                                             "possible non-zero-loop pedigree.");
        }

        // check to see if mom and dad have been seen before
        auto id = edge(dad, mom, pedigree_graph);

        if (!id.second) {
            add_edge(dad, mom, EdgeType::Spousal, pedigree_graph);
            //Connect dad-mom to make a dependend trio
        }
        std::cout << "===after spoesal: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
        // add the meiotic edges
        add_edge(mom, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);
        add_edge(dad, child, {EdgeType::Meiotic, 1.0f}, pedigree_graph);
        std::cout << "===after add_edge: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
        std::cout << "===C_DM: " << child << "\t" << dad << "\t" << mom << std::endl;
        // Process newick file
        int res = newick::parse(row.sample_tree, child, pedigree_graph);
        std::cout << "===after parse: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
        if (res == 0) {
            // this line has a blank somatic line, so use the name from the pedigree
            vertex_t v = add_vertex(DNG_SM_PREFIX + pedigree.name(child), pedigree_graph);
            add_edge(child, v, {EdgeType::Mitotic, 1.0f}, pedigree_graph);
            std::cout << "===add res==0: " << child << "\t" << v << "\t" << pedigree.name(child) << std::endl;
        } else if (res == -1) {
            throw std::runtime_error(
                    "unable to parse somatic data for individual '" +
                    pedigree.name(child) + "'.");
        }
        std::cout << "Loop END: V E:" << num_vertices(pedigree_graph) << "\t" << num_edges(pedigree_graph) << std::endl;
    }
    PrintEdges("Parsed io::pedigree, before cleas_vertex dummy", pedigree_graph);
//    vertex_t dummy_index = 0; //TODO: maybe const in header?
    clear_vertex(DUMMY_INDEX, pedigree_graph);
}


void dng::PedigreeV2::SimplifyPedigree(dng::Graph &pedigree_graph, const PropEdgeType &edge_types,
                                       const PropEdgeLength &lengths) {

//property_map<Graph, edge_type_t>::type edge_types = get(boost::edge_type, pedigree_graph);
//property_map<Graph, edge_length_t>::type lengths = get(boost::edge_length, pedigree_graph);

//    auto edge_types = get(boost::edge_type, pedigree_graph);
//    auto lengths = get(edge_length, pedigree_graph);
//    auto labels = get(vertex_label, pedigree_graph);
//    auto groups = get(vertex_group, pedigree_graph);
//    auto families = get(edge_family, pedigree_graph);

//    std::cout << edge_types[0] << std::endl;
//    testE(edge_types, pedigree_graph);

    for (vertex_t w = first_library_; w > first_founder_; --w) {
        vertex_t v = w - 1;
        size_t children = 0, ancestors = 0, spouses = 0;

        auto rng = out_edges(v, pedigree_graph);

        for (auto it = rng.first; it != rng.second; ++it) {
            edge_t e = *it;
            if (edge_types[e] == EdgeType::Spousal) {
                spouses += 1;
            } else if (target(e, pedigree_graph) > v) {
                children += 1;
            } else {
                ancestors += 1;
            }
        }
        if (children == 0) {
            // this node has no descendants
            clear_vertex(v, pedigree_graph);
        } else if (children >= 2 || spouses != 0) {
            /*noop*/;
        }
        else {
            edge_t eN[3];// TODO: How "wrong" does the graph has to be to hove ancestor>2? should be impossible (at least logically)
            vertex_t vN[3];//
            for (size_t j = 0; j <= ancestors; ++j) {

                eN[j] = *(rng.first + j);
                vN[j] = target(eN[j], pedigree_graph);
            }
            int child_inedx = ancestors;
            for (size_t j = 0; j < ancestors; ++j) {

                if (vN[child_inedx ] < vN[j]) {
                    boost::swap(vN[child_inedx ], vN[j]);
                    boost::swap(eN[child_inedx ], eN[j]);
                }

                add_edge(vN[j], vN[child_inedx ], {edge_types[eN[j]], lengths[eN[j]] + lengths[eN[child_inedx ]]},
                         pedigree_graph);

            }
            clear_vertex(v, pedigree_graph);
            std::cout << "CleanV2 anc==" << ancestors << " remove: " << v << std::endl;
        }

    }
    PrintEdges("After SimplifyPedigree", pedigree_graph);

}

void dng::PedigreeV2::PrintEdges(std::string prefix, const dng::Graph &pedigree_graph){

//    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    PropVertexIndex index = get(boost::vertex_index, pedigree_graph);

    typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
//    std::pair<edge_iter, edge_iter> ep;
    edge_iter ei2, ei_end2;


    std::cout << prefix << ": V: " << num_vertices(pedigree_graph) << "\tE: " << num_edges(pedigree_graph) <<std::endl;
    for (tie(ei2, ei_end2) = edges(pedigree_graph); ei2 != ei_end2; ++ei2){
        std::cout << "(" << index[source(*ei2, pedigree_graph)]
        << "," << index[target(*ei2, pedigree_graph)] << ") ";
    }
    std::cout <<" ==END==" <<  std::endl;
    std::cout <<"Founder, Non_F, Lib, Somatic: " << first_founder_ << "\t" << first_nonfounder_ << "\t" << first_library_ << "\t" << first_somatic_<< std::endl;
    std::cout << std::endl;


}

void dng::PedigreeV2::EraseRemovedLibraries(dng::ReadGroups &rgs, std::vector<size_t> &node_ids) {

    std::size_t num_libraries = rgs.libraries().size();
    std::vector<std::string> bad_libraries;
    bad_libraries.reserve(num_libraries);
    auto it = rgs.libraries().begin();
    for(size_t u = first_library_; u < node_ids.size(); ++u, ++it) {
        if(node_ids[u] != -1) {
            continue;
        }
        bad_libraries.push_back(*it);
    }
    rgs.EraseLibraries(bad_libraries);


}


void dng::PedigreeV2::ConnectSomaticToLibraries(dng::Graph &pedigree_graph, const ReadGroups &rgs,
                                                const PropVertexLabel &labels){

    const size_t STRLEN_DNG_SM_PREFIX = strlen(DNG_SM_PREFIX);

    for(vertex_t v = (vertex_t) first_somatic_; v < first_library_; ++v) {
        if (labels[v].empty()) {
            continue;
        }
        //using orders from vcf files rgs

        auto r = rgs.data().get<rg::sm>().equal_range(labels[v].c_str() + STRLEN_DNG_SM_PREFIX);
        for (; r.first != r.second; ++r.first) {
            vertex_t w = first_library_ + rg::index(rgs.libraries(), r.first->library);
            if (!edge(v, w, pedigree_graph).second) {
                add_edge(v, w, {EdgeType::Library, 1.0f}, pedigree_graph);
            }
        }
    }
}


void dng::PedigreeV2::CreateFamiliesInfo(dng::Graph &pedigree_graph, family_labels_t &family_labels, std::vector<vertex_t> &pivots){
//
//
//    //LOGIC:: only need 2 variables??
//
//    family_labels_t family_labels(num_families);
//
//    vector<vertex_t> pivots(num_families, dummy_index);
//
//    CreateFamiliesInfo(pedigree_graph, family_labels, pivots);

    auto groups = get(boost::vertex_group, pedigree_graph);
    auto families = get(boost::edge_family, pedigree_graph);


    boost::graph_traits<Graph>::edge_iterator ei, ei_end;


    // Calculate the connected components.  This defines independent sections
    // of the graph.
    std::size_t num_groups = connected_components(pedigree_graph, groups);

    // Calculate the biconnected components and articulation points.
    // This defines "nuclear" families and pivot individuals.
    // Nodes which have no edges will not be part of any family.
    std::vector<vertex_t> articulation_vertices;
    std::size_t num_families = biconnected_components(pedigree_graph, families,
                                                      back_inserter(articulation_vertices)).first;



    family_labels = family_labels_t (num_families);
    pivots = std::vector<vertex_t> (num_families, DUMMY_INDEX);




    // Determine which edges belong to which nuclear families.
//    typedef vector<vector<graph_traits<Graph>::edge_descriptor>> family_labels_t;
//    family_labels_t family_labels(num_families);
    for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
        family_labels[families[*ei]].push_back(*ei);
    }

    // Determine the last family in each group.  All singleton groups will have
    // a value of -1 since they have no family assignment.
    //XXX: strictly size_t != -1, is it good to overflow trick here?
    typedef std::deque<std::size_t> root_families_t;
    root_families_t root_families(num_groups, -1);
    for(std::size_t f = 0; f < family_labels.size(); ++f) {
        // last one wins
        auto first_edge = family_labels[f][0];
        auto src_vertex = source(first_edge, pedigree_graph);
        root_families[groups[src_vertex]] = f;

    }


    // Identify the pivot for each family.
    // The pivot will be the last art. point that has an edge in
    // the group.  The pivot of the last group doesn't matter.
    for(auto a : articulation_vertices) {
        boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
        for(tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
            // Just overwrite existing value so that the last one wins.
            pivots[families[*ei]] = a;
        }
    }


    // Root Pivots are special
    for(auto f : root_families) {
        if(f == -1) { // skip all singleton groups
            continue;
        }
        pivots[f] = DUMMY_INDEX;
    }


//    std::cout << "\n=======================\n" <<
//    "num_goups: " << num_groups << "\tnum_fam: " << num_families << std::endl;
//    for (int j = 0; j < num_vertices(pedigree_graph) ; ++j) {
//        std::cout << "groups["<< j << "] ->" << groups[j]<< std::endl;
//    }
//    for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
//        std::cout << "Edge in family: " << *ei << "\t" << families[*ei] << std::endl;
//    }
//    std::cout << "articulation_vertices: ";
//    for (auto a : articulation_vertices) {
//        std::cout << a << " ";
//    }
//    std::cout << ""<< std::endl;



//    for (int k = 0; k < root_families.size(); ++k) {
//        std::cout << "root_families: [group= " << k << "] -> " <<  root_families[k] << std::endl;
//    }
//
//
//    std::cout << "articulation_vertices: ";
//    for (auto a : articulation_vertices) {
//        std::cout << a << " ";
//    }
//    std::cout << ""<< std::endl;


//    {
//        graph_traits<Graph>::edge_iterator ei, ei_end;
//        for(tie(ei, ei_end) = edges(pedigree_graph); ei != ei_end; ++ei) {
//            std::cout <<"E:     "<< *ei << std::endl;
//        }
//    }
//    {
//        for(auto a : articulation_vertices) {
//            graph_traits<Graph>::in_edge_iterator ei, ei_end;
//            for (tie(ei, ei_end) = in_edges(a, pedigree_graph); ei != ei_end; ++ei) {
//                std::cout << "In_E:  " << *ei << "\t" << a<< std::endl;
//            }
//        }
//    }
//    {
//        for(auto a : articulation_vertices) {
//            graph_traits<Graph>::out_edge_iterator ei, ei_end;
//            for (tie(ei, ei_end) = out_edges(a, pedigree_graph); ei != ei_end; ++ei) {
//                std::cout << "Out_E: " << * ei << "\t" << a<< std::endl;
//            }
//        }
//    }
//    for (int m = 0; m < pivots.size(); ++m) {
//        std::cout << "pivots: [family= " << m << "] -> " << pivots[m] << std::endl;
//    }



    for (int l = 0; l < family_labels.size(); ++l) {
        std::cout << "Families : " << l << "\tPivots: " << pivots[l] << "\t";
        for (auto f : family_labels[l]) {
            std::cout << f << " ";
        }
        std::cout  << std::endl;
    }


}


void dng::PedigreeV2::CreatePeelingOps(dng::Graph &pedigree_graph, const std::vector<size_t> &node_ids,
                                       family_labels_t &family_labels, std::vector<vertex_t> &pivots){

    auto edge_types = get(boost::edge_type, pedigree_graph);
    auto lengths = get(boost::edge_length, pedigree_graph);


    // Resize the information in the pedigree
    transitions_.resize(num_nodes_);
    for(std::size_t i = first_founder_; i < first_nonfounder_; ++i) {
        transitions_[i] = {TransitionType::Founder, static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 0};
    }

//    index = get(vertex_index, pedigree_graph);
    // Detect Family Structure and pivot positions
    for(std::size_t k = 0; k < family_labels.size(); ++k) {
        std::cout <<"\nStart family: " << k << std::endl;
        auto &family_edges = family_labels[k];

        for (auto a : family_edges) {
            std::cout << "edges: " << a << std::endl;
        }
        // Sort edges based on type and target
        boost::sort(family_edges, [&](edge_t x, edge_t y) -> bool { return
                (edge_types(x) < edge_types(y)) &&
                (target(x, pedigree_graph) < target(y, pedigree_graph)); });
        for (auto a : family_edges) {
            std::cout << "Sorted edges: " << a << std::endl;
        }
        // Find the range of the parent types
        auto pos = boost::find_if(family_edges, [&](edge_t x) -> bool {
            return (edge_types(x) != EdgeType::Spousal); });
        size_t num_parent_edges = distance(family_edges.begin(), pos);

        std::cout << "family_edge.size(): " << family_edges.size() <<
        "\tnum_parent_E: " << num_parent_edges << "\t" <<
        "\tpos!=(EdgeType::Spousal): "<< *pos <<
        "\tpivot_V: "<< pivots[k] << std::endl;
        // Check to see what type of graph we have
        if(num_parent_edges == 0) {
            // If we do not have a parent-child single branch,
            // we can't construct the pedigree.
            // TODO: Write Error message
            if(family_edges.size() != 1) {
//                return false;
                throw std::runtime_error("Unable to construct peeler for pedigree;  do not have a parent-child single branch");

            }
            // Create a mitotic peeling operation.
            size_t parent = node_ids[source(*pos, pedigree_graph)];
            size_t child = node_ids[target(*pos, pedigree_graph)];
            std::cout << "=numParentEdge==0: parent " << parent << "\tChild " << child << std::endl;
            TransitionType tt = (edge_types(*pos) == EdgeType::Library) ?
                                TransitionType::Library : TransitionType::Somatic;
            transitions_[child] = {tt, parent, static_cast<size_t>(-1), lengths[*pos], 0};
            std::cout << "===pivots[k]: " << pivots[k] <<
            "\t\tnode_id:" << node_ids[pivots[k]] <<  std::endl;
            family_members_.push_back({parent, child});
            if(node_ids[pivots[k]] == child) {
                peeling_ops_.push_back(peel::op::DOWN);
                std::cout << "=====ADD OP: down" << std::endl;
            } else {
                peeling_ops_.push_back(peel::op::UP);
                std::cout << "=====ADD OP: up" << std::endl;
                if(pivots[k] == DUMMY_INDEX) {
                    roots_.push_back(parent);
                }
            }

        } else if(num_parent_edges == 1) {
            // If this family contains no children, skip it
            if(pos == family_edges.end()) {
                continue;
            }
            // We have a nuclear family with 1 or more children
            size_t dad = node_ids[source(family_edges.front(), pedigree_graph)];
            size_t mom = node_ids[target(family_edges.front(), pedigree_graph)];
            std::cout << "==numParentEdge==1: dad: " << dad << "\tmom: " << mom << std::endl;

            family_members_.push_back({dad, mom});
            auto &family_members = family_members_.back();

            while(pos != family_edges.end()) {
                vertex_t child = node_ids[target(*pos, pedigree_graph)];
                transitions_[child] = { TransitionType::Germline, dad, mom,
                                        lengths[*pos], lengths[*(pos + 1)]
                };
                family_members.push_back(child); // Child
                std::cout << "==Add child to family: " <<child<< std::endl;

                // child edges come in pairs
                ++pos;
                assert(node_ids[target(*pos, pedigree_graph)] == child);
                ++pos;
            }
            if(node_ids[pivots[k]] == node_ids[DUMMY_INDEX]) {
                // A family without a pivot is a root family
                peeling_ops_.push_back(peel::op::TOFATHER);
                std::cout << "=====ADD OP: toFather Root" << std::endl;
                roots_.push_back(family_members[0]);
            } else {
                auto pivot_pos = boost::range::find(family_members, node_ids[pivots[k]]);
                size_t p = distance(family_members.begin(), pivot_pos);
                std::cout << "===pivot_pos: "<< *pivot_pos << "," << node_ids[pivots[k]]<< "\tdistance_pivot_to_family_member_begin(): " << p << std::endl;
                if(p == 0) {
                    peeling_ops_.push_back(peel::op::TOFATHER);
                    std::cout << "=====ADD OP: toFather" << std::endl;
                } else if(p == 1) {
                    peeling_ops_.push_back(peel::op::TOMOTHER);
                    std::cout << "=====ADD OP: toMother" << std::endl;
                } else if(p == 2) {
                    peeling_ops_.push_back(peel::op::TOCHILD);
                    std::cout << "=====ADD OP: toChild" << std::endl;
                } else {
                    peeling_ops_.push_back(peel::op::TOCHILD);
                    std::cout << "=====ADD OP: toChild swap" << std::endl;
                    boost::swap(family_members[p], family_members[2]);
                }
            }

        } else {
            throw std::runtime_error("Unable to construct peeler for pedigree; Not a zero-loop pedigree");
            // TODO: write error message
        }
    }


}