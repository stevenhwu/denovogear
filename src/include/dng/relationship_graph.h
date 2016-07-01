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

#pragma once
#ifndef DNG_PEDIGREE_H
#define DNG_PEDIGREE_H

#include <functional>
#include <cmath>
#include <array>

#include <dng/matrix.h>
#include <dng/io/ped.h>
#include <dng/newick.h>
#include <dng/read_group.h>
#include <dng/peeling.h>
#include <dng/detail/unit_test.h>

#define DEBUG_VERBOSE 1

namespace dng {

class RelationshipGraph {
public:

    typedef boost::property_map<Graph, boost::edge_type_t>::type PropEdgeType;
    typedef boost::property_map<Graph, boost::edge_length_t>::type PropEdgeLength;
    typedef boost::property_map<Graph, boost::vertex_label_t>::type PropVertexLabel;
    typedef boost::property_map<Graph, boost::vertex_group_t>::type PropVertexGroup;
    typedef boost::property_map<Graph, boost::vertex_index_t>::type PropVertexIndex;

    typedef boost::property_map<Graph, boost::vertex_index_t>::type IndexMap;

    typedef std::vector<std::vector<boost::graph_traits<dng::Graph>::edge_descriptor>> family_labels_t;


    enum class InheritancePattern : int {
        AUTOSOMAL = 0,
        DEFAULT = 0,
        MATERNAL = 1,
        PATERNAL = 2,
        X_LINKED = 3,
        Y_LINKED = 4,
        W_LINKED = 5,
        Z_LINKED = 6,

//        autosomal (the default)
//        xlinked (females have 2 copies, males have 1; males transmit to daughters, not to sons)
//        ylinked (males have 1 copy, only transmits it to sons)
//        wlinked (females have 1 copy, only transmited to daughters)
//        zlinked (males have 2 copies, females have 1; females transmit to sons, not to daughters)
//        maternal (transmitted by mother to child)
//        paternal (transmitter by father to child)
    };

    //TODO: struct FamilyInfo/Family structure.
    //Op1: A struct to record info in each family. family_t and ops
    //Op2: Another struct to group families together, include pivots and root?
    enum class FamilyType : int {
        PAIR = 2,
        TRIO = 3
    };

    enum class TransitionType {
        Founder, Germline, Somatic, Library
    };

    struct FamilyInfo {
        FamilyType family_type;
        family_labels_t family_labels;//(num_families);
        std::vector<vertex_t> pivots;//(num_families, dummy_index);
    };

    struct transition_t {
        TransitionType type;
        std::size_t parent1;
        std::size_t parent2;
        double length1;
        double length2;
    };

    bool Construct(const io::Pedigree &pedigree, dng::ReadGroups &rgs,
                   double mu, double mu_somatic, double mu_library);

    //TODO(SW): Eventually replace with this, or pass inheritance with a different method
    bool Construct(const io::Pedigree &pedigree, dng::ReadGroups &rgs,
            const InheritancePattern &pattern, double mu, double mu_somatic,
            double mu_library);


    double PeelForwards(peel::workspace_t &work,
                        const TransitionVector &mat) const {
        if(work.dirty_lower) {
            work.CleanupFast();
        }
        // Peel pedigree one family at a time
        for(std::size_t i = 0; i < peeling_functions_.size(); ++i) {
            (*peeling_functions_[i])(work, family_members_[i], mat);
        }
        // Sum over roots
        double ret = 0.0;
        for(auto r : roots_) {
            ret += log((work.lower[r] * work.upper[r]).sum());
        }
        return ret;
    }

    double PeelBackwards(peel::workspace_t &work,
                         const TransitionVector &mat) const {
        double ret = 0.0;
        // Divide by the log-likelihood
        for(auto r : roots_) {
            double sum = (work.lower[r] * work.upper[r]).sum();
            ret += log(sum);
            sum = sqrt(sum);
            work.lower[r] /= sum;
            work.upper[r] /= sum;
        }

        for(std::size_t i = peeling_reverse_functions_.size(); i > 0; --i) {
            (*peeling_reverse_functions_[i - 1])(work, family_members_[i - 1], mat);
        }
        work.dirty_lower = true;
        return ret;
    }

    peel::workspace_t CreateWorkspace() const {
        peel::workspace_t work;
        work.Resize(num_nodes_);
        work.founder_nodes = std::make_pair(first_founder_, first_nonfounder_);
        work.germline_nodes = std::make_pair(first_founder_, first_somatic_);
        work.somatic_nodes = std::make_pair(first_somatic_, first_library_);
        work.library_nodes = std::make_pair(first_library_, num_nodes_);
        return work;
    }

    void PrintMachine(std::ostream &os);
    void PrintTable(std::ostream &os);
    void PrintStates(std::ostream &os, double scale = 0.0);

    std::vector<std::string> BCFHeaderLines() const;

    const std::vector<transition_t> &transitions() const { return transitions_; }

    const std::vector<std::string> &labels() const { return labels_; }

    size_t num_nodes() const { return num_nodes_; }
    std::pair<size_t, size_t> library_nodes() const { return {first_library_, num_nodes_}; }


protected:

    // node structure:
    // founder germline, non-founder germline, somatic, library
    std::size_t num_nodes_{0};        // total number of nodes
    std::size_t first_founder_{0};    // start of founder germline
    std::size_t first_nonfounder_{0}; // start of non-founder germline
    std::size_t first_somatic_{0};    // start of somatic nodes
    std::size_t first_library_{0};    // start of libraries

    std::vector<std::size_t> roots_;

    // Pedigree Structure
    std::vector<std::string> labels_;
    std::vector<transition_t> transitions_;

    // The original, simplified peeling operations
    std::vector<decltype(peel::op::NUM)> peeling_ops_;
    // The modified, "faster" operations
    std::vector<decltype(peel::op::NUM)> peeling_functions_ops_;
    // Array of functions that will be called to perform the peeling
    std::vector<peel::function_t> peeling_functions_;
    std::vector<peel::function_t> peeling_reverse_functions_;

    // The arguments to a peeling operation
    std::vector<peel::family_members_t> family_members_;

    void ConstructPeelingMachine();

//PR_NOTE(SW): New functions here
    void SetupFirstNodeIndex(const io::Pedigree &pedigree);

    void ParseIoPedigree(dng::Graph &pedigree_graph,
            const dng::io::Pedigree &pedigree);

    void AddLibrariesFromReadGroups(dng::Graph &pedigree_graph,
            const dng::ReadGroups &rgs, PropVertexLabel &labels);

    void ConnectSomaticToLibraries(dng::Graph &pedigree_graph,
            const ReadGroups &rgs, const PropVertexLabel &labels);

    void UpdateEdgeLengths(dng::Graph &pedigree_graph, double mu,
            double mu_somatic, double mu_library);

    void SimplifyPedigree(dng::Graph &pedigree_graph,
            const PropEdgeType &edge_types, const PropEdgeLength &lengths);

    void UpdateLabelsNodeIds(dng::Graph &pedigree_graph, dng::ReadGroups rgs,
            std::vector<size_t> &node_ids);

    void EraseRemovedLibraries(dng::ReadGroups &rgs,
            std::vector<size_t> &node_ids);

    void CreateFamiliesInfo(dng::Graph &pedigree_graph,
            family_labels_t &family_labels, std::vector<vertex_t> &pivots);

    void CreatePeelingOps(dng::Graph &pedigree_graph,
            const std::vector<size_t> &node_ids, family_labels_t &family_labels,
            std::vector<vertex_t> &pivots);

private:
    void PrintDebugEdges(const std::string &prefix,
            const dng::Graph &pedigree_graph);

    const vertex_t DUMMY_INDEX = 0;

    DNG_UNIT_TEST(test_pedigree_inspect);
    DNG_UNIT_TEST(test_parse_io_pedigree);
    DNG_UNIT_TEST(test_add_lib_from_rgs);
};

}; // namespace dng


#endif // DNG_PEDIGREE_H
