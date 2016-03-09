/*
 * Copyright (c) 2015 Reed A. Cartwright
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
#ifndef DNG_PEELING_H
#define DNG_PEELING_H

#include <vector>
#include <utility>
#include <map>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/fill.hpp>
#include <boost/range/algorithm/fill_n.hpp>

#include <dng/matrix.h>
#include <dng/detail/unit_test.h>
#include <dng/likelihood.h>




#define WORKSPACE_T_MULTIPLE_UPPER_LOWER(workspace, index) \
    (workspace.upper[index] * workspace.lower[index])

#define WORKSPACE_T_KRONECKER_PRODUCT_PARENTS(workspace, dad, mom) \
    kroneckerProduct(WORKSPACE_T_MULTIPLE_UPPER_LOWER(work, dad).matrix(),  \
                     WORKSPACE_T_MULTIPLE_UPPER_LOWER(work, mom).matrix() )


namespace dng {
namespace peel {

namespace op {
enum {
    UP, DOWN, TOFATHER, TOMOTHER, TOCHILD,
    UPFAST, DOWNFAST, TOFATHERFAST, TOMOTHERFAST,
    TOCHILDFAST,
    NUM // Total number of possible forward operations
};

#ifdef DEBUG_PEDIGREE
    std::map<decltype(peel::op::NUM), std::string> map_enum_string {
            {peel::op::UP, "UP"},
            {peel::op::DOWN, "DOWN"},
            {peel::op::TOFATHER, "TOFATHER"},
            {peel::op::TOMOTHER, "TOMOTHER"},
            {peel::op::TOCHILD, "TOCHILD"},
            {peel::op::UPFAST, "UPFAST"},
            {peel::op::DOWNFAST, "DOWNFAST"},
            {peel::op::TOFATHERFAST, "TOFATHERFAST"},
            {peel::op::TOMOTHERFAST, "TOMOTHERFAST"},
            {peel::op::TOCHILDFAST, "TOCHILDFAST"}
    };
#endif

} // namespace dng::peel::op



//XXX: Maybe this should turn into a class, too many functions
    struct workspace_t {
        IndividualVector upper; // Holds P(~Descendent_Data & G=g)
        IndividualVector lower; // Holds P( Descendent_Data | G=g)
        ParentVector super; // Holds P(~Descendent_Data & G=g) for parent nodes

        bool dirty_lower = false;

        double forward_result;

        // Temporary data used by some peeling ops
        PairedGenotypeArray paired_buffer;

        // Information about the pedigree
        std::size_t num_nodes = 0;
        std::pair<std::size_t, std::size_t> founder_nodes;
        std::pair<std::size_t, std::size_t> germline_nodes;
        std::pair<std::size_t, std::size_t> somatic_nodes;
        std::pair<std::size_t, std::size_t> library_nodes;

        // Resize the workspace to fit a pedigree with sz nodes
        void Resize(std::size_t sz) {
            num_nodes = sz;
            upper.resize(num_nodes);
            super.resize(num_nodes);
            lower.assign(num_nodes, DNG_INDIVIDUAL_BUFFER_ONES);
            paired_buffer.resize(100, 1);
            dirty_lower = false;
        }

        // Cleanup after a backwards peeling algorithm
        void Cleanup() {
            boost::fill(lower, DNG_INDIVIDUAL_BUFFER_ONES);
            dirty_lower = false;
        }

        // Cleanup after a backwards peeling algorithm
        // Do not update libraries since they might be set in another operation
        void CleanupFast() {
            // TODO: create a check that sees if this has been done before the
            // forward algorithm.

            boost::fill_n(lower, somatic_nodes.second, DNG_INDIVIDUAL_BUFFER_ONES);
            dirty_lower = false;
        }

        void SetFounders(const GenotypeArray &prior) {
            assert(founder_nodes.first <= founder_nodes.second);
            std::fill(upper.begin() + founder_nodes.first,
                      upper.begin() + founder_nodes.second, prior);
        }

        template<typename T>
        void SetLibraries(const T &range) {
            assert(boost::size(range) == library_nodes.second - library_nodes.first);
            boost::copy(range, lower.begin() + library_nodes.first);
        }

        // calculate genotype likelihoods and store in the lower library vector, maybe/maybe not
        double SetGenotypeLikelihood(
                dng::genotype::DirichletMultinomialMixture &genotype_likelihood,
                const std::vector<depth_t> &depths, std::size_t ref_index) {
            double scale = 0.0, stemp;
            for (std::size_t u = 0; u < depths.size(); ++u) {
                std::tie(lower[library_nodes.first + u], stemp) =
                        genotype_likelihood(depths[u], ref_index);
                scale += stemp;
            }
            return scale;

        }


        inline dng::GenotypeArray multiply_upper_lower(std::size_t index) {

            return (upper[index] * lower[index]);
        }

    };


typedef std::vector<std::size_t> family_members_t;

enum Parents : std::size_t {Father=0, Mother=1};
//enum class Parents : std::size_t {Father=0, Mother=1};


// utility
dng::PairedGenotypeArray sum_over_children(workspace_t &work,
                                           const family_members_t &family,
                                           const TransitionVector &mat);

dng::PairedGenotypeArray sum_over_children(workspace_t &work,
                                           const family_members_t &family,
                                           const TransitionVector &mat,
                                           int first_child_index);


// Core operations
dng::GenotypeArray up_core(workspace_t &work, const family_members_t &family,
                           const TransitionVector &mat);


dng::GenotypeArray to_parent_core(workspace_t &work, const family_members_t &family,
                                  const TransitionVector &mat,
                                  const Parents other_parent);


        // Basic peeling operations
void up(workspace_t &work, const family_members_t &family,
        const TransitionVector &mat);
void down(workspace_t &work, const family_members_t &family,
          const TransitionVector &mat);
void to_father(workspace_t &work, const family_members_t &family,
               const TransitionVector &mat);
void to_mother(workspace_t &work, const family_members_t &family,
               const TransitionVector &mat);
void to_child(workspace_t &work, const family_members_t &family,
              const TransitionVector &mat);

// Fast versions that can be used on "dirty" workspaces
void up_fast(workspace_t &work, const family_members_t &family,
             const TransitionVector &mat);
void down_fast(workspace_t &work, const family_members_t &family,
               const TransitionVector &mat);
void to_father_fast(workspace_t &work, const family_members_t &family,
                    const TransitionVector &mat);
void to_mother_fast(workspace_t &work, const family_members_t &family,
                    const TransitionVector &mat);
void to_child_fast(workspace_t &work, const family_members_t &family,
                   const TransitionVector &mat);

// Reverse versions that can be used for backwards algorithms
void up_reverse(workspace_t &work, const family_members_t &family,
                const TransitionVector &mat);
void down_reverse(workspace_t &work, const family_members_t &family,
                  const TransitionVector &mat);
void to_father_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionVector &mat);
void to_mother_reverse(workspace_t &work, const family_members_t &family,
                       const TransitionVector &mat);
void to_child_reverse(workspace_t &work, const family_members_t &family,
                      const TransitionVector &mat);




typedef decltype(&down) function_t;

struct info_t {
    bool writes_lower;
    int writes_to;
};

constexpr info_t info[op::NUM] = {
    /* Up           */ {true,  0},
    /* Down         */ {false, 1},
    /* ToFather     */ {true,  0},
    /* ToMother     */ {true,  1},
    /* ToChild      */ {false, 2},
    /* UpFast       */ {true,  0},
    /* DownFast     */ {false, 1},
    /* ToFatherFast */ {true,  0},
    /* ToMotherFast */ {true,  1},
    /* ToChildFast  */ {false, 2}
};

// TODO: Write test case to check that peeling ops are in the right order.
constexpr function_t functions[op::NUM] = {
    &up, &down, &to_father, &to_mother, &to_child,
    &up_fast, &down_fast, &to_father_fast, &to_mother_fast,
    &to_child_fast
};

constexpr function_t reverse_functions[op::NUM] = {
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse,
    &up_reverse, &down_reverse, &to_father_reverse,
    &to_mother_reverse, &to_child_reverse
};

} // namespace dng::peel
} // namespace dng

#endif // DNG_PEELING_H
