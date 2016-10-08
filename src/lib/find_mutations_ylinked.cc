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


#include <dng/find_mutations_ylinked.h>


using namespace dng;

FindMutationsYLinked::~FindMutationsYLinked() {
	// TODO Auto-generated destructor stub
}


FindMutationsYLinked::FindMutationsYLinked(double min_prob,
        const RelationshipGraph &graph, params_t params)
        : FindMutationsAbstract(min_prob, graph, params)
    {

    SetupPopulationPriorHaploid();
    SetupTransitionMatrix();

#if CALCULATE_ENTROPY == 1
    std::cerr << "Entropy will not be calculated for Y-linked model!!" << std::endl;
#endif

}

// Returns true if a mutation was found and the record was modified
bool FindMutationsYLinked::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;

    assert(stats != nullptr);

    //TODO(SW): Eventually, MutationStats this will replace all stats_t
    MutationStats mutation_stats(min_prob_);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_full_.lower[work_full_.library_nodes.first + u], stemp) =
                genotype_likelihood_.OperateHaploid(depths[u], ref_index);
        scale += stemp;
    }
    // Set the prior probability of the founders given the reference
    work_full_.SetFounders(genotype_prior_[ref_index]);
    work_nomut_ = work_full_;

    bool is_mup_less_threshold = CalculateMutationProb(mutation_stats);

    //TODO(SW): HACK: for unittest
    stats-> mup = mutation_stats.mup_;
    if (is_mup_less_threshold) {
        return false;
    }

    relationship_graph_.PeelBackwards(work_full_, full_transition_matrices_);

    mutation_stats.SetGenotypeLikelihoods(work_full_, depths.size());
    mutation_stats.SetScaledLogLikelihood(scale);
    mutation_stats.SetPosteriorProbabilities(work_full_);
//
    //PR_NOTE(SW): Reassign mutation_stasts back to stats_t. This avoid major changes in call.cc at this stage
    //Remove this section after all PR are completed
    stats-> mup = mutation_stats.mup_;
    double pmut = stats->mup; //rest of the code will work
    stats-> lld = mutation_stats.lld_;
//    stats-> llh = mutation_stats.llh_;
    stats-> genotype_likelihoods = mutation_stats.genotype_likelihoods_;
    stats-> posterior_probabilities = mutation_stats.posterior_probabilities_;
    //End Remove this section


    //PR_NOTE(SW): Rest of the ylinked model will be add back after unittest is completed
//    double mux = 0.0;
//    event_.assign(work_full_.num_nodes, 0.0);
//    for(std::size_t i = work_full_.founder_nodes.second; i < work_full_.num_nodes; ++i) {
//        mux += (work_full_.super[i] * (mean_matrices_[i] *
//                                  work_full_.lower[i].matrix()).array()).sum();
//        event_[i] = (work_full_.super[i] * (posmut_transition_matrices_[i] *
//                                       work_full_.lower[i].matrix()).array()).sum();
//        event_[i] = event_[i] / pmut;
//    }
//    stats->mux = mux;
//
//    stats->node_mup.resize(work_full_.num_nodes, hts::bcf::float_missing);
//    for(size_t i = work_full_.founder_nodes.second; i < work_full_.num_nodes; ++i) {
//        stats->node_mup[i] = static_cast<float>(event_[i]);
//    }
////
////    /**** Forward-Backwards with no-mutation ****/
////
////    // TODO: Better to use a separate workspace???
//    relationship_graph_.PeelForwards(work_nomut_, nomut_transition_matrices_);
//    relationship_graph_.PeelBackwards(work_nomut_, nomut_transition_matrices_);
//    event_.assign(work_nomut_.num_nodes, 0.0);
//    double total = 0.0, entropy = 0.0, max_coeff = -1.0;
//    size_t dn_loc = 0, dn_col = 0, dn_row = 0;
//    for(std::size_t i = work_nomut_.founder_nodes.second; i < work_nomut_.num_nodes; ++i) {
//		Eigen::ArrayXXd mat = (work_nomut_.super[i].matrix()
//				* work_nomut_.lower[i].matrix().transpose()).array()
//				* onemut_transition_matrices_[i].array();
//
//        std::size_t row, col;
//        double mat_max = mat.maxCoeff(&row, &col);
//        if(mat_max > max_coeff) {
//            max_coeff = mat_max;
//            dn_row  = row;
//            dn_col = col;
//            dn_loc = i;
//        }
//        event_[i] = mat.sum();
//        entropy += (mat.array() == 0.0).select(mat.array(),
//                                               mat.array() * mat.log()).sum();
//        total += event_[i];
//    }
//    // Calculate P(only 1 mutation)
//    const double pmut1 = total * (1.0 - pmut);
//    stats->mu1p = pmut1;
////
//    // Output statistics for single mutation only if it is likely
//    if(pmut1 / pmut >= min_prob_) {
//        stats->has_single_mut = true;
//        for(std::size_t i = work_nomut_.founder_nodes.second; i < work_nomut_.num_nodes; ++i) {
//            event_[i] = event_[i] / total;
//        }
//
//        // Calculate entropy of mutation location
//        entropy = (-entropy / total + log(total)) / M_LN2;
//        entropy /= max_entropies_[ref_index];
//        stats->dnc = std::round(100.0 * (1.0 - entropy));
//
//        stats->dnq = lphred<int32_t>(1.0 - (max_coeff / total), 255);
//        stats->dnl = relationship_graph_.labels()[dn_loc];
//        if(relationship_graph_.transitions()[dn_loc].type == RelationshipGraph::TransitionType::Germline) {
//            stats->dnt = &meiotic_diploid_mutation_labels[dn_row][dn_col][0];
//        } else {
//            stats->dnt = &mitotic_diploid_mutation_labels[dn_row][dn_col][0];
//        }
//
//        stats->node_mu1p.resize(work_nomut_.num_nodes, hts::bcf::float_missing);
//        for(size_t i = work_nomut_.founder_nodes.second; i < work_nomut_.num_nodes; ++i) {
//            stats->node_mu1p[i] = static_cast<float>(event_[i]);
//        }
//    } else {
//        stats->has_single_mut = false;
//    }
    return true;
}



void FindMutationsYLinked::SetupTransitionMatrix(){

    for(size_t child = 0; child < work_full_.num_nodes; ++child) {
        auto trans = relationship_graph_.transitions()[child];

        if (trans.type == RelationshipGraph::TransitionType::Somatic
                || trans.type == RelationshipGraph::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);
            full_transition_matrices_[child] = mitosis_haploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_haploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_haploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_haploid_mean_matrix(orig);
        }
    }


}
