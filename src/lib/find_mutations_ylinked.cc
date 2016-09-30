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
#include <iostream>
#include <dng/depths.h>
using namespace dng;

FindMutationsYLinked::~FindMutationsYLinked() {
	// TODO Auto-generated destructor stub
}


FindMutationsYLinked::FindMutationsYLinked(double min_prob,
        const RelationshipGraph &graph, params_t params)
        : FindMutationsAbstract(min_prob, graph, params)
//    relationship_graph_(graph), min_prob_(min_prob),
//    params_(params), genotype_likelihood_{params.params_a, params.params_b},
//    work_full_(graph.CreateWorkspace()), work_nomut_(graph.CreateWorkspace())
    {
    std::cout << "FMY:FMY" << std::endl;

    SetupTransitionMatrix();

//    // Extract relevant subsets of matrices
////    for(
//    dng::TransitionVector transition_matrices_;
//[[deprecated]] size_t color = 0;
////            color < dng::pileup::AlleleDepths::type_info_table_length; ++color) {
//        int width = 4;//dng::pileup::AlleleDepths::type_info_gt_table.width;
//        // Resize our subsets to the right width
//        transition_matrices_.resize(full_transition_matrices_.size());
//        dng::TransitionVector transition_matrices_2 (3);
//        dng::TransitionVector transition_matrices_3 {3};
////        explicit dng::TransitionVector transition_matrices_4 {3};
////        std::cout << transition_matrices_.size()
////                << transition_matrices_2.size()
////                << transition_matrices_3.size()
//////                << transition_matrices_4.size()
////                << std::endl; //233
//
//        // enumerate over all children
//        for(size_t child = 0; child < work_full_.num_nodes; ++child) {
//            auto trans = relationship_graph_.transitions()[child];
//            if(trans.type == RelationshipGraph::TransitionType::Germline) {
//                // resize transition matrix to w*w,w
//                transition_matrices_[child].resize(width*width,width);
//                // Assume column major order which is the default
//                for(int a = 0; a < width; ++a) {
//                    int ga = dng::pileup::AlleleDepths::type_info_gt_table[40].indexes[a];
//                    for(int b = 0; b < width; ++b) {
//                        int gb = dng::pileup::AlleleDepths::type_info_gt_table[40].indexes[b];
//                        // copy correct value from the full matrix to the subset matrix
//                        int x = a*width+b;
//                        int gx = ga*10+gb;
//                        for(int y = 0; y < width; ++y) {
//                            int gy = dng::pileup::AlleleDepths::type_info_gt_table[40].indexes[y];
//                            // copy correct value from the full matrix to the subset matrix
//                            transition_matrices_[child](x,y) = full_transition_matrices_[child](gx,gy);
//                        }
//                    }
//                }
//            } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
//                      trans.type == RelationshipGraph::TransitionType::Library) {
//                transition_matrices_[child].resize(SIZE4, SIZE4);
//                // Assume column major order which is the default
//                int cc = 0;
//                for (int var = 0; var < 10; ++var) {
//                    for (int var2 = 0; var2 < 10; ++var2) {
//                        full_transition_matrices_[child](var,var2) = cc++;
//                    }
//
//                }
//
//                for(int x = 0; x < SIZE4; ++x) {
//                    int gx = MAP_4_TO_10[x];
//                    for(int y = 0; y < SIZE4; ++y) {
//                        int gy = MAP_4_TO_10[y];
//                        full_transition_matrices_[child](x,y) = full_transition_matrices_[child](gx,gy);
//                    }
//                }
//                std::cout << full_transition_matrices_[child] << "\n\n" << std::endl;
//                std::cout << full_transition_matrices_[child].block(0,0,4,4) << "\n\n" << std::endl;
//                std::cout << full_transition_matrices_[child].topLeftCorner(4,4) << "\n\n" << std::endl;
//                std::cout << transition_matrices_[child] << "\n\n" << std::endl;
////                full_transition_matrices_[child] = full_transition_matrices_[child].topLeftCorner(4,4);
////                full_transition_matrices_[child] = topLeftCorner(full_transition_matrices_[child], 4,4);
//                full_transition_matrices_[child] = full_transition_matrices_[child].block<4,4>(0,0).eval() ;
//                std::cout << full_transition_matrices_[child] << "\n\n" << std::endl;
//            } else {
//                transition_matrices_[child] = {};
//            }
//        }
////    }

//        std::exit(93);
        for (int i = 0; i < 5; ++i) {
            Resize10To4(genotype_prior_[i]);
        }
//    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq,
//                                             {params_.ref_weight, 0, 0, 0});
//       genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq,
//                                             {0, params_.ref_weight, 0, 0});
//       genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq,
//                                             {0, 0, params_.ref_weight, 0});
//       genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq,
//                                             {0, 0, 0, params_.ref_weight});
//       genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq,
//                                             {0, 0, 0, 0});




#if CALCULATE_ENTROPY == 1
    std::cerr << "Entropy will not be calculated for Y-linked model!!" << std::endl;
#endif

    keep_library_index_ = graph.keep_library_index();
}

// Returns true if a mutation was found and the record was modified
bool FindMutationsYLinked::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;

    std::cout << "FMY()" << std::endl;
    assert(stats != nullptr);

    //TODO(SW): Eventually, MutationStats this will replace all stats_t
    MutationStats mutation_stats(min_prob_);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    std::cout << "DepthsSize: " << depths.size() << std::endl;
//    for(std::size_t u = 0; u < keep_library_index_.size(); ++u) {
//            std::cout << u << "\t" << work_full_.library_nodes.first + u
//                    << "\t" << keep_library_index_[u] << std::endl;
//            for (int i = 0; i < 4; ++i) {
////               std::cout << depths[keep_library_index_[u]].counts[i] << " ";
//            }
//
//            std::cout << "" << std::endl;
//            std::tie(work_full_.lower[work_full_.library_nodes.first + u], stemp) =
//                genotype_likelihood_(depths[keep_library_index_[u]], ref_index);
//            scale += stemp;
//    }
//    std::cout << "\n" << std::endl;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::cout << u << "\tLowerIndex: " << work_full_.library_nodes.first + u << "\t" <<   std::endl;
        for (int i = 0; i < 4; ++i) {
           std::cout << depths[u].counts[i] << " ";
        } std::cout << "" << std::endl;

        std::tie(work_full_.lower[work_full_.library_nodes.first + u], stemp) =
            genotype_likelihood_(depths[u], ref_index);

        Resize10To4(work_full_.lower[work_full_.library_nodes.first + u]);

        scale += stemp;
        std::cout.precision(10);
        std::cout << "Genotype:" << u << " "
                << work_full_.library_nodes.first + u << "\n"
//                << work_full_.lower[work_full_.library_nodes.first + u]
                << std::endl;


    }
//
    // Set the prior probability of the founders given the reference
    work_full_.SetFounders(genotype_prior_[ref_index]);
    work_nomut_ = work_full_;

    bool is_mup_less_threshold = CalculateMutationProb(mutation_stats);

    if (is_mup_less_threshold) {
        return false;
    }
//    std::exit(133);
//    relationship_graph_.PeelBackwards(work_full_, full_transition_matrices_);
//
//    mutation_stats.SetGenotypeLikelihoods(work_full_, depths.size());
//    mutation_stats.SetScaledLogLikelihood(scale);
//    mutation_stats.SetPosteriorProbabilities(work_full_);
//
//    //PR_NOTE(SW): Reassign mutation_stasts back to stats_t. This avoid major changes in call.cc at this stage
//    //Remove this section after all PR are completed
//    stats-> mup = mutation_stats.mup_;
//    double pmut = stats->mup; //rest of the code will work
//    stats-> lld = mutation_stats.lld_;
////    stats-> llh = mutation_stats.llh_;
//    stats-> genotype_likelihoods = mutation_stats.genotype_likelihoods_;
//    stats-> posterior_probabilities = mutation_stats.posterior_probabilities_;
//    //End Remove this section
//
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
//
//    /**** Forward-Backwards with no-mutation ****/
//
//    // TODO: Better to use a separate workspace???
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
//
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
    std::cout << "FMY:SetupTransitionMatrix" << std::endl;
    for(size_t child = 0; child < work_full_.num_nodes; ++child) {
        auto trans = relationship_graph_.transitions()[child];
        std::cout << "Child: " << child << "\tTransType: " << (int) trans.type << std::endl;

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
//        else { //if(trans.type == RelationshipGraph::TransitionType::Germline) {
//            full_transition_matrices_[child] = {};
//            nomut_transition_matrices_[child] = {};
//            posmut_transition_matrices_[child] = {};
//            onemut_transition_matrices_[child] = {};
//            mean_matrices_[child] = {};
//        }
    }


    //

    //Original below
//    for(size_t child = 0; child < work_full_.num_nodes; ++child) {
//        auto trans = relationship_graph_.transitions()[child];
//        std::cout << "Child: " << child << "\tTransType: " << (int) trans.type << std::endl;
//        if(trans.type == RelationshipGraph::TransitionType::Germline) {
//            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
//            auto mom = f81::matrix(trans.length2, params_.nuc_freq);
//
//            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
//
//            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
//            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
//                                                 nomut_transition_matrices_[child];
//            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
//
//            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
//        } else if(trans.type == RelationshipGraph::TransitionType::Somatic ||
//                  trans.type == RelationshipGraph::TransitionType::Library) {
//            auto orig = f81::matrix(trans.length1, params_.nuc_freq);
//            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
//            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
//            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
//                                                 nomut_transition_matrices_[child];
//            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
//            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
//        } else {
////            full_transition_matrices_[child] = {};
////            nomut_transition_matrices_[child] = {};
////            posmut_transition_matrices_[child] = {};
////            onemut_transition_matrices_[child] = {};
////            mean_matrices_[child] = {};
//        }
//    }

}
