//
// Created by steven on 2/22/16.
//


#include <dng/find_mutation_x.h>


FindMutationsXLinked::FindMutationsXLinked(const RelationshipGraph &ship_graph,
                                           const FindMutationParams &params)
        : AbstractFindMutations(ship_graph, params) {

    using namespace dng;

    for (size_t child = 0; child < work_nomut_.num_nodes; ++child) {
        auto trans = ship_graph.transitions()[child];
        std::cout << "Node:" << child << "\t" <<
                "\tType:" << (int) trans.type <<
                "\tSex:" << (int) trans.sex <<
                //                trans.length1 << "\t" << trans.length2 << "\t" <<
                "\tP: " <<trans.parent1 << "_" << trans.parent2 << "\t" <<
                std::endl;

        if (trans.type == RelationshipGraph::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
        } else if (trans.type == RelationshipGraph::TransitionType::Somatic ||
                   trans.type == RelationshipGraph::TransitionType::Library) {
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
        } else {
            full_transition_matrices_[child] = {};
            nomut_transition_matrices_[child] = {};
            posmut_transition_matrices_[child] = {};
            onemut_transition_matrices_[child] = {};
            mean_matrices_[child] = {};
        }

    }

}


void FindMutationsXLinked::run(const std::vector<depth_t> &depths,
                               int ref_index, MutationStats &mutation_stats) {

    double scale = work_nomut_.SetGenotypeLikelihood(genotype_likelihood_, depths,
                                                     ref_index);
    // Set the prior probability of the founders given the reference
    work_nomut_.SetFounders(genotype_prior_[ref_index]);
    work_full_ = work_nomut_; //TODO: full test on copy assignment operator

    bool is_mup_less_threshold = CalculateMutationProb(mutation_stats);
    if (is_mup_less_threshold) {

    }
    ship_graph.PeelBackwards(work_full_, full_transition_matrices_);

    mutation_stats.SetScaledLogLikelihood(scale);
    mutation_stats.SetGenotypeLikelihoods(work_full_, depths.size());
    mutation_stats.SetPosteriorProbabilities(work_full_);

    mutation_stats.CalculateExpectedMutation(work_full_, mean_matrices_);
    mutation_stats.CalculateNodeMutation(work_full_, posmut_transition_matrices_);
    CalculateDenovoMutation(mutation_stats);

}



bool AbstractFindMutations::CalculateMutationProb(MutationStats &mutation_stats) {

    // Calculate log P(Data, nomut ; model)
    ship_graph.PeelForwards(work_nomut_, nomut_transition_matrices_);

    /**** Forward-Backwards with full-mutation ****/
    // Calculate log P(Data ; model)
    ship_graph.PeelForwards(work_full_, full_transition_matrices_);

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    bool is_mup_less_threshold = mutation_stats.CalculateMutationProb(work_nomut_,
                                                                      work_full_);

    return is_mup_less_threshold;
}


void AbstractFindMutations::CalculateDenovoMutation(MutationStats &mutation_stats) {

    ship_graph.PeelBackwards(work_nomut_, nomut_transition_matrices_);
//    mutation_stats.CalculateDenovoMutation(work_nomut_, onemut_transition_matrices_,
//                                           ship_graph);

}

