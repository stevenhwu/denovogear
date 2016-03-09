//
// Created by steven on 2/22/16.
//


#include <dng/find_mutation_x.h>

FindMutationsXLinked::FindMutationsXLinked(const Pedigree &pedigree,
                                           FindMutationParams &params)
        : AbstractFindMutations(pedigree, params) {

    using namespace dng;

//    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
//    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
//    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
//    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
//    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    for(size_t child = 0; child < work_nomut_.num_nodes; ++child) {

        auto trans = pedigree.transitions()[child];
//        std::cout << "Node:" << child << "\t" << "\t" << trans.length1 << "\t" << trans.length2 <<
//                "\tType:" << (int) trans.type << std::endl;

        if(trans.type == Pedigree::TransitionType::Germline) {
            auto dad = f81::matrix(trans.length1, params_.nuc_freq);
            auto mom = f81::matrix(trans.length2, params_.nuc_freq);

            full_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom);
            nomut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = meiosis_diploid_matrix(dad, mom, 1);
            mean_matrices_[child] = meiosis_diploid_mean_matrix(dad, mom);
        } else if(trans.type == Pedigree::TransitionType::Somatic ||
                  trans.type == Pedigree::TransitionType::Library) {
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

    double scale = work_nomut_.SetGenotypeLikelihood(genotype_likelihood_, depths, ref_index);

    // Set the prior probability of the founders given the reference
    work_nomut_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data, nomut ; model)
//        const double logdata_nomut = pedigree_.PeelForwards(work_, nomut_transition_matrices_);

    /**** Forward-Backwards with full-mutation ****/
    // Calculate log P(Data ; model)
    const double logdata = pedigree_.PeelForwards(work_nomut_, full_transition_matrices_);


    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
//        bool is_mup_lt_threshold = mutation_stats.set_mutation_prob(logdata_nomut, logdata);

//
//        bool is_mup_lt_threshold = calculate_mutation_prob(mutation_stats);
//        if (is_mup_lt_threshold) {
//            return false;
//        }
//
//        mutation_stats.set_scaled_log_likelihood(scale);
//        mutation_stats.set_genotype_likelihood(work_, depths.size());
}
