//
// Created by steven on 2/3/16.
//
#include <iostream>


#include <dng/mutation_stats.h>



MutationStats::MutationStats(double min_prob) : min_prob_(min_prob) { }


bool MutationStats::CalculateMutationProb(const dng::peel::workspace_t &work_nomut,
                                          const dng::peel::workspace_t &work_full) {
    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    logdata_nomut_ = work_nomut.forward_result;
    logdata_ = work_full.forward_result;
    mup_ = -std::expm1(logdata_nomut_ - logdata_) ;
    return mup_ < min_prob_;
}


void MutationStats::SetScaledLogLikelihood(double scale) {
    lld_ = (logdata_ + scale) / M_LN10;
    llh_ = logdata_ / M_LN10;

}

void MutationStats::SetGenotypeLikelihoods(
        const dng::peel::workspace_t &workspace, const int depth_size) {

    genotype_likelihoods_.resize(workspace.num_nodes);
    for(std::size_t u = 0; u < depth_size; ++u) {
        std::size_t pos = workspace.library_nodes.first + u;
        genotype_likelihoods_[pos] = workspace.lower[pos].log() / M_LN10;
        //TODO: Eigen 3.3 might have log10()
    }
}



void MutationStats::SetPosteriorProbabilities(
        const dng::peel::workspace_t &workspace) {

    posterior_probabilities_.resize(workspace.num_nodes);
    for(std::size_t i = 0; i < workspace.num_nodes; ++i) {
        posterior_probabilities_[i] = WORKSPACE_T_MULTIPLE_UPPER_LOWER(
                workspace, i);
        posterior_probabilities_[i] /= posterior_probabilities_[i].sum();
    }

}

void MutationStats::SetExactlyOneMutation(double total){
    mu1p_ = total * (1.0 - mup_);
    has_single_mut_ = (mu1p_ / mup_) >= min_prob_;
}



void MutationStats::CalculateExpectedMutation(dng::peel::workspace_t &work_full,
                                              dng::TransitionVector &mean_matrices){
    mux_ = 0.0;
    for(size_t i = work_full.founder_nodes.second; i < work_full.num_nodes; ++i) {
        mux_ += (work_full.super[i] * (mean_matrices[i] *
                                      work_full.lower[i].matrix()).array()).sum();
    }

};

void MutationStats::CalculateNodeMutation(dng::peel::workspace_t &work_full,
                                          dng::TransitionVector &posmut_transition_matrices) {
    std::vector<double> event;
    event.assign(work_full.num_nodes, 0.0);
    for (size_t i = work_full.founder_nodes.second; i < work_full.num_nodes; ++i) {
        event[i] = (work_full.super[i] * (posmut_transition_matrices[i] *
                                          work_full.lower[i].matrix()).array()).sum();
        event[i] = event[i] / mup_;
    }

    SetNodeMup(event, work_full.founder_nodes.second);


}

void MutationStats::CalculateDenovoMutation(dng::peel::workspace_t &work_nomut,
                                            dng::TransitionVector &onemut_transition_matrices,
                                            const dng::Pedigree &pedigree_) {
    std::vector<double> event;
    event.assign(work_nomut.num_nodes, 0.0);
    double total = 0.0, entropy = 0.0, max_coeff = -1.0;
    size_t dn_loc = 0, dn_col = 0, dn_row = 0;

    for (std::size_t i = work_nomut.founder_nodes.second; i < work_nomut.num_nodes; ++i) {
        Eigen::ArrayXXd mat = (work_nomut.super[i].matrix() *
                               work_nomut.lower[i].matrix().transpose()).array() *
                              onemut_transition_matrices[i].array();

        std::size_t row, col;
        double mat_max = mat.maxCoeff(&row, &col);
        if (mat_max > max_coeff) {
            max_coeff = mat_max;
            dn_row = row;
            dn_col = col;
            dn_loc = i;
        }
        event[i] = mat.sum();
        total += event[i];

#if CALCULATE_ENTROPY == true
        entropy += (mat.array() == 0.0).select(mat.array(),
                                                       mat.array() * mat.log()).sum();
#endif

    }
    SetExactlyOneMutation(total);

    if (has_single_mut_) {
        SetNodeMu1p(event, total, work_nomut.founder_nodes.second);

        dnq_ = dng::util::lphred<int32_t>(1.0 - (max_coeff / total), 255);
        dnl_ = pedigree_.labels()[dn_loc];
        if (pedigree_.transitions()[dn_loc].type == dng::Pedigree::TransitionType::Germline) {
            dnt_ = &dng::meiotic_diploid_mutation_labels[dn_row][dn_col][0];
        } else {
            dnt_ = &dng::mitotic_diploid_mutation_labels[dn_row][dn_col][0];
        }

#if CALCULATE_ENTROPY == true
        // Calculate entropy of mutation location
                entropy = (-entropy / total + log(total)) / M_LN2;
                entropy /= max_entropies_[ref_index];
                mutation_stats.dnc = std::round(100.0 * (1.0 - entropy));
#endif

    }
}




void MutationStats::SetGenotypeRelatedStats(const int (&acgt_to_refalt_allele)[5],
                                            const int (&refalt_to_acgt_allele)[5],
                                            const uint32_t n_alleles,
                                            const std::size_t ref_index,
                                            const std::size_t num_nodes,
                                            const std::size_t library_start) {

// Construct numeric genotypes
    int numeric_genotype[10][2] = {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0},
                                   {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}
    };
    for(int i = 0; i < 10; ++i) {
        int n1 = acgt_to_refalt_allele[dng::folded_diploid_nucleotides[i][0]];
        int n2 = acgt_to_refalt_allele[dng::folded_diploid_nucleotides[i][1]];
        if(n1 > n2) {
            numeric_genotype[i][0] = hts::bcf::encode_allele_unphased(n2);
            numeric_genotype[i][1] = hts::bcf::encode_allele_unphased(n1);
        } else {
            numeric_genotype[i][0] = hts::bcf::encode_allele_unphased(n1);
            numeric_genotype[i][1] = hts::bcf::encode_allele_unphased(n2);
        }
    }
    // Link VCF genotypes to our order
    int genotype_index[15];
    for(int i = 0, k = 0; i < n_alleles; ++i) {
        int n1 = refalt_to_acgt_allele[i];
        for(int j = 0; j <= i; ++j, ++k) {
            int n2 = refalt_to_acgt_allele[j];
            genotype_index[k] = (j == 0 && ref_index == 4) ?
                                -1 : dng::folded_diploid_genotypes_matrix[n1][n2];
        }
    }

    // Calculate sample genotypes
//    std::vector<int32_t> best_genotypes(2 * num_nodes);
//    std::vector<int32_t> genotype_qualities(num_nodes);
    int gt_count = n_alleles * (n_alleles + 1) / 2;
//    std::vector<float> gp_scores(gt_count*num_nodes );
    best_genotypes_.resize(2 * num_nodes);
    genotype_qualities_.resize(num_nodes);
    gp_scores_.resize( gt_count * num_nodes);

    for(size_t i = 0, k = 0; i < num_nodes; ++i) {
        size_t pos;
        double d = posterior_probabilities_[i].maxCoeff(&pos);
        best_genotypes_[2 * i] = numeric_genotype[pos][0];
        best_genotypes_[2 * i + 1] = numeric_genotype[pos][1];
        genotype_qualities_[i] = dng::util::lphred<int32_t>(1.0 - d, 255);
        // If either of the alleles is missing set quality to 0
        if(hts::bcf::allele_is_missing({best_genotypes_[2 * i]}) ||
           hts::bcf::allele_is_missing({best_genotypes_[2 * i + 1]})) {
            genotype_qualities_[i] = 0;
        }
        for(int j = 0; j < gt_count; ++j) {
            int n = genotype_index[j];
            gp_scores_[k++] = (n == -1) ? 0.0 : posterior_probabilities_[i][n];
        }
    }

    // Sample Likelihoods
//    std::vector<float> gl_scores.(gt_count *num_nodes, hts::bcf::float_missing);
    gl_scores_.resize(gt_count *num_nodes, hts::bcf::float_missing);
    for(size_t i = library_start, k = library_start * gt_count; i < num_nodes;
        ++i) {
        for(int j = 0; j < gt_count; ++j) {
            int n = genotype_index[j];
            gl_scores_[k++] = (n == -1) ? hts::bcf::float_missing :
                              genotype_likelihoods_[i][n];
        }
    }


}



void MutationStats::RecordBasicStats(hts::bcf::Variant &record){

    record.info("MUP", mup_);
    record.info("LLD", lld_);
    record.info("LLH", llh_);
    record.info("MUX", mux_); //OutputLevel::Complete
    record.info("MU1P", mu1p_); //OutputLevel::single


};

void MutationStats::RecordSingleMutationStats(hts::bcf::Variant &record){

    if(has_single_mut_) {
        record.info("DNT", dnt_);
        record.info("DNL", dnl_);
        record.info("DNQ", dnq_);
//        record.info("DNC", dnc_);
        record.samples("MU1P", node_mu1p_);
    }


};

void MutationStats::RecordGenotypeStats(hts::bcf::Variant &record){
    record.sample_genotypes(best_genotypes_);
    record.samples("GQ", genotype_qualities_);
    record.samples("GP", gp_scores_);
    record.samples("GL", gl_scores_);

}

void MutationStats::SetNodeMup(const std::vector<double> &event,
                               std::size_t first_nonfounder_index) {
    SetNodeCore(node_mup_, event, first_nonfounder_index);
}

void MutationStats::SetNodeMu1p(std::vector<double> &event, double total,
                                std::size_t first_nonfounder_index) {
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        event[i] = event[i] / total;
    }
    SetNodeCore(node_mu1p_, event, first_nonfounder_index);
}

void MutationStats::SetNodeCore(std::vector<float> &stats,
                                const std::vector<double> &event,
                                std::size_t first_nonfounder_index) {
    stats.resize(event.size(), hts::bcf::float_missing);
    for (std::size_t i = first_nonfounder_index; i < event.size(); ++i) {
        stats[i] = static_cast<float>(event[i]);
    }
    //HOWTO: use std::copy with cast??    std::copy(event.begin()+first_nonfounder_index, event.end(), stats.begin() );
}

