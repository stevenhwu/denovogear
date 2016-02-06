//
// Created by steven on 2/3/16.
//



#include "find_mutation.h"
//TODO: fix include to <*> later

//using namespace dng::task;
//using namespace dng;

// Build a list of all of the possible contigs to add to the vcf header
std::vector<std::pair<std::string, uint32_t>> parse_contigs(const bam_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    std::vector<std::pair<std::string, uint32_t>> contigs;
    uint32_t n_targets = hdr->n_targets;
    for(size_t a = 0; a < n_targets; a++) {
        if(hdr->target_name[a] == nullptr) {
            continue;
        }
        contigs.emplace_back(hdr->target_name[a], hdr->target_len[a]);
    }
    return contigs;
}


// VCF header lacks a function to get sequence lengths
// So we will extract the contig lines from the input header
std::vector<std::string> extract_contigs(const bcf_hdr_t *hdr) {
    if(hdr == nullptr)
        return {};
    // Read text of header
    int len;
    std::unique_ptr<char[], void(*)(void *)> str{bcf_hdr_fmt_text(hdr, 0, &len), free};
    if(!str)
        return {};
    std::vector<std::string> contigs;

    // parse ##contig lines
    const char *text = str.get();
    if(strncmp(text, "##contig=", 9) != 0) {
        text = strstr(text, "\n##contig=");
    } else {
        text = text - 1;
    }
    const char *end;
    for(; text != nullptr; text = strstr(end, "\n##contig=")) {
        for(end = text + 10; *end != '\n' && *end != '\0'; ++end)
            /*noop*/;
        if(*end != '\n') {
            return contigs;    // bad header, return what we have.
        }
        contigs.emplace_back(text + 1, end);
    }

    return contigs;
}



FindMutations::FindMutations(double min_prob, const Pedigree &pedigree,
                             params_t params) :
        pedigree_{pedigree}, min_prob_{min_prob},
        params_(params), genotype_likelihood_{params.params_a, params.params_b},
        work_(pedigree.CreateWorkspace()) {

    using namespace dng;
    std::cout << "theta: " << params_.theta << std::endl;

    for(auto s: params_.nuc_freq)
        std::cout << s << ' ';
    std::cout << "<- fresq" << std::endl;

    // Use a parent-independent mutation model, which produces a
    // beta-binomial
    genotype_prior_[0] = population_prior(params_.theta, params_.nuc_freq, {params_.ref_weight, 0, 0, 0});
    genotype_prior_[1] = population_prior(params_.theta, params_.nuc_freq, {0, params_.ref_weight, 0, 0});
    genotype_prior_[2] = population_prior(params_.theta, params_.nuc_freq, {0, 0, params_.ref_weight, 0});
    genotype_prior_[3] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, params_.ref_weight});
    genotype_prior_[4] = population_prior(params_.theta, params_.nuc_freq, {0, 0, 0, 0});

    // Calculate mutation expectation matrices
    full_transition_matrices_.assign(work_.num_nodes, {});
    nomut_transition_matrices_.assign(work_.num_nodes, {});
    posmut_transition_matrices_.assign(work_.num_nodes, {});
    onemut_transition_matrices_.assign(work_.num_nodes, {});
    mean_matrices_.assign(work_.num_nodes, {});

    for(size_t child = 0; child < work_.num_nodes; ++child) {


        auto trans = pedigree.transitions()[child];
        std::cout << "Node:" << child << "\t" << "\t" << trans.length1 << "\t" << trans.length2 << "\t";
//                                       << std::endl;
        if(trans.type == Pedigree::TransitionType::Germline) {
            std::cout << "==Germline.\n";
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
            std::cout << "==Somatic/Library.\n";
            auto orig = f81::matrix(trans.length1, params_.nuc_freq);

            full_transition_matrices_[child] = mitosis_diploid_matrix(orig);
            nomut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 0);
            posmut_transition_matrices_[child] = full_transition_matrices_[child] -
                                                 nomut_transition_matrices_[child];
            onemut_transition_matrices_[child] = mitosis_diploid_matrix(orig, 1);
            mean_matrices_[child] = mitosis_diploid_mean_matrix(orig);
        } else {
            std::cout << "==Other.\n";
            full_transition_matrices_[child] = {};
            nomut_transition_matrices_[child] = {};
            posmut_transition_matrices_[child] = {};
            onemut_transition_matrices_[child] = {};
            mean_matrices_[child] = {};
        }

    }
    std::cout << "" << std::endl;
    //Calculate max_entropy based on having no data
    for(int ref_index = 0; ref_index < 5; ++ref_index) {
        work_.SetFounders(genotype_prior_[ref_index]);

        pedigree_.PeelForwards(work_, nomut_transition_matrices_);
        pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
        event_.assign(work_.num_nodes, 0.0);
        double total = 0.0, entropy = 0.0;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                                   work_.lower[i].matrix().transpose()).array() *
                                  onemut_transition_matrices_[i].array();
            total += mat.sum();
            entropy += (mat.array() == 0.0).select(mat.array(),
                                                   mat.array() * mat.log()).sum();
        }
        // Calculate entropy of mutation location
        max_entropies_[ref_index] = (-entropy / total + log(total)) / M_LN2;
        std::cout << "Ref: " << ref_index << "\t" << max_entropies_[ref_index] << std::endl;
    }


    std::cout << "END FM const" << std::endl;


}

// Returns true if a mutation was found and the record was modified
bool FindMutations::operator()(const std::vector<depth_t> &depths,
                               int ref_index, stats_t *stats) {
    using namespace std;
    using namespace hts::bcf;
    using dng::util::lphred;
    using dng::util::phred;


    MutationStats mutation_stats(min_prob_);
    assert(stats != nullptr);

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
                genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }


    std::cout << "FM() depth: " << "\t" << depths.size() << std::endl;
    for (auto d : depths) {
        auto cc = d.counts;
        std::cout << cc[0] << " "<< cc[1] << " "<< cc[2] << " "<< cc[3] << " "<< std::endl;
    }

    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    // Calculate log P(Data, nomut ; model)
    const double logdata_nomut = pedigree_.PeelForwards(work_, nomut_transition_matrices_);

//    std::cout << "exp lower\n" << work_.lower[0] << std::endl;
    /**** Forward-Backwards with full-mutation ****/

    // Calculate log P(Data ; model)
    const double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    const double pmut = -std::expm1(logdata_nomut - logdata) ;//+ 10000;
    std::cout << logdata_nomut << "\t" << logdata << "\t" << pmut << "\t" << min_prob_ << std::endl;
    // Skip this site if it does not meet lower probability threshold

    stats->mup = pmut;
    stats->lld = (logdata + scale) / M_LN10;
    stats->llh = logdata / M_LN10;
    mutation_stats.store_mup(logdata_nomut, logdata);

//    if(mutation_stats.check_mutation_prob_lt_threshold() ){
//        return false;
//    }
//    if(pmut < min_prob_) {
//        return false;
//    }


    // Copy genotype likelihoods
    std::cout << "Copy Genotype likelihood" << std::endl;
    mutation_stats.store_genotype_likelihood(depths, work_);
    stats->genotype_likelihoods.resize(work_.num_nodes);
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::size_t pos = work_.library_nodes.first + u;
        stats->genotype_likelihoods[pos] = work_.lower[pos].log() / M_LN10;
    }

    mutation_stats.store_scaled_log_likelihood(scale);
    size_t library_start = work_.library_nodes.first; //XXX: Can delete?
    stats->mup = pmut;
    stats->lld = (logdata + scale) / M_LN10;
    stats->llh = logdata / M_LN10;

    mutation_stats.store_posterior_probabilities(work_);
    // Peel Backwards with full-mutation
    std::cout << "Peel Backwards with full-mutation" << std::endl;
    //TODO: After refactor, how to prevent calling something multiple times??
    pedigree_.PeelBackwards(work_, full_transition_matrices_);
//    pedigree_.PeelForwards(work_, nomut_transition_matrices_);
//    pedigree_.PeelBackwards(work_, nomut_transition_matrices_);

    // Calculate statistics after Forward-Backwards
    stats->posterior_probabilities.resize(work_.num_nodes);
    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
//        std::cout << stats->posterior_probabilities[i].prod() << "\t";
    }
//    std::cout << std::endl;




    double mux = 0.0;
    event_.assign(work_.num_nodes, 0.0);
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        mux += (work_.super[i] * (mean_matrices_[i] * work_.lower[i].matrix()).array()).sum();
        event_[i] =  (work_.super[i] * (posmut_transition_matrices_[i] * work_.lower[i].matrix()).array()).sum();
        event_[i] = event_[i] / pmut;
    }
    stats->mux = mux;

    stats->node_mup.resize(work_.num_nodes, hts::bcf::float_missing);
    for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        stats->node_mup[i] = static_cast<float>(event_[i]);
    }

    /**** Forward-Backwards with no-mutation ****/
    std::cout << "/**** Forward-Backwards with no-mutation ****/" << std::endl;
    // TODO: Better to use a separate workspace???
    pedigree_.PeelForwards(work_, nomut_transition_matrices_);
    pedigree_.PeelBackwards(work_, nomut_transition_matrices_);
    event_.assign(work_.num_nodes, 0.0);
    double total = 0.0, entropy = 0.0, max_coeff = -1.0;
    size_t dn_loc = 0, dn_col = 0, dn_row = 0;
    for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
        Eigen::ArrayXXd mat = (work_.super[i].matrix() *
                               work_.lower[i].matrix().transpose()).array() *
                              onemut_transition_matrices_[i].array();
        std::size_t row, col;
        double mat_max = mat.maxCoeff(&row, &col);
        if(mat_max > max_coeff) {
            max_coeff = mat_max;
            dn_row  = row;
            dn_col = col;
            dn_loc = i;
        }
        event_[i] = mat.sum();
        entropy += (mat.array() == 0.0).select(mat.array(),
                                               mat.array() * mat.log()).sum();
        total += event_[i];
    }
    // Calculate P(only 1 mutation)
    const double pmut1 = total * (1.0 - pmut);
    stats->mu1p = pmut1;

    std::cout << pmut1 << "\t" << pmut << std::endl;

    // Output statistics for single mutation only if it is likely
    if(pmut1 / pmut >= min_prob_) {
        stats->has_single_mut = true;
        for(std::size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            event_[i] = event_[i] / total;
        }

        // Calculate entropy of mutation location
        entropy = (-entropy / total + log(total)) / M_LN2;
        entropy /= max_entropies_[ref_index];
        stats->dnc = std::round(100.0 * (1.0 - entropy));

        stats->dnq = lphred<int32_t>(1.0 - (max_coeff / total), 255);
        stats->dnl = pedigree_.labels()[dn_loc];
        if(pedigree_.transitions()[dn_loc].type == Pedigree::TransitionType::Germline) {
            stats->dnt = &meiotic_diploid_mutation_labels[dn_row][dn_col][0];
        } else {
            stats->dnt = &mitotic_diploid_mutation_labels[dn_row][dn_col][0];
        }

        stats->node_mu1p.resize(work_.num_nodes, hts::bcf::float_missing);
        for(size_t i = work_.founder_nodes.second; i < work_.num_nodes; ++i) {
            stats->node_mu1p[i] = static_cast<float>(event_[i]);
        }
    } else {
        stats->has_single_mut = false;
    }
    return true;
}


// Returns true if a mutation was found and the record was modified
bool FindMutations::calculate_mutation(const std::vector<depth_t> &depths,
                               int ref_index, MutationStats &mutation_stats) {

    using namespace std;
    using namespace hts::bcf;
    using dng::util::lphred;
    using dng::util::phred;


//    assert(stats != nullptr); //TODO: Implement similar check or reset for MutationStats

    // calculate genotype likelihoods and store in the lower library vector
    double scale = 0.0, stemp;
    for(std::size_t u = 0; u < depths.size(); ++u) {
        std::tie(work_.lower[work_.library_nodes.first + u], stemp) =
                genotype_likelihood_(depths[u], ref_index);
        scale += stemp;
    }
    // Set the prior probability of the founders given the reference
    work_.SetFounders(genotype_prior_[ref_index]);

    calculate_mup(mutation_stats);

    if(mutation_stats.check_mutation_prob_lt_threshold() ){
        return false;
    }

    mutation_stats.store_scaled_log_likelihood(scale);
    mutation_stats.store_genotype_likelihood(depths, work_);


    mutation_stats.store_posterior_probabilities(work_);
    // Peel Backwards with full-mutation
    std::cout << "Peel Backwards with full-mutation" << std::endl;
    pedigree_.PeelBackwards(work_, full_transition_matrices_);

    // Calculate statistics after Forward-Backwards

//    stats->posterior_probabilities.resize(work_.num_nodes);
//    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
//        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
//        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
//    }



    // Peel Backwards with full-mutation
    std::cout << "Peel Backwards with full-mutation" << std::endl;
    pedigree_.PeelBackwards(work_, full_transition_matrices_);
    mutation_stats.store_posterior_probabilities(work_);

    size_t library_start = work_.library_nodes.first;

//    // Calculate statistics after Forward-Backwards
//    stats->posterior_probabilities.resize(work_.num_nodes);
//    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
//        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
//        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
//    }
//    double mux = 0.0;

    calculate_posterior_probabilities(mutation_stats);

//    // Peel Backwards with full-mutation
//    std::cout << "Peel Backwards with full-mutation" << std::endl;
//    pedigree_.PeelBackwards(work_, full_transition_matrices_);
//
//    // Calculate statistics after Forward-Backwards
//
//    stats->posterior_probabilities.resize(work_.num_nodes);
//    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
//        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
//        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
//    }
//


    return true;

}

void FindMutations::calculate_mup(MutationStats &mutation_stats) {
    // Calculate log P(Data, nomut ; model)
    const double logdata_nomut = pedigree_.PeelForwards(work_, nomut_transition_matrices_);

    /**** Forward-Backwards with full-mutation ****/
    // Calculate log P(Data ; model)
    const double logdata = pedigree_.PeelForwards(work_, full_transition_matrices_);

    // P(mutation | Data ; model) = 1 - [ P(Data, nomut ; model) / P(Data ; model) ]
    mutation_stats.store_mup(logdata_nomut, logdata);
    std::cout << logdata_nomut << "\t" << logdata << "\t" << mutation_stats.get_mup() << "\t" << min_prob_ << std::endl;

}


void FindMutations::calculate_posterior_probabilities(MutationStats &mutation_stats) {

    pedigree_.PeelBackwards(work_, full_transition_matrices_);

    mutation_stats.store_posterior_probabilities(work_);
//    stats->posterior_probabilities.resize(work_.num_nodes);
//    for(std::size_t i = 0; i < work_.num_nodes; ++i) {
//        stats->posterior_probabilities[i] = work_.upper[i] * work_.lower[i];
//        stats->posterior_probabilities[i] /= stats->posterior_probabilities[i].sum();
//    }
}