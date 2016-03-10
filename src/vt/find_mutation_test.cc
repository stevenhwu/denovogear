//
// Created by steven on 1/11/16.
//
//#define BOOST_TEST_MODULE "VT_TEST"

#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iosfwd>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <sstream>
#include <string>

#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm/replace.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/algorithm/string.hpp>

#include <dng/task/call.h>
#include <dng/pedigree.h>
#include <dng/fileio.h>
#include <dng/pileup.h>
#include <dng/read_group.h>
#include <dng/likelihood.h>
#include <dng/seq.h>
#include <dng/utility.h>
#include <dng/hts/bcf.h>
#include <dng/hts/extra.h>
#include <dng/vcfpileup.h>
#include <dng/mutation.h>
#include <dng/stats.h>

#include <htslib/faidx.h>
#include <htslib/khash.h>
#include <dng/app.h>

#include "version.h"



#include <dng/relationship_graph.h>

#include <utils/vcf_utils.h>
#include <utils/assert_utils.h>
#include <utils/boost_utils.h>
#include <dng/find_mutation.h>
#include <include/dng/workspace.h>


//#include <boost/test/unit_test.hpp>

#include <dng/find_mutation_x.h>

using namespace dng::task;
using namespace dng;


// The main loop for dng-call application
// argument_type arg holds the processed command line arguments
int dng::task::Call::operator()(dng::task::Call::argument_type &arg) {
    using namespace std;
    using namespace hts::bcf;
    using dng::utility::lphred;
    using dng::utility::phred;

    std::cout << arg.mu << "\t" << arg.mu_somatic << "\t" << arg.mu_library << std::endl;

    // Parse pedigree from file
    dng::io::Pedigree ped;

    if (!arg.ped.empty()) {
        ifstream ped_file(arg.ped);
        if (!ped_file.is_open()) {
            throw std::runtime_error(
                    "unable to open pedigree file '" + arg.ped + "'.");
        }
        ped.Parse(utils::istreambuf_range(ped_file));
    } else {
        throw std::runtime_error("pedigree file was not specified.");
    }

    // Open Reference
    unique_ptr<char[], void (*)(void *)> ref{nullptr, free};
    int ref_sz = 0, ref_target_id = -1;
    unique_ptr<faidx_t, void (*)(faidx_t *)> fai{nullptr, fai_destroy};
    if (!arg.fasta.empty()) {
        fai.reset(fai_load(arg.fasta.c_str()));
        if (!fai)
            throw std::runtime_error("unable to open faidx-indexed reference file '"
                                     + arg.fasta + "'.");
    }

    // Parse Nucleotide Frequencies
    std::array<double, 4> freqs;
    {
        auto f = utility::parse_double_list(arg.nuc_freqs, ',', 4);
        std::copy(f.first.begin(), f.first.end(), &freqs[0]);
    }

    // quality thresholds
    int min_qual = arg.min_basequal;
    double min_prob = arg.min_prob;

    // Open input files
    dng::ReadGroups rgs;
    vector<hts::File> indata;
    vector<hts::bam::File> bamdata;
    vector<hts::bcf::File> bcfdata;
    for (auto &&str : arg.input) {
        indata.emplace_back(str.c_str(), "r");
        if (indata.back().is_open()) {
            continue;
        }
        throw std::runtime_error("unable to open input file '" + str + "'.");
    }

    // Check to see if all inputs are of the same type
    const htsFormatCategory cat = indata[0].format().category;
    for (auto &&f : indata) {
        if (f.format().category == cat) {
            continue;
        }
        throw std::runtime_error("mixing sequence data and variant data as input is not supported.");
    }

    // Begin writing VCF header
    auto out_file = vcf_get_output_mode(arg);
    hts::bcf::File vcfout(out_file.first.c_str(), out_file.second.c_str());
    vcf_add_header_text(vcfout, arg);

    if (cat == sequence_data) {
        // Wrap input in hts::bam::File
        for (auto &&f : indata) {
            bamdata.emplace_back(std::move(f), arg.region.c_str(), arg.fasta.c_str(),
                                 arg.min_mapqual, arg.header.c_str());
        }

        // Read header from first file
        const bam_hdr_t *h = bamdata[0].header();

        // Add contigs to header
        for (auto &&contig : parse_contigs(h)) {
            vcfout.AddContig(contig.first.c_str(), contig.second);
        }

        // Add each genotype/sample column
        rgs.ParseHeaderText(bamdata, arg.rgtag);
    } else if (cat == variant_data) {
        cout << "CALL VCF DATA\n" << endl;
        bcfdata.emplace_back(std::move(indata[0]));
        // Read header from first file
        const bcf_hdr_t *h = bcfdata[0].header();



        // Add contigs to header
        for (auto &&contig : extract_contigs(h)) {
            vcfout.AddHeaderMetadata(contig.c_str());
//            cout << contig.c_str() << endl;
        }
        // Add each genotype/sample column
        rgs.ParseSamples(bcfdata[0]);
    } else {
        throw runtime_error("unsupported file category.");
    }

    // Construct peeling algorithm from parameters and pedigree information
    dng::Pedigree pedigree;
    cout << "MU: " << arg.mu << "\t" << arg.mu_somatic << "\t" << arg.mu_library << std::endl;
    if (!pedigree.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library)) {
        throw std::runtime_error("Unable to construct peeler for pedigree; "
                                         "possible non-zero-loop pedigree.");
    }
    if (arg.gamma.size() < 2) {
        throw std::runtime_error("Unable to construct genotype-likelihood model; "
                                         "Gamma needs to be specified at least twice to change model from default.");
    }

    for (auto &&line : pedigree.BCFHeaderLines()) {
        vcfout.AddHeaderMetadata(line.c_str());
//        cout << line.c_str() << endl;
    }
//    pedigree.PrintTable(std::cout);
//    pedigree.PrintMachine(std::cout);
//    for(auto && line : pedigree.BCFHeaderLines()) {
//        cout << line.c_str() << endl;
//    }
//    pedigree.PrintStates(std::cout);

    dng::PedigreeV2 pedigree2;
    pedigree2.Construct(ped, rgs, arg.mu, arg.mu_somatic, arg.mu_library);
//    auto compare = pedigree.Equal(pedigree2);
//    std::cout << "Pedigree == PedigreeV2? " <<  std::boolalpha << compare << std::endl;



    /*{
        std::cout << "Founder, nonfounder, somatic, library, num_nodes" << std::endl;
        cout << pedigree.first_founder_ << "\t" << //endl;    // start of founder germline
        pedigree.first_nonfounder_ << "\t" <<  //endl; // start of non-founder germline
        pedigree.first_somatic_ << "\t" <<  //endl;    // start of somatic nodes
        pedigree.first_library_ << "\t" <<  //endl;    // start of libraries
        pedigree.num_nodes_ << endl;        // total number of nodes

        cout << "Roots\t" << pedigree.roots_.size() << endl;
        // Pedigree Structure
        cout << "Labels: " << pedigree.labels_.size() << endl;
        cout << pedigree.transitions_.size() << endl;


        // The original, simplified peeling operations
        cout << pedigree.peeling_ops_.size() << endl;
        for (int j = 0; j < pedigree.peeling_ops_.size(); ++j) {
            int op_index = pedigree.peeling_ops_[j];
            std::cout << "op: " << pedigree.peeling_ops_[j] << "\t" <<
            pedigree.peeling_functions_ops_[j] << "\t" <<
            peel::op::map_enum_string[pedigree.peeling_functions_ops_[j]] << std::endl;
        }

        // The arguments to a peeling operation
        cout << pedigree.family_members_.size() << endl;
    }*/

//    FindMutations calculate{min_prob, pedigree,
//                            {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]}};
    int ref_index = 2;
    std::vector<depth_t> read_depths{3};
    uint16_t cc[3][4] = {{0, 1, 25, 29},
                         {0, 0, 57, 0},
                         {0, 0, 76, 1}};
    for (int j = 0; j < 3; ++j) {
        std::copy(cc[j], cc[j] + 4, read_depths[j].counts);
    }


    FindMutations::FindMutationParams test_param_1
            {arg.theta, freqs, arg.ref_weight, arg.gamma[0], arg.gamma[1]};
    // Pileup data
    FindMutations find_mutation{min_prob, pedigree, test_param_1};
//    test_basic_parameterts(find_mutation);
//    test_prior(find_mutation);
////
//    test_full_transition(find_mutation);
////
//    test_operator(find_mutation);

    workspace_t tt;// = pedigree.CreateWorkspace();
    tt.Resize(10);
    tt.Cleanup();



    for(std::size_t u = 0; u < 10; ++u) {
        tt.lower[tt.library_nodes.first + u] = GenotypeArray::Random();
    }
    workspace_t t2;
    t2 = tt;

    std::cout << tt.lower[0]<< std::endl;
    std::cout << "\n"<< std::endl;
    std::cout << t2.lower[0]<< std::endl;

    std::cout << "\n"<< std::endl;
    tt.lower[0] = GenotypeArray::Random();
    std::cout << tt.lower[0]<< std::endl;
    std::cout << "\n"<< std::endl;
    std::cout << t2.lower[0]<< std::endl;

    FindMutations fm_x{0.01, pedigree, {arg.theta, freqs, arg.ref_weight,
                                        arg.gamma[0], arg.gamma[1]}};
//    FindMutationsXLinked fm_x{pedigree, find_mutation_params};
//
//    fm_x.init();
//    MutationStats stats = {0.001};
//    fm_x.run(read_depths, ref_index, stats);
//
    std::cout << "ready to exit: " << EXIT_SUCCESS << std::endl;
    return EXIT_SUCCESS;

}





int main(int argc, char *argv[]) {


//BOOST_AUTO_TEST_CASE(Test_FM){


    try {
//        return CallApp(argc, argv)();
//        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;

//        char *argv[] = {"test", "-p", "arg2", NULL};
//        int argc = sizeof(argv) / sizeof(char*) - 1;
        int argc=4;
        char *argv[argc+1];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";
        argv[2] = (char*) "testdata/sample_5_3/ceu.ped"; //"pedFile";
        argv[3] = (char*) "testdata/sample_5_3/test1.vcf"; //test1.bam
        argv[2] = (char*) "testDataSW/ceu_M12.ped"; //"pedFile";
        argv[3] = (char*) "testDataSW/test_M12.vcf"; //test1.bam
        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
        a();

    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
    return EXIT_FAILURE;

}

//
// get unknow error message
//unknown location(0): fatal error in "Test_FM": signal: illegal operand; address of failing instruction: 0x0061ef82
//
//*** 1 failure detected in test suite "VT_TEST"

//BOOST_AUTO_TEST_CASE(Test_FM){
//    try {
////        return CallApp(argc, argv)();
////        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
//
////        char *argv[] = {"test", "-p", "arg2", NULL};
////        int argc = sizeof(argv) / sizeof(char*) - 1;
//        int argc=4;
//        char *argv[argc+1];
//        char const *p = "abc";
//        argv[0] = (char*) "test";
//        argv[1] = (char*) "-p";
//        argv[2] = "testdata/sample_5_3/ceu.ped"; //"pedFile";
//        argv[3] = "testdata/sample_5_3/test1.vcf"; //test1.bam
//        dng::CommandLineApp<dng::task::Call> a (argc, argv) ;
//        std::cout << "DIV" << std::endl;
//        BOOST_CHECK_EQUAL(0.01, 0.011);
//        a();
//    } catch(std::exception &e) {
//        std::cerr << e.what() << std::endl;
//    }
////    return EXIT_FAILURE;
//
//}
//
