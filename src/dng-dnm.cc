#include "denovogear.h"
#include "pedParser.h"

#include <iostream>
#include <string>
#include <vector>

#include <dng/hts/extra.h>
#include <dng/app.h>
#include <dng/task/dnm.h>
#include "htslib/synced_bcf_reader.h"

using namespace std;
using namespace dng::task;

int RD_cutoff = 10; // cutoff for the read depth filter
double PP_cutoff = 0.0001; // posterior probability cutoff

double indel_mrate =
    0; // indel mutation prior is based on a log linear model unless specified by the user.
double snp_mrate = 1e-8; // snp mutation prior
double poly_rate = 1e-3; // polymorphism rate - used in prior calculations
double pair_mrate = 1e-9; // mutation prior for paired samples
double mu_scale = 1.0; // scaling factor for indel priors



int callMakeSNPLookup(lookup_table_t &tgtSNP, lookup_snp_t &lookupSNP,
                      std::string model) {
    std::vector<string> tmp;
    for(int l = 0; l < 10; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 100; l++) {
        tgtSNP.push_back(tmp);
    }


    if(model == "auto") { makeSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else if(model == "XS") { makeXSSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else if(model == "XD") { makeXDSNPLookup(snp_mrate, poly_rate, tgtSNP, lookupSNP); }
    else {
        cerr << endl << "Invalid model for SNP lookup, exiting.";
        exit(1);
    }

    cerr << "\nCreated SNP lookup table\n";
    cerr << " First mrate: " << lookupSNP.mrate(0,
            0) << " last: " << lookupSNP.mrate(99, 9) << endl;
    cerr << " First code: " << lookupSNP.code(0,
            0) << " last: " << lookupSNP.code(99, 9) << endl;
    cerr << " First target string: " << tgtSNP[0][0] << " last: " << tgtSNP[99][9]
         << endl;
    cerr << " First tref: " << lookupSNP.tref(0,
            0) << " last: " << lookupSNP.tref(99, 9) << endl;
    return 0;
}

int callMakeINDELLookup(lookup_table_t &tgtIndel, lookup_indel_t &lookupIndel,
                        std::string model) {
    std::vector<std::string> tmp;
    for(int l = 0; l < 3; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 9; l++) {
        tgtIndel.push_back(tmp);
    }


    if(model == "auto") { makeIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else if(model == "XS") { makeXSIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else if(model == "XD") { makeXDIndelLookup(poly_rate, tgtIndel, lookupIndel); }
    else {
        std::cerr << std::endl << "Invalid model for INDEL lookup, exiting.";
        exit(1);
    }

    std::cerr << "\nCreated indel lookup table";
    std::cerr << " First code: " << lookupIndel.code(0,
              0) << " last: " << lookupIndel.code(8, 2) << std::endl;
    std::cerr << " First target string: " << tgtIndel[0][0] << " last: " <<
              tgtIndel[8][2] << std::endl;
    std::cerr << " First prior: " << lookupIndel.priors(0,
              0) << " last: " << lookupIndel.priors(8, 2) << std::endl;
    return 0;
}

int callMakePairedLookup(lookup_table_t &tgtPair, lookup_pair_t &lookupPair) {
    std::vector<std::string> tmp;
    for(int l = 0; l < 10; l++) {
        tmp.push_back("NA");
    }
    for(int l = 0; l < 10; l++) {
        tgtPair.push_back(tmp);
    }
    makePairedLookup(pair_mrate, tgtPair, lookupPair);
    std::cerr << "\nCreated paired lookup table" << std::endl;
    std::cerr << " First target string: " << tgtPair[0][0] << " last: " <<
              tgtPair[9][9] << std::endl;
    std::cerr << " First prior " << lookupPair.priors(0,
              0) << " last: " << lookupPair.priors(9, 9) << std::endl;
    return 0;
}

void writeVCFHeader(hts::bcf::File &vcfout, std::string &bcf_file,
                    std::string &ped_file, std::vector<std::string> &samples) {

    vcfout.AddHeaderMetadata("##fileformat=VCFv4.1");
    vcfout.AddHeaderMetadata("##source=", PACKAGE_STRING);
    vcfout.AddHeaderMetadata("##input_bcf=", bcf_file);
    vcfout.AddHeaderMetadata("##input_ped=", ped_file);
    vcfout.AddHeaderMetadata("##cutoff_read_depth=", RD_cutoff);
    vcfout.AddHeaderMetadata("##cutoff_posterior_probability=", PP_cutoff);
    vcfout.AddHeaderMetadata("##mutation_rate_snp=", snp_mrate);
    vcfout.AddHeaderMetadata("##mutation_rate_indel=", indel_mrate);
    vcfout.AddHeaderMetadata("##mutation_rate_paired=", pair_mrate);
    vcfout.AddHeaderMetadata("##mutation_rate_polymorphism=", poly_rate);
    vcfout.AddHeaderMetadata("##mutation_rate_indel_scaling_constant=", mu_scale);

#ifndef NEWVCFOUT
    // Old format, doesn't work well with htslib
    vcfout.AddHeaderMetadata("##INFO=<ID=RD_MOM,Number=1,Type=Integer,Description=\"Read depth for MOM\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=RD_DAD,Number=1,Type=Integer,Description=\"Read depth for DAD\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=RD_NORMAL,Number=1,Type=Integer,Description=\"Read depth for Normal sample in Paired Sample Analysis\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQ_MOM,Number=1,Type=Integer,Description=\"Mapping quality for MOM\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQ_DAD,Number=1,Type=Integer,Description=\"Mapping quality for DAD\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=MQ_NORMAL,Number=1,Type=Integer,Description=\"Mapping quality for Normal sample in Paired Sample Analysis\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=NULL_CONFIG(child/mom/dad),Number=1,Type=String,Description=\"NULL trio configuration\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=PP_NULL,Number=1,Type=Float,Description=\"Posterior probability for the NULL configuration\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=ML_NULL,Number=1,Type=Float,Description=\"Maximum Likelihood for the NULL configuration\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=ML_DNM,Number=1,Type=Float,Description=\"Maximum Likelihood for the DNM configuration\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=SNPcode,Number=1,Type=Integer,Description=\"Code that shows whether DNM configuration shown is monomorphic or contains variation\">");
    vcfout.AddHeaderMetadata("##INFO=<ID=INDELcode,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
    vcfout.AddHeaderMetadata("##INFO=<ID=PairSNPcode,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
    vcfout.AddHeaderMetadata("##INFO=<ID=code,Number=1,Type=Integer,Description=\" \">"); //// Needs Detail
    vcfout.AddHeaderMetadata("##FORMAT=<ID=DNM_CONFIG(child/mom/dad),Number=1,Type=String,Description=\"DNM trio configuration\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=PP_DNM,Number=1,Type=Float,Description=\"Posterior probability for the DNM configuration\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read Depth of the child\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality of the child\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=RD_T,Number=1,Type=Integer,Description=\"Read Depth of the tumor sample\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ_T,Number=1,Type=Integer,Description=\"Mapping quality of the normal sample\">");
    vcfout.AddSample(samples[0].c_str());
#else
    // Newer header, produces more standarized VCF output
    vcfout.AddHeaderMetadata("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Read Depth\">");
    vcfout.AddHeaderMetadata("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mapping quality\">");
    for(std::string sample : samples) {
        vcfout.AddSample(sample.c_str());
    }
#endif
    vcfout.WriteHeader();
}



int DNM::operator()(std::string &model, DNM::argument_type &arg) {
    if(model != "auto" && model != "XS" && model != "XD") {
        throw std::runtime_error("Invalid model option " + model +
                                 ". Use auto, XS, or XD.");
    }

    std::string input_file;
    if(!arg.bcf.empty()) {
    	if(!arg.vcf.empty()) {
    		throw std::runtime_error("ERROR ! Attempting to use both a vcf and bcf file, can only use one. Exiting !");
    	}
    	input_file = arg.bcf;
    }
    else if(!arg.vcf.empty()) {
    	input_file = arg.vcf;
    }
    else {
    	throw std::runtime_error("ERROR ! please specify input variant file with --bcf command! Exiting!");
    }

    if(arg.ped.empty()) {
        throw std::runtime_error("ERROR ! No PED file! Exiting!");
    }

    snp_mrate = arg.snp_mrate;
    indel_mrate = arg.indel_mrate;
    poly_rate = arg.poly_rate;
    pair_mrate = arg.pair_mrate;
    mu_scale = arg.mu_scale;
    PP_cutoff = arg.pp_cutoff;
    RD_cutoff = arg.rd_cutoff;

    // Create SNP lookup
    lookup_snp_t lookupSNP;
    lookup_table_t tgtSNP;
    callMakeSNPLookup(tgtSNP, lookupSNP, model);

    // Create INDEL lookup
    lookup_indel_t lookupIndel;
    lookup_table_t tgtIndel;
    callMakeINDELLookup(tgtIndel, lookupIndel, model);

    // Create paired lookup
    lookup_pair_t lookupPair;
    lookup_table_t tgtPair;
    callMakePairedLookup(tgtPair, lookupPair);

    // TODO: Use Reed's PED parser
    // Iterate each position of BCF file
    Trio *trios;
    Pair *pairs;
    int trio_count = 0, pair_count = 0;
    parse_ped(arg.ped, &trios, &pairs, trio_count, pair_count);

    //create output vcf -- assumes theres only one trio/pair
    std::vector<std::string> samples;
    if(trio_count > 0) {
        samples.push_back(trios[0].cID);
        samples.push_back(trios[0].mID);
        samples.push_back(trios[0].fID);
    } else {
        samples.push_back(pairs[0].tumorID);
        samples.push_back(pairs[0].normalID);
    }

    std::vector<hts::bcf::File> output_vcf;
    if(!arg.write.empty()) {
        output_vcf.emplace_back(arg.write.c_str(), "w");
        writeVCFHeader(output_vcf[0], input_file, arg.ped, samples);
    }

    qcall_t mom_snp, dad_snp, child_snp;
    indel_t mom_indel, dad_indel, child_indel;
    pair_t tumor, normal;
    //bcf_hdr_t *hout, *hin;
    //int tid, begin, end;

    int snp_total_count = 0, snp_pass_count = 0;
    int indel_total_count = 0, indel_pass_count = 0;
    int pair_total_count = 0, pair_pass_count = 0;


    bcf_srs_t *rec_reader = bcf_sr_init(); // used for iterating each rec in BCF/VCF

    // Open region if specified
    if(!arg.region.empty()) {
        int ret = bcf_sr_set_regions(rec_reader, arg.region.c_str(), 0);
        if(ret == -1) {
            throw std::runtime_error("no records in the query region " + arg.region);
        }
    }

    // Initialize the record reader to iterate through the BCF/VCF input
    int ret = bcf_sr_add_reader(rec_reader, input_file.c_str());
    if(ret == 0) {
        int errnum = rec_reader->errnum;
        switch(errnum) {
        case not_bgzf:
            throw std::runtime_error("Input file type does not allow for region searchs. Exiting!");
            break;
        case idx_load_failed:
            throw std::runtime_error("Unable to load query region, no index. Exiting!");
            break;
        case file_type_error:
            throw std::runtime_error("Could not load filetype. Exiting!");
            break;
        default:
            throw std::runtime_error("Could not load input sequence file into htslib. Exiting!");
        };
    }

    // Get the header (should be only one input file)
    const bcf_hdr_t *hdr = bcf_sr_get_header(rec_reader, 0);

    
    while(bcf_sr_next_line(rec_reader)) {
        bcf1_t *rec = bcf_sr_get_line(rec_reader, 0);
        int j = 0;
        int flag = 0;
        for(j = 0; j < trio_count; j++) {
            bcf_unpack(rec, BCF_UN_STR);
            int is_indel = bcf_2qcall(hdr, rec, trios[j],
                                      &mom_snp, &dad_snp, &child_snp,
                                      &mom_indel, &dad_indel, &child_indel,
                                      flag);

            if(is_indel == 0) {
                snp_total_count++;
                trio_like_snp(child_snp, mom_snp, dad_snp, flag,
                              tgtSNP, lookupSNP, output_vcf,//vcf_out,
                              PP_cutoff, RD_cutoff, snp_pass_count);
            } else if(is_indel == 1) {

                indel_total_count++;
                trio_like_indel(&child_indel, &mom_indel, &dad_indel, flag,
                                tgtIndel, lookupIndel, mu_scale, output_vcf,
                                PP_cutoff, RD_cutoff, indel_pass_count, indel_mrate);

            } else if(is_indel < 0) {
                printf("\n BCF PARSING ERROR - Trios!  %d\n Exiting !\n", is_indel);
                exit(1);
            }
        }

        // PROCESS  PAIRS
        if(model == "auto") { // paired sample model not developed for XS, XD yet
            for(j = 0; j < pair_count; j++) {
                int is_indel = bcf2Paired(hdr, rec, pairs[j], &tumor, &normal, flag);
                if(is_indel == 0) {
                    pair_total_count++;
                    pair_like(tumor, normal, tgtPair, lookupPair, flag, output_vcf,
                              PP_cutoff, RD_cutoff, pair_pass_count);
                } else if(is_indel < 0) {
                    printf("\n BCF PARSING ERROR - Paired Sample!  %d\n Exiting !\n", is_indel);
                    exit(1);
                }
            }
        }
    }

    cerr << endl << "Total number of SNP sites interrogated: " << snp_total_count;
    cerr << endl << "Total number of SNP sites passing read-depth filters: " <<
         snp_pass_count;
    cerr << endl << "Total number of INDEL sites interrogated: " <<
         indel_total_count;
    cerr << endl << "Total number of INDEL sites passing read-depth filters: " <<
         indel_pass_count;
    cerr << endl << "Total number of Paired sample sites interrogated: " <<
         pair_total_count;
    cerr << endl <<
         "Total number of Paired sample sites passing read-depth filters: " <<
         pair_pass_count;

    std::cerr << std::endl << "Done !" << std::endl;
    exit(0);
}



/****************************************************************************************
int DNM::operator()(std::string &model, DNM::argument_type &arg) {
  std::cout << "HERE with " << model << std::endl;
}
*/

typedef dng::CommandLineApp<dng::task::DNM> App;

class DNGApp : App {
public:
    std::string model;

    DNGApp(int argc, char *argv[]) : App(argc, argv) {
        if(argc > 1) {
            model = argv[1];
        }
    }

    int operator()() {
        using namespace std;
        if(arg.version) {
            return CmdVersion();
        }
        if(arg.help || arg.input.empty()) {
            return CmdHelp();
        }
        return task_(model, arg);
    }

protected:
    virtual int CmdHelp() const {
        string usage_name(arg.run_name);
        if(usage_name.substr(0, 4) == "dng-") {
            usage_name[3] = ' ';
        }
        cerr << "Usage:" << std::endl
             << "Autosomes:" << std::endl
             << "\tdng dnm auto --bcf bcf_f --ped ped_f [OR] dng dnm auto --vcf vcf_f --ped ped_f\n"
             << "X chromosome in male offspring:\n"
             << "\tdng dnm XS --bcf bcf_f --ped ped_f [OR] dng dnm XS --vcf vcf_f --ped ped_f\n"
             << "X chromosome in female offspring:\n"
             << "\tdng dnm XD --bcf bcf_f --ped ped_f [OR] dng dnm XD --vcf vcf_f --ped ped_f\n"
             << endl;
        cerr << ext_desc_ << endl;
        return EXIT_SUCCESS;
    }

};

int main(int argc, char *argv[]) {
    try {
        return DNGApp(argc, argv)();
        //return CallApp(argc, argv)();
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }

    return EXIT_FAILURE;
}
