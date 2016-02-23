/*
 * Copyright (c) 2010, 2011 Genome Research Ltd.
 * Copyright (c) 2012, 2013 Donald Conrad and Washington University in St. Louis
 * Authors: Donald Conrad <dconrad@genetics.wustl.edu>,
 * Avinash Ramu <aramu@genetics.wustl.edu>
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

#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include "parser.h"
//#include "newmatap.h"
//#include "newmatio.h"
#include <Eigen/KroneckerProduct>
#include "lookup.h"

#define MIN_READ_DEPTH_INDEL 10

// parameters for the indel mutation rate model
const float INSERTION_SLOPE = -0.2994;
const float INSERTION_INTERCEPT = -22.8689;
const float DELETION_SLOPE = -0.2856;
const float DELETION_INTERCEPT = -21.9313;


using namespace std;

// Calculate DNM and Null PP
void trio_like_indel(indel_t *child, indel_t *mom, indel_t *dad, int flag,
                     lookup_table_t &tgtIndel, lookup_indel_t &lookupIndel,
                     double mu_scale, vector<hts::bcf::File> &vcfout, double pp_cutoff,
                     int RD_cutoff, int &n_site_pass, double user_indel_mrate) {
    // Read depth filter
    if(child->depth < RD_cutoff ||
            mom->depth < RD_cutoff || dad->depth < RD_cutoff) {
        return;
    }

    n_site_pass += 1;
    //Real a[3];
    Real maxlike_null, maxlike_denovo, pp_null, pp_denovo, denom;
    Matrix M(1, 3);
    Matrix C(3, 1);
    Matrix D(3, 1);
    Matrix P(3, 3);
    Matrix F(9, 3);
    Matrix L(9, 3);
    Matrix T(9, 3);
    Matrix DN(9, 3);
    Matrix PP(9, 3);
    //int i, j, k, l;
    int coor = child->pos;
    char ref_name[50];
    strcpy(ref_name, child->chr);  // Name of the reference sequence

    // this is the case where child has more than 1 indel alt allele, ignore for now
    if(strstr(child->alt, ",") != 0) {
        return;
    }

    bool is_insertion = false; // insertion/deletion event
    int len_diff = strlen(mom->ref_base) - strlen(mom->alt); // diff b/w alt, ref
    if(len_diff < 0) {
        is_insertion = true;
        len_diff = -len_diff;
    } else if(len_diff > 0) {
        is_insertion = false;
    }

    // mu_scale is the variable used to scale the mutation rate prior
    double log_indel_mrate, indel_mrate, new_indel_mrate;
    if(is_insertion) {
        log_indel_mrate = (INSERTION_SLOPE * len_diff + INSERTION_INTERCEPT);
    } else {
        log_indel_mrate = (DELETION_SLOPE * len_diff + DELETION_INTERCEPT);
    }
    indel_mrate = exp(log_indel_mrate);// antilog of log ratio

    //use user supplied indel mrate if this is not the default value.
    if(user_indel_mrate != 0) {
        indel_mrate = user_indel_mrate;
    }

    indel_mrate = mu_scale * indel_mrate;

    for(int j = 0; j < 9; j++) {
        for(int l = 0; l < 3; l++) {
            // hit is 0, 1 or 2 (number of indels in the trio config)
            new_indel_mrate = pow(indel_mrate, lookupIndel.hit(j, l));
            lookupIndel.mrate(j, l) = new_indel_mrate;
        }
    }

    //Load likelihood vectors
    for(int j = 0; j < 3; ++j) {
        M(0, j) = pow(10, -mom->lk[j] / 10.);
    }
    for(int j = 0; j < 3; ++j) {
        D(j, 0) = pow(10, -dad->lk[j] / 10.);
    }
    for(int j = 0; j < 3; ++j) {
        C(j, 0) = pow(10, -child->lk[j] / 10.);
    }

    P = kroneckerProduct(M, D);
    F = kroneckerProduct(P, C);

    // combine with transmission probs
    T = F.cwiseProduct(lookupIndel.tp);

    // combine with priors
    L = T.cwiseProduct(lookupIndel.priors);

    // combine with mutation rate
    DN = L.cwiseProduct(lookupIndel.mrate);

    // Find max likelihood of null configuration
    PP = DN.cwiseProduct(lookupIndel.norm);
    int i, j;
    maxlike_null = PP.maxCoeff(&i, &j);


    //Find max likelihood of de novo trio configuration
    PP = DN.cwiseProduct(lookupIndel.denovo);
    int k, l;
    maxlike_denovo = PP.maxCoeff(&k, &l);

    //make proper posterior probs
    denom = DN.sum();
    pp_denovo = maxlike_denovo / denom; // denovo posterior probability
    pp_null = 1 - pp_denovo; // null posterior probability

    // Check for PP cutoff
    if(pp_denovo > pp_cutoff) {

        //remove ",X" from alt, helps with VCF op.
        string alt = mom->alt;
        size_t posX = alt.find(",X");
        size_t posN = alt.find(",N");
        if(posX != std::string::npos) {
            alt.replace(posX, 2, "");
        }
        else if(posN != std::string::npos) {
            alt.replace(posN, 2, "");
        }

        cout << "DENOVO-INDEL CHILD_ID: " << child->id;       
        cout << " chr: " << ref_name;
        cout << " pos: " << coor; 
        cout << " ref: " << mom->ref_base;
        cout << " alt: " << mom->alt;
        cout << " maxlike_null: " << maxlike_null;
        cout << " pp_null: " << pp_null;
        cout << " tgt_null(child/mom/dad): " << tgtIndel[i][j];
        cout << " snpcode: " << lookupIndel.snpcode(i, j);
        cout << " code: " << lookupIndel.code(i, j);
        cout << " maxlike_dnm: " << maxlike_denovo;
        cout << " pp_dnm: " << pp_denovo;       
        cout << " tgt_dnm(child/mom/dad): " << tgtIndel[k][l];
        cout << " lookup: " << lookupIndel.code(k, l);
        cout << " flag: " << flag;
        cout << " READ_DEPTH child: " << child->depth; 
        cout << " dad: " << dad->depth;
        cout << " mom: " << mom->depth;
        cout << " MAPPING_QUALITY child: " << child->rms_mapQ;
        cout << " dad: " << dad->rms_mapQ;
        cout << " mom: " << mom->rms_mapQ;
        cout << endl;
        

        if(!vcfout.empty()) {
            auto rec = vcfout[0].InitVariant();
            rec.target(ref_name);
            rec.position(coor);
            rec.alleles(std::string(mom->ref_base) + "," + alt);
            rec.quality(0);
            rec.filter("PASS");
#ifndef NEWVCFOUT
            // Old vcf output format
            rec.info("RD_MOM", mom->depth);
            rec.info("RD_DAD", dad->depth);
            rec.info("MQ_MOM", mom->rms_mapQ);
            rec.info("MQ_DAD", dad->rms_mapQ);
            rec.info("INDELcode", static_cast<float>(lookupIndel.snpcode(i, j)));
            rec.samples("NULL_CONFIG(child/mom/dad)", std::vector<std::string> {tgtIndel[i - 1][j - 1]});
            rec.samples("PP_NULL", std::vector<float> {static_cast<float>(pp_null)});
            rec.samples("DNM_CONFIG(child/mom/dad)", std::vector<std::string> {tgtIndel[k - 1][l - 1]});
            rec.samples("PP_DNM", std::vector<float> {static_cast<float>(pp_denovo)});
            rec.samples("RD", std::vector<int32_t> {child->depth});
            rec.samples("MQ", std::vector<int32_t> {child->rms_mapQ});
#else

            // Newer, more accurate VCF output format
            rec.info("INDELcode", static_cast<float>(lookupIndel.snpcode(i, j)));
            rec.info("PP_NULL", static_cast<float>(pp_null));
            rec.info("PP_DNM", static_cast<float>(pp_denovo));
            rec.samples("RD", std::vector<int32_t> {child->depth, mom->depth, dad->depth});
            rec.samples("MQ", std::vector<int32_t> {child->rms_mapQ, mom->rms_mapQ, dad->rms_mapQ});
            std::vector<std::string> configs;
            boost::split(configs, tgtIndel[i - 1][j - 1], boost::is_any_of("/"));
            rec.samples("NULL_CONFIG", configs);
            boost::split(configs, tgtIndel[k - 1][l - 1], boost::is_any_of("/"));
            rec.samples("DNM_CONFIG", configs);
#endif
            vcfout[0].WriteRecord(rec);
        }
    }
}
