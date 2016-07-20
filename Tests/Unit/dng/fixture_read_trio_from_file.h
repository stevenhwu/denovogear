/*
 * Copyright (c) 2016 Steven H. Wu
 * Authors:  Steven H. Wu <stevenwu@asu.edu>
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
#ifndef DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H
#define DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H

#include <fstream>

#include <dng/task/call.h>
#include <dng/hts/bcf.h>

using namespace dng;

// Helper function that mimics boost::istream_range
template<class Elem, class Traits> inline boost::iterator_range<
        std::istreambuf_iterator<Elem, Traits>> istreambuf_range(
        std::basic_istream<Elem, Traits> &in) {
    return boost::iterator_range<std::istreambuf_iterator<Elem, Traits>>(
            std::istreambuf_iterator<Elem, Traits>(in),
            std::istreambuf_iterator<Elem, Traits>());
}


struct ReadTrioFromFile {

    std::string fixture;

    dng::io::Pedigree io_pedigree;
    dng::ReadGroups rgs;

    typedef dng::task::Call task_type;
    struct arg_t : public task_type::argument_type {
        bool help;
        bool version;
        std::string arg_file;

        std::string run_name;
        std::string run_path;
    } arg;

    ReadTrioFromFile(std::string s = "ReadTrioFromFile") : fixture(s) {
        BOOST_TEST_MESSAGE("set up fixture: " << fixture);

        po::options_description ext_desc, int_desc;
        po::positional_options_description pos_desc;
        po::variables_map vm;

        int argc=4;
        char *argv[argc];
        argv[0] = (char*) "test";
        argv[1] = (char*) "-p";

        std::string ped_filename (TESTDATA_DIR);
        ped_filename.append("/sample_5_3/ceu.ped");
        argv[2] = (char*) ped_filename.data();

        std::string vcf_filename = TESTDATA_DIR;
        vcf_filename.append("/sample_5_3/test1.vcf");
        argv[3] = (char*) vcf_filename.data();

        add_app_args(ext_desc, static_cast<typename task_type::argument_type &>(arg));
        int_desc.add_options()
                ("input", po::value< std::vector<std::string> >(&arg.input), "input files")
                ;
        int_desc.add(ext_desc);
        pos_desc.add("input", -1);
        po::store(po::command_line_parser(argc, argv)
                          .options(int_desc).positional(pos_desc).run(), vm);
        po::notify(vm);

        // Parse pedigree from file

        std::ifstream ped_file(arg.ped);
        io_pedigree.Parse(istreambuf_range(ped_file));


        std::vector<hts::File> indata;
        std::vector<hts::bcf::File> bcfdata;
        for (auto &&str : arg.input) {
            indata.emplace_back(str.c_str(), "r");
            if (indata.back().is_open()) {
                continue;
            }
            throw std::runtime_error("unable to open input file '" + str + "'.");
        }
        bcfdata.emplace_back(std::move(indata[0]));
        rgs.ParseSamples(bcfdata[0]);

    }

    ~ReadTrioFromFile() {
        BOOST_TEST_MESSAGE("tear down fixture: " << fixture);
    }


};


#endif //DENOVOGEAR_FIXTURE_READ_TRIO_FROM_FILE_H
