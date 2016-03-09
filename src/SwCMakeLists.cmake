
#####################################################
####### STEVEN CUSTOM CMAKE
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/utils")
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")
ADD_EXECUTABLE(sw_dng_call
        sw_dng_call.cc
        lib/likelihood.cc lib/newick.cc lib/pedigree.cc lib/peeling.cc
        lib/mutation.cc lib/stats.cc

        vt/calculateProbs.cc
        utils/vcf_utils.h
        utils/assert_utils.h
        lib/find_mutation.cc
        lib/mutation_stats.cc
        lib/workspace.cc
        )

TARGET_LINK_LIBRARIES(sw_dng_call
        HTSLIB::HTSLIB Threads::Threads
        EIGEN3::EIGEN3
        Boost::PROGRAM_OPTIONS Boost::FILESYSTEM Boost::SYSTEM
        )


ADD_EXECUTABLE(sw_test_find_mutation 
        vt/find_mutation_test.cc
        lib/pedigree.cc lib/peeling.cc lib/likelihood.cc lib/newick.cc lib/mutation.cc lib/stats.cc
        lib/pedigree_v2.cc

        utils/vcf_utils.h
        utils/assert_utils.h
        lib/find_mutation.cc
        lib/mutation_stats.cc
        lib/workspace.cc
        )

TARGET_LINK_LIBRARIES(sw_test_find_mutation
        Threads::Threads
        Boost::PROGRAM_OPTIONS
        Boost::FILESYSTEM
        Boost::SYSTEM
        Boost::UNIT_TEST_FRAMEWORK

        EIGEN3::EIGEN3
        HTSLIB::HTSLIB

        )


IF(DEVEL_MODE)
  TARGET_LINK_LIBRARIES(sw_test_find_mutation Boost::TIMER)
ENDIF()

#################################################################################
## Compile unit tests
#include_directories(../src/include)
#
#foreach (test vt/peeling_test.cpp)
#    get_filename_component(testName ${test} NAME_WE)
#    add_executable(${testName} EXCLUDE_FROM_ALL ${test}
#            lib/peeling.cc
##            vt/find_mutation_getter.h
##            lib/likelihood.cc lib/newick.cc lib/pedigree.cc
##            lib/mutation.cc lib/stats.cc
##            vt/vcf_helper.h
##            vt/assert_helper.h
#            )
#    target_link_libraries(${testName}
#            Threads::Threads
#            Boost::PROGRAM_OPTIONS
#            Boost::FILESYSTEM
#            Boost::SYSTEM
#            Boost::UNIT_TEST_FRAMEWORK
#            EIGEN3::EIGEN3
#            HTSLIB::HTSLIB
#            )
#    add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
#    add_test(${testName}_run ${testName})
#    set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
#
#endforeach (test)

#
## Make a test for each unit test file located in test/ dir
#enable_testing()
##file(GLOB_RECURSE TESTS test_*[cc|cpp])
##foreach(test ${TESTS})
#SET(testSrc "vt/find_mutation_test.cpp")
#  get_filename_component(testName ${testSrc} NAME_WE)
#  add_executable(${testName} ${testSrc}
#          lib/likelihood.cc lib/newick.cc lib/pedigree.cc lib/peeling.cc
#          lib/mutation.cc lib/stats.cc
#          vt/vcf_helper.h vt/find_mutation.h)
#  target_link_libraries(${testName}
#          Threads::Threads
#          Boost::PROGRAM_OPTIONS
#          Boost::FILESYSTEM
#          Boost::SYSTEM
#          Boost::UNIT_TEST_FRAMEWORK
#          EIGEN3::EIGEN3
#          HTSLIB::HTSLIB
#          )
#  add_test(${testName}_build "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target ${testName})
#  add_test(${testName}_run ${testName})
#  set_tests_properties(${testName}_run PROPERTIES DEPENDS ${testName}_build)
#
##  set_target_properties(${testName} PROPERTIES
##            RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/testBin)
##  add_test(NAME ${testName}A
##           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/testBin
##           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/testBin/${testName} )
##endforeach(test)