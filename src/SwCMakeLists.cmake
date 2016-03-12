#####################################################
####### STEVEN CUSTOM CMAKE
INCLUDE_DIRECTORIES(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")

SET(DNG_LIBRARIES
        Threads::Threads
        Boost::PROGRAM_OPTIONS
        Boost::FILESYSTEM
        Boost::SYSTEM
        EIGEN3::EIGEN3
        HTSLIB::HTSLIB)


## A lazy hack! which might be against whatever CMAKE recommend
FILE(GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/lib/*cc")
LIST(REMOVE_ITEM SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/lib/call.cc")

ADD_EXECUTABLE(sw_dng_call
        dng-call_test.cc
        vt/calculateProbs.cc
        ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(sw_dng_call ${DNG_LIBRARIES})

ADD_EXECUTABLE(sw_test_find_mutation
        vt/find_mutation_test.cc
        ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(sw_test_find_mutation  ${DNG_LIBRARIES})

IF (DEVEL_MODE)
    TARGET_LINK_LIBRARIES(sw_test_find_mutation Boost::TIMER)
ENDIF ()

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

SET(DEBUG_VERBOSE 0)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g  -O0")
IF(DEBUG_VERBOSE GREATER 0)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wvla -Wall -Wextra -pedantic") #-pedantic-errors
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wextra -Waggregate-return -Wcast-align") #-pedantic-errors
ENDIF()
IF(DEBUG_VERBOSE GREATER 1)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wcast-qual  -Wchar-subscripts  -Wcomment " )
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wdisabled-optimization ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wfloat-equal  -Wformat  -Wformat=2 ")#-Werror
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wformat-nonliteral -Wformat-security  ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wformat-y2k ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wimplicit  -Wimport  -Winit-self  -Winline ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Winvalid-pch   ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wlong-long -Wmissing-braces ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wmissing-field-initializers -Wmissing-format-attribute   ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wmissing-include-dirs ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wpacked -Wparentheses  -Wpointer-arith ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wredundant-decls -Wreturn-type ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wsequence-point   -Wsign-compare  -Wstack-protector ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wswitch-enum -Wtrigraphs  -Wuninitialized ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunknown-pragmas  -Wunused -Wunreachable-code")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunused-function  -Wunused-label  -Wunused-parameter ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wunused-value  -Wunused-variable  -Wvariadic-macros ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wvolatile-register-var  -Wwrite-strings ")
ENDIF()
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wextra-semi ")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wconversion -Wsign-conversion ")

##Everything
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic") #forget about c++98

##turn off some for now, should add back to final code

IF(DEBUG_VERBOSE LESS 2)
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-deprecated-declarations -Wno-deprecated")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unused-private-field -Wno-unused-variable -Wno-unused-private-field")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unused-parameter -Wno-unused-function")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-unreachable-code")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-conversion -Wno-sign-conversion -Wno-missing-noreturn") #Way too many from Eigen/Boost
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-float-equal -Wno-sign-compare")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-shadow")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-padded")

    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-global-constructors -Wno-exit-time-destructors")
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}  -Wno-weak-vtables -Wno-missing-prototypes")


ENDIF()







##

##

