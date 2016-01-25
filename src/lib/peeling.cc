/*
 * Copyright (c) 2014-2015 Reed A. Cartwright
 * Authors:  Reed A. Cartwright <reed@cartwrig.ht>
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

#include <dng/peeling.h>
#include <iostream>

bool DEBUG_PEELING = false;
//#define DEBUG_PEELING true;
//#ifdef DEBUG_PEELING
//    std::cout << "==" << __FUNCTION__ << "==" << std::endl;
//#endif

[[deprecated]] dng::GenotypeArray dng::peel::multiply_upper_lower(workspace_t &work, std::size_t index){

    return (work.upper[index] * work.lower[index]);
}

[[deprecated]]
dng::GenotypeArray dng::peel::multiply_lower_upper(workspace_t &work, std::size_t index){
    return multiply_upper_lower(work, index);
}


dng::PairedGenotypeArray dng::peel::kroneckerProductDadMom(workspace_t &work, std::size_t dad, std::size_t mom){

    dng::PairedGenotypeArray result = kroneckerProduct(
            work.multiply_upper_lower(dad).matrix(),
            work.multiply_upper_lower(mom).matrix() );
    return result;

}


// Family Order: Father, Mother, Child1, Child2, ...
dng::PairedGenotypeArray dng::peel::sum_over_children(workspace_t &work, const family_members_t &family,
                                                      const TransitionVector &mat) {
    PairedGenotypeArray buffer = sum_over_children(work, family, mat, 2);
    return buffer;
}

// Family Order: Father, Mother, Child1, Child2, ...
dng::PairedGenotypeArray dng::peel::sum_over_children(workspace_t &work, const family_members_t &family,
                                                      const TransitionVector &mat, int first_child_index) {
//    assert(family.size() >= 3);
//    assert(family.size() >= first_child_index);
//    PairedGenotypeArray buffer = PairedGenotypeArray::Ones(100,1);
    PairedGenotypeArray buffer = (mat[family[first_child_index]] * work.lower[family[first_child_index]].matrix()).array();
    for(std::size_t i = first_child_index+1; i < family.size(); ++i) {
        buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    return buffer;
}


// Family Order: Parent, Child
dng::GenotypeArray dng::peel::up_core(workspace_t &work, const family_members_t &family,
                                      const TransitionVector &mat) {

    assert(family.size() == 2);
    auto child = family[1];
    dng::GenotypeArray geno_array = (mat[child] * work.lower[child].matrix()).array();
    return geno_array;
}


// Family Order: Father, Mother, Child1, Child2, ...
dng::GenotypeArray dng::peel::to_parent_core(workspace_t &work, const family_members_t &family,
                                             const TransitionVector &mat, const Parents to_parent) {

    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];

    work.paired_buffer = sum_over_children(work, family, mat, 2);
    work.paired_buffer.resize(10, 10);

    GenotypeArray parent_array;
    if(to_parent == Parents::Father) {
        parent_array = (work.paired_buffer.matrix() * work.multiply_upper_lower(mom).matrix() ).array();
    }

    else if (to_parent == Parents::Mother){
        parent_array = (work.paired_buffer.matrix().transpose() * work.multiply_upper_lower(dad).matrix() ).array();
    }
    else{
        //ERROR: Not possible until we add more enum!!
    }


    return parent_array;
}


// Family Order: Parent, Child
void dng::peel::down(workspace_t &work, const family_members_t &family,
                     const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }

    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

    work.upper[child] = (mat[child].transpose() * work.multiply_upper_lower(parent).matrix() ).array();
}

// Family Order: Parent, Child
void dng::peel::down_fast(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }

    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.upper[child] = (mat[child].transpose() * work.upper[parent].matrix()).array();

}

// Family Order: Parent, Child
void dng::peel::up(workspace_t &work, const family_members_t &family,
                   const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }

    auto parent = family[0];
    work.lower[parent] *= up_core(work, family, mat);
}

// Family Order: Parent, Child
void dng::peel::up_fast(workspace_t &work, const family_members_t &family,
                        const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    auto parent = family[0];
    work.lower[parent] = up_core(work, family, mat);
}



// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    auto dad = family[0];
    work.lower[dad] *= to_parent_core(work, family, mat, Parents::Father);

}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father_fast(workspace_t &work, const family_members_t &family,
                               const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    auto dad = family[0];
    work.lower[dad] =  to_parent_core(work, family, mat, Parents::Father);

}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    auto mom = family[1];
    work.lower[mom] *= to_parent_core(work, family, mat, Parents::Mother);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_fast(workspace_t &work,
                               const family_members_t &family, const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    auto mom = family[1];
    work.lower[mom] = to_parent_core(work, family, mat, Parents::Mother);
}




// Family Order: Father, Mother, Child, Child2, ....
void dng::peel::to_child(workspace_t &work, const family_members_t &family,
                         const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() >= 4);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Parents
    work.paired_buffer = kroneckerProductDadMom(work, dad, mom).array();
    work.paired_buffer *= sum_over_children(work, family, mat, 3);
    work.upper[child] = (mat[child].transpose() * work.paired_buffer.matrix()).array();

}

// Family Order: Father, Mother, CHild
void dng::peel::to_child_fast(workspace_t &work, const family_members_t &family,
                              const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() == 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Parents
    work.paired_buffer = kroneckerProductDadMom(work, dad, mom).array();

    work.upper[child] = (mat[child].transpose() * work.paired_buffer.matrix()).array();
}

// Family Order: Parent, Child
void dng::peel::down_reverse(workspace_t &work, const family_members_t &family,
                             const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

    work.super[child] = work.upper[parent];
    work.lower[parent] *= (mat[child] * work.lower[child].matrix()).array();
}

// Family Order: Parent, Child
// prevent divide by zero errors by adding a minor offset to one of the calculations
void dng::peel::up_reverse(workspace_t &work, const family_members_t &family,
                           const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

    work.super[child] = work.upper[parent] * (work.lower[parent] /
                        ((mat[child] * work.lower[child].matrix()).array() +
                         DNG_INDIVIDUAL_BUFFER_MIN));
    work.upper[child] = (mat[child].transpose() *
                         work.super[child].matrix()).array();
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father_reverse(workspace_t &work,
                                  const family_members_t &family, const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    // work.paired_buffer will contain P(child data | mom & dad)
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    // Calculate P(mom-only data & mom = g)
    IndividualVector::value_type mom_v = work.upper[mom] * work.lower[mom];

    // Calculate P(dad-only data & dad = g)
    work.paired_buffer.resize(10, 10);
    IndividualVector::value_type dad_v = work.upper[dad] * (work.lower[dad] /
                                         ((work.paired_buffer.matrix() * mom_v.matrix()).array() +
                                          DNG_INDIVIDUAL_BUFFER_MIN));

    // Calculate P(dependent data | mom = g)
    work.lower[mom] *= (work.paired_buffer.matrix().transpose() *
                        dad_v.matrix()).array();
    work.paired_buffer.resize(100, 1);

    // Calculate P(data & dad & mom)
    work.paired_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();

    for(std::size_t i = 2; i < family.size(); ++i) {
        auto child = family[i];
        work.super[child] = work.paired_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_reverse(workspace_t &work,
                                  const family_members_t &family, const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    IndividualVector::value_type dad_v = work.upper[dad] * work.lower[dad];

    work.paired_buffer.resize(10, 10);
    IndividualVector::value_type mom_v = work.upper[mom] * (work.lower[mom] /
                                         ((work.paired_buffer.matrix().transpose() * dad_v.matrix()).array() +
                                          DNG_INDIVIDUAL_BUFFER_MIN));
    work.lower[dad] *= (work.paired_buffer.matrix() * mom_v.matrix()).array();
    work.paired_buffer.resize(100, 1);

    work.paired_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();

    for(std::size_t i = 2; i < family.size(); ++i) {
        auto child = family[i];
        work.super[child] = work.paired_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_child_reverse(workspace_t &work,
                                 const family_members_t &family, const TransitionVector &mat) {
    if(DEBUG_PEELING){ std::cout << "==" << __FUNCTION__ << "==" << std::endl; }
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Sum over children
    work.paired_buffer = (mat[child] * work.lower[child].matrix()).array();
    for(std::size_t i = 3; i < family.size(); i++) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    // Update Parents
    IndividualVector::value_type dad_v = work.upper[dad] * work.lower[dad];
    IndividualVector::value_type mom_v = work.upper[mom] * work.lower[mom];
    work.paired_buffer.resize(10, 10);
    work.lower[dad] *= (work.paired_buffer.matrix() * mom_v.matrix()).array();
    work.lower[mom] *= (work.paired_buffer.matrix().transpose() *
                        dad_v.matrix()).array();
    work.paired_buffer.resize(100, 1);

    // Update Siblings
    work.paired_buffer *= kroneckerProduct(dad_v.matrix(), mom_v.matrix()).array();
    work.super[child] = work.paired_buffer / ((mat[child] *
                        work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
    for(std::size_t i = 3; i < family.size(); ++i) {
        child = family[i];
        work.super[child] = work.paired_buffer / ((mat[child] *
                            work.lower[child].matrix()).array() + DNG_INDIVIDUAL_BUFFER_MIN);
        work.upper[child] = mat[child].transpose() * work.super[child].matrix();
    }
}





// Family Order: Father, Mother, Child1, Child2, ...
//[[deprecated]]
dng::GenotypeArray dng::peel::to_father_core(workspace_t &work, const family_members_t &family,
                                             const TransitionVector &mat) {

    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];

//    work.paired_buffer = sum_over_children(work, family, mat);
//    // Include Mom
//    work.paired_buffer.resize(10, 10);
//    PairedGenotypeArray other_parent = multiply_upper_lower(work, mom);
//    GenotypeArray dad_array = (work.paired_buffer.matrix() * other_parent.matrix() ).array();
    GenotypeArray dad_array = to_parent_core(work, family, mat, Parents::Father);

    return dad_array;
}
// Family Order: Father, Mother, Child1, Child2, ...
[[deprecated]]
dng::GenotypeArray dng::peel::to_mother_core(workspace_t &work, const family_members_t &family,
                                             const TransitionVector &mat) {

    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];

//    work.paired_buffer = sum_over_children(work, family, mat);
//    // Include Mom
//    work.paired_buffer.resize(10, 10);
//    PairedGenotypeArray other_parent = multiply_upper_lower(work, dad);
//    GenotypeArray mom_array = (work.paired_buffer.matrix().transpose() * other_parent.matrix() ).array();
    GenotypeArray  mom_array = to_parent_core(work, family, mat, Parents::Mother);
    return mom_array;
}



