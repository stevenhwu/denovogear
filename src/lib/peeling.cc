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

dng::GenotypeArray dng::peel::multiply_upper_lower(workspace_t &work, int index){

    return (work.upper[index] * work.lower[index]);
}

dng::GenotypeArray dng::peel::multiply_lower_upper(workspace_t &work, int index){

    return (work.lower[index] * work.upper[index]);
}

// Family Order: Parent, Child
void dng::peel::down(workspace_t &work, const family_members_t &family,
                     const TransitionVector &mat) {
    std::cout << " down " << std::endl;
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.upper[child] = (mat[child].transpose() * (work.upper[parent] *
                         work.lower[parent]).matrix()).array();
}

// Family Order: Parent, Child
void dng::peel::down_fast(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    std::cout << " down fast" << std::endl;
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.upper[child] = (mat[child].transpose() *
                         work.upper[parent].matrix()).array();
}

dng::PairedGenotypeArray dng::peel::sum_over_child(workspace_t &work, const family_members_t &family,
                                        const TransitionVector &mat) {
    assert(family.size() >= 3);
    auto first_child = family[2];

    PairedGenotypeArray buffer = (mat[first_child] * work.lower[first_child].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    return buffer;
}

dng::GenotypeArray dng::peel::up_core(workspace_t &work, const family_members_t &family,
    const TransitionVector &mat) {

    assert(family.size() == 2);
//    auto parent = family[0];
    auto child = family[1];

    dng::GenotypeArray geno_array = (mat[child] * work.lower[child].matrix()).array();
    return geno_array;
}
// Family Order: Parent, Child
void dng::peel::up(workspace_t &work, const family_members_t &family,
                   const TransitionVector &mat) {
    std::cout << " up " << std::endl;
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];

//    std::cout << work.lower[parent] << "\n" << std::endl;
//    std::cout << work.lower[child] << "\n" << std::endl;
//    auto a = (mat[child] * work.lower[child].matrix());
//    auto b = up_core(work, family, mat);
//    std::cout << "A:\n" << a << std::endl;
//    auto c = a;
//    auto d = a.array();
//    std::cout << "B:\n" << b << std::endl;
//    std::cout << "C:\n" << c << std::endl;
//    std::cout << "D:\n" << d << std::endl;
//    work.lower[parent] *= a;

    work.lower[parent] *= (mat[child] * work.lower[child].matrix()).array();
//    work.lower[parent] *= up_core(work, family, mat); //TODO: replace with this


//    std::cout << work.lower[parent] << "\n" << std::endl;
//    work.lower[parent] = (mat[child] * work.lower[child].matrix()).array();
//    std::cout << work.lower[parent] << "\n" << std::endl;
}

// Family Order: Parent, Child
void dng::peel::up_fast(workspace_t &work, const family_members_t &family,
                        const TransitionVector &mat) {
    std::cout << " up fast " << std::endl;
    assert(family.size() == 2);
    auto parent = family[0];
    auto child = family[1];
    work.lower[parent] = (mat[child] * work.lower[child].matrix()).array();
//    work.lower[parent] = up_core(work, family, mat); //TODO: replace with this


//    std::cout << work.lower[parent] << std::endl;
//    for (int i = 0; i < 5; ++i) {
//        std::cout << "lower: "<<i<<"\n"<<work.lower[i] << "\n" << std::endl;
//    }
//    for (int i = 0; i < 5; ++i) {
//        std::cout << "upper: "<<i<<"\n"<<work.upper[i] << "\n" <<std::endl;
//    }
}


// Family Order: Father, Mother, Child1, Child2, ...
dng::GenotypeArray dng::peel::to_father_core(workspace_t &work, const family_members_t &family,
                                        const TransitionVector &mat) {
    std::cout << " to father core: " << family[0] << "\t" << family[1] << "\t" << family[2] << std::endl;
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];

    work.paired_buffer = sum_over_child(work, family, mat); //TODO: replace with this
    // Include Mom
    work.paired_buffer.resize(10, 10);
    PairedGenotypeArray mom_array = multiply_upper_lower(work, mom); //TODO: replace with this
    GenotypeArray dad_array = (work.paired_buffer.matrix() * mom_array.matrix() ).array();//TODO: replace with this
    return dad_array;
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    std::cout << " to father: " << family[0] << "\t" << family[1] << "\t" << family[2] << std::endl;
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
//    work.paired_buffer = sum_over_child(work, family, mat);

    // Include Mom
    work.paired_buffer.resize(10, 10);

//    std::cout << work.paired_buffer.sum() << "\n" << std::endl;
//    std::cout << (work.upper[mom] * work.lower[mom]).sum() << std::endl;

    work.lower[dad] *= (work.paired_buffer.matrix() * (work.upper[mom] *
                        work.lower[mom]).matrix()).array();


    std::cout << ((work.paired_buffer.matrix() * (work.upper[mom] * work.lower[mom]).matrix()).array()) << std::endl;
    std::cout<< work.lower[dad] << std::endl;

    work.paired_buffer.resize(100, 1); //Might not need this, from the website: Assignment is the action of copying a matrix into another, using operator=. Eigen resizes the matrix on the left-hand side automatically so that it matches the size of the matrix on the right-hand size. For example:


//    work.lower[dad] *= to_father_core(work, family, mat); //TODO: replace with this
}


// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_father_fast(workspace_t &work,
                               const family_members_t &family, const TransitionVector &mat) {
    std::cout << " to father fast: " << family[0] << "\t" << family[1] << "\t" << family[2] << std::endl;
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
//    work.paired_buffer = sum_over_child(work, family, mat); //TODO: replace with this

    // Include Mom
    work.paired_buffer.resize(10, 10);
    work.lower[dad] = (work.paired_buffer.matrix() * (work.upper[mom] *
                       work.lower[mom]).matrix()).array();

    std::cout <<  work.paired_buffer.sum() << std::endl;
    std::cout <<  work.upper[mom].sum() << std::endl;
    std::cout <<  work.lower[mom].sum() << std::endl;
    std::cout <<  work.lower[dad].sum() << std::endl;

    work.paired_buffer.resize(100, 1);


}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother(workspace_t &work, const family_members_t &family,
                          const TransitionVector &mat) {
    std::cout << " to mother " << std::endl;
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    //    work.paired_buffer = sum_over_child(work, family, mat); //TODO: replace with this

    // Include Dad
    work.paired_buffer.resize(10, 10);
    work.lower[mom] *= (work.paired_buffer.matrix().transpose() * (work.upper[dad] *
                        work.lower[dad]).matrix()).array();
    work.paired_buffer.resize(100, 1);
}

// Family Order: Father, Mother, Child1, Child2, ...
void dng::peel::to_mother_fast(workspace_t &work,
                               const family_members_t &family, const TransitionVector &mat) {
    std::cout << " to mother fast " << std::endl;
    assert(family.size() >= 3);
    auto dad = family[0];
    auto mom = family[1];
    // Sum over children
    work.paired_buffer = (mat[family[2]] * work.lower[family[2]].matrix()).array();
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }
    //    work.paired_buffer = sum_over_child(work, family, mat); //TODO: replace with this

    // Include Dad
    work.paired_buffer.resize(10, 10);
    work.lower[mom] = (work.paired_buffer.matrix().transpose() * (work.upper[dad] *
                       work.lower[dad]).matrix()).array();
    work.paired_buffer.resize(100, 1);
}

// Family Order: Father, Mother, Child, Child2, ....
void dng::peel::to_child(workspace_t &work, const family_members_t &family,
                         const TransitionVector &mat) {
    std::cout << " to child " << std::endl;
    assert(family.size() >= 4);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    // Parents
    work.paired_buffer = kroneckerProduct((work.lower[dad] *
                                           work.upper[dad]).matrix(),
                                          (work.lower[mom] * work.upper[mom]).matrix()).array();
    // Sum over fullsibs
    for(std::size_t i = 3; i < family.size(); ++i) {
        work.paired_buffer *= (mat[family[i]] * work.lower[family[i]].matrix()).array();
    }

    work.upper[child] = (mat[child].transpose() *
                         work.paired_buffer.matrix()).array();
}

// Family Order: Father, Mother, CHild
void dng::peel::to_child_fast(workspace_t &work, const family_members_t &family,
                              const TransitionVector &mat) {
    std::cout << " to child fast " << std::endl;
    assert(family.size() == 3);
    auto dad = family[0];
    auto mom = family[1];
    auto child = family[2];
    work.upper[child] = (mat[child].transpose() * kroneckerProduct(
                             (work.lower[dad] * work.upper[dad]).matrix(),
                             (work.lower[mom] * work.upper[mom]).matrix())).array();
}

// Family Order: Parent, Child
void dng::peel::down_reverse(workspace_t &work, const family_members_t &family,
                             const TransitionVector &mat) {
    std::cout << " down_reverse " << std::endl;
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
    std::cout << " up reverse" << std::endl;
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
    std::cout << "to father  reverse" << std::endl;
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
    std::cout << "to mother  reverse" << std::endl;
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
    std::cout << "to child  reverse" << std::endl;
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
