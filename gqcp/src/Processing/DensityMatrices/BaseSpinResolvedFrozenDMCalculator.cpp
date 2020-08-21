// This file is part of GQCG-GQCP.
//
// Copyright (C) 2017-2020  the GQCG developers
//
// GQCG-GQCP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-GQCP is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-GQCP.  If not, see <http://www.gnu.org/licenses/>.

#include "Processing/DensityMatrices/BaseSpinResolvedFrozenDMCalculator.hpp"

#include "Utilities/linalg.hpp"


namespace GQCP {

/*
 *  CONSTRUCTORS
 */

/**
 *  @param rdm_builder                  shared pointer to active (non-frozen core) RDM builder
 *  @param X                            the number of frozen orbitals
 */
BaseSpinResolvedFrozenDMCalculator::BaseSpinResolvedFrozenDMCalculator(const std::shared_ptr<BaseSpinResolvedDMCalculator> rdm_builder, const size_t X) :
    BaseSpinResolvedDMCalculator(),
    active_rdm_builder {std::move(rdm_builder)},
    X {X} {}


/*
 * PUBLIC OVERRIDDEN METHODS
 */

/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 1-RDMs given a coefficient vector
 */
SpinResolvedOneDM<double> BaseSpinResolvedFrozenDMCalculator::calculate1RDMs(const VectorX<double>& x) const {

    auto K = this->onvBasis()->numberOfOrbitals();

    OneDM<double> D_aa = OneDM<double>::Zero(K, K);
    OneDM<double> D_bb = OneDM<double>::Zero(K, K);

    auto K_active = K - this->X;

    // Set the values for the frozen orbital
    for (size_t i = 0; i < this->X; i++) {
        D_aa(i, i) = 1;
        D_bb(i, i) = 1;
    }

    SpinResolvedOneDM<double> sub_1rdms = this->active_rdm_builder->calculate1RDMs(x);

    // Incorporate the submatrices from the active space
    D_aa.block(this->X, this->X, K_active, K_active) += sub_1rdms.alpha();
    D_bb.block(this->X, this->X, K_active, K_active) += sub_1rdms.beta();

    return SpinResolvedOneDM<double>(D_aa, D_bb);
};


/**
 *  @param x        the coefficient vector representing the wave function
 *
 *  @return all 2-RDMs given a coefficient vector
 */
SpinResolvedTwoDM<double> BaseSpinResolvedFrozenDMCalculator::calculate2RDMs(const VectorX<double>& x) const {

    auto K = this->onvBasis()->numberOfOrbitals();

    TwoDM<double> d_aaaa {K};
    d_aaaa.setZero();

    TwoDM<double> d_aabb {K};
    d_aabb.setZero();

    TwoDM<double> d_bbaa {K};
    d_bbaa.setZero();

    TwoDM<double> d_bbbb {K};
    d_bbbb.setZero();


    SpinResolvedOneDM<double> one_rdms = this->active_rdm_builder->calculate1RDMs(x);
    auto D_aa = one_rdms.alpha();
    auto D_bb = one_rdms.beta();

    // Implement frozen RDM formulas
    for (size_t p = 0; p < this->X; p++) {  // iterate over frozen orbitals

        // RDM Overlap between frozen and active space:
        //      frozen orbital indices (p) must always have one annihilation and one creation index (always occupied)
        //      values are dictated by the 'active' orbital indices and correspond to that of the active 1RDMs
        //      Hence we start adding the 1RDMs starting from index 'X' the number frozen orbitals
        d_aaaa.addBlock<0, 1>(D_aa, this->X, this->X, p, p);
        d_aaaa.addBlock<2, 3>(D_aa, p, p, this->X, this->X);
        d_aaaa.addBlock<2, 1>(-D_aa, p, this->X, this->X, p);
        d_aaaa.addBlock<3, 0>(-D_aa, this->X, p, p, this->X);

        d_bbbb.addBlock<0, 1>(D_bb, this->X, this->X, p, p);
        d_bbbb.addBlock<2, 3>(D_bb, p, p, this->X, this->X);
        d_bbbb.addBlock<2, 1>(-D_bb, p, this->X, this->X, p);
        d_bbbb.addBlock<3, 0>(-D_bb, this->X, p, p, this->X);

        d_aabb.addBlock<2, 3>(D_bb, p, p, this->X, this->X);
        d_aabb.addBlock<0, 1>(D_aa, this->X, this->X, p, p);

        d_bbaa.addBlock<2, 3>(D_aa, p, p, this->X, this->X);
        d_bbaa.addBlock<0, 1>(D_bb, this->X, this->X, p, p);


        // Set the values for the frozen orbital
        d_bbaa(p, p, p, p) = 1;
        d_aabb(p, p, p, p) = 1;

        for (size_t q = p + 1; q < this->X; q++) {  // iterate over frozen orbitals

            d_aaaa(p, p, q, q) = 1;
            d_aaaa(q, q, p, p) = 1;
            d_aaaa(p, q, q, p) = -1;
            d_aaaa(q, p, p, q) = -1;

            d_bbbb(p, p, q, q) = 1;
            d_bbbb(q, q, p, p) = 1;
            d_bbbb(p, q, q, p) = -1;
            d_bbbb(q, p, p, q) = -1;

            d_aabb(p, p, q, q) = 1;
            d_bbaa(p, p, q, q) = 1;

            d_aabb(q, q, p, p) = 1;
            d_bbaa(q, q, p, p) = 1;
        }
    }


    // Incorporate the 2-RDM subblocks into the total 2RDMs
    SpinResolvedTwoDM<double> sub_2rdms = this->active_rdm_builder->calculate2RDMs(x);

    d_aaaa.addBlock(sub_2rdms.two_rdm_aaaa, this->X, this->X, this->X, this->X);
    d_bbbb.addBlock(sub_2rdms.two_rdm_bbbb, this->X, this->X, this->X, this->X);
    d_aabb.addBlock(sub_2rdms.two_rdm_aabb, this->X, this->X, this->X, this->X);
    d_bbaa.addBlock(sub_2rdms.two_rdm_bbaa, this->X, this->X, this->X, this->X);

    return SpinResolvedTwoDM<double>(d_aaaa, d_aabb, d_bbaa, d_bbbb);
};


}  // namespace GQCP
