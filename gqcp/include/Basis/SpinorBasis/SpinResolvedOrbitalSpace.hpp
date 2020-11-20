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

#pragma once


#include "Basis/SpinorBasis/OrbitalSpace.hpp"
#include "QuantumChemical/SpinResolvedBase.hpp"


namespace GQCP {


/**
 *  A type that encapsulates alpha and beta Orbital Spaces.
 */
class SpinResolvedOrbitalSpace:
    public SpinResolvedBase<OrbitalSpace, SpinResolvedOrbitalSpace> {


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SpinResolvedBase`'s constructors.
    using SpinResolvedBase<OrbitalSpace, SpinResolvedOrbitalSpace>::SpinResolvedBase;


    /*
     *  MARK: Public methods
     */

    /**
     *  Create an implicit mathematical object that can serve as the representation of a object with the given occupation types, from the dense representation of the matrix slice.
     *  Each occupation type can belong to either the alpha or beta spin component, initializing a mixed spin mathematical object.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit matrix
     * 
     *  @param row_type                     the spinor occupation type for the rows
     *  @param row_spin                     the spin component for the rows
     *  @param column_type                  the spinor occupation type for the columns
     *  @param column_spin                  the spin component for the columns
     *  @param M                            the dense representation of the matrix slice
     * 
     *  @return an implicit matrix slice, according to the given occupation types
     */
    template <typename Scalar>
    ImplicitMatrixSlice<Scalar> createMixedRepresentableObjectFor(const OccupationType row_type, const Spin row_spin,
                                                                  const OccupationType column_type, const Spin column_spin, const MatrixX<Scalar>& M) const {

        // Prepare the necessary members for ImplicitMatrixSlice.
        // The spin component of each axis determines from which of the two internal OrbitalSpace members the orbitals of the chosen occupation type are used.
        // Determine the spin component and dimension of the rows.
        std::vector<size_t> row_indices;
        if (row_spin == Spin::alpha) {
            row_indices = this->alpha().indices(row_type);
        } else {
            row_indices = this->beta().indices(row_type);
        }

        // Determine the spin component and dimension of the columns.
        std::vector<size_t> column_indices;
        if (column_spin == Spin::alpha) {
            column_indices = this->alpha().indices(column_type);
        } else {
            column_indices = this->beta().indices(column_type);
        }

        return ImplicitMatrixSlice<Scalar>::FromIndices(row_indices, column_indices, M);
    }


    /**
     *  Create an implicit mathematical object that can serve as the representation of a object with the given occupation types, from the dense tensor representation of a slice.
     *  Each occupation type can belong to either the alpha or beta spin component, initializing a mixed spin mathematical object.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit tensor
     * 
     *  @param axis1_type                   the spinor occupation type for the first tensor axis
     *  @param axis1_spin                   the spin component for the first tensor axis
     *  @param axis2_type                   the spinor occupation type for the second tensor axis
     *  @param axis2_spin                   the spin component for the second tensor axis
     *  @param axis3_type                   the spinor occupation type for the third tensor axis
     *  @param axis3_spin                   the spin component for the third tensor axis
     *  @param axis4_type                   the spinor occupation type for the fourth tensor axis
     *  @param axis4_spin                   the spin component for the fourth tensor axis
     *  @param T                            the dense representation of the tensor slice
     * 
     *  @return an implicit rank-four tensor slice, according to the given occupation types
     */
    template <typename Scalar>
    ImplicitRankFourTensorSlice<Scalar> createMixedRepresentableObjectFor(const OccupationType axis1_type, const Spin axis1_spin, const OccupationType axis2_type, const Spin axis2_spin,
                                                                          const OccupationType axis3_type, const Spin axis3_spin, const OccupationType axis4_type, const Spin axis4_spin, const Tensor<Scalar, 4>& T) const {

        // Prepare the necessary members for ImplicitRankFourTensor.
        // The spin component of each axis determines from which of the two internal OrbitalSpace members the orbitals of the chosen occupation type are used.
        // Determine the spin component and dimension of the first axis.
        std::vector<size_t> axis1_indices;
        if (axis1_spin == Spin::alpha) {
            axis1_indices = this->alpha().indices(axis1_type);
        } else {
            axis1_indices = this->beta().indices(axis1_type);
        }

        // Determine the spin component and dimension of the second axis.
        std::vector<size_t> axis2_indices;
        if (axis2_spin == Spin::alpha) {
            axis2_indices = this->alpha().indices(axis2_type);
        } else {
            axis2_indices = this->beta().indices(axis2_type);
        }

        // Determine the spin component and dimension of the third axis.
        std::vector<size_t> axis3_indices;
        if (axis3_spin == Spin::alpha) {
            axis3_indices = this->alpha().indices(axis3_type);
        } else {
            axis3_indices = this->beta().indices(axis3_type);
        }

        // Determine the spin component and dimension of the fourth axis.
        std::vector<size_t> axis4_indices;
        if (axis4_spin == Spin::alpha) {
            axis4_indices = this->alpha().indices(axis4_type);
        } else {
            axis4_indices = this->beta().indices(axis4_type);
        }

        return ImplicitRankFourTensorSlice<Scalar>::FromIndices(axis1_indices, axis2_indices, axis3_indices, axis4_indices, T);
    }


    /**
     *  Create an implicit, zero-initialized, mathematical object that can serve as the representation of a object with the given occupation types.
     *  Each occupation type can belong to either the alpha or beta spin component, initializing a mixed spin mathematical object.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit matrix
     * 
     *  @param row_type                     the spinor occupation type for the rows
     *  @param row_spin                     the spin component for the rows
     *  @param column_type                  the spinor occupation type for the columns
     *  @param column_spin                  the spin component for the columns
     * 
     *  @return a zero-initialized implicit matrix slice, according to the given occupation types
     */
    template <typename Scalar>
    ImplicitMatrixSlice<Scalar> initializeMixedRepresentableObjectFor(const OccupationType row_type, const Spin row_spin,
                                                                      const OccupationType column_type, const Spin column_spin) const {

        // Prepare the necessary members for the other method.
        size_t rows {};
        if (row_spin == Spin::alpha) {
            rows = this->alpha().numberOfOrbitals(row_type);
        } else {
            rows = this->beta().numberOfOrbitals(row_type);
        }

        size_t columns {};
        if (column_spin == Spin::alpha) {
            columns = this->alpha().numberOfOrbitals(column_type);
        } else {
            columns = this->beta().numberOfOrbitals(column_type);
        }

        const MatrixX<Scalar> M = MatrixX<Scalar>::Zero(rows, columns);

        return this->createMixedRepresentableObjectFor<Scalar>(row_type, column_type, M);
    }


    /**
     *  Create an implicit, zero-initialized, mathematical object that can serve as the representation of a object with the given occupation types, from the dense tensor representation of a slice.
     *  Each occupation type can belong to either the alpha or beta spin component, initializing a mixed spin mathematical object.
     * 
     *  @tparam Scalar                      the scalar type of the elements of the implicit tensor
     * 
     *  @param axis1_type                   the spinor occupation type for the first tensor axis
     *  @param axis1_spin                   the spin component for the first tensor axis
     *  @param axis2_type                   the spinor occupation type for the second tensor axis
     *  @param axis2_spin                   the spin component for the second tensor axis
     *  @param axis3_type                   the spinor occupation type for the third tensor axis
     *  @param axis3_spin                   the spin component for the third tensor axis
     *  @param axis4_type                   the spinor occupation type for the fourth tensor axis
     *  @param axis4_spin                   the spin component for the fourth tensor axis
     * 
     *  @return a zero-initialized implicit mixed spin rank-four tensor slice, according to the given occupation types and spin components. 
     */
    template <typename Scalar>
    ImplicitRankFourTensorSlice<Scalar> initializeMixedRepresentableObjectFor(const OccupationType axis1_type, const Spin axis1_spin, const OccupationType axis2_type, const Spin axis2_spin,
                                                                              const OccupationType axis3_type, const Spin axis3_spin, const OccupationType axis4_type, const Spin axis4_spin) const {

        // Prepare the necessary members for the other method.
        // The spin component of each axis determines from which of the two internal OrbitalSpace members the orbitals of the chosen occupation type are used.
        // We need to use static_cast to determine the Tensor dimension.

        // Determine the spin component and dimension of the first axis.
        size_t axis1_dimension {};
        if (axis1_spin == Spin::alpha) {
            axis1_dimension = static_cast<long>(this->alpha().numberOfOrbitals(axis1_type));
        } else {
            axis1_dimension = static_cast<long>(this->beta().numberOfOrbitals(axis1_type));
        }

        // Determine the spin component and dimension of the second axis.
        size_t axis2_dimension {};
        if (axis2_spin == Spin::alpha) {
            axis2_dimension = static_cast<long>(this->alpha().numberOfOrbitals(axis2_type));
        } else {
            axis2_dimension = static_cast<long>(this->beta().numberOfOrbitals(axis2_type));
        }

        // Determine the spin component and dimension of the second axis.
        size_t axis3_dimension {};
        if (axis3_spin == Spin::alpha) {
            axis3_dimension = static_cast<long>(this->alpha().numberOfOrbitals(axis3_type));
        } else {
            axis3_dimension = static_cast<long>(this->beta().numberOfOrbitals(axis3_type));
        }

        // Determine the spin component and dimension of the second axis.
        size_t axis4_dimension {};
        if (axis4_spin == Spin::alpha) {
            axis4_dimension = static_cast<long>(this->alpha().numberOfOrbitals(axis4_type));
        } else {
            axis4_dimension = static_cast<long>(this->beta().numberOfOrbitals(axis4_type));
        }

        // Zero-initialize the tensor of the specified dimensions.
        Tensor<Scalar, 4> T {axis1_dimension, axis2_dimension, axis3_dimension, axis4_dimension};
        T.setZero();

        return this->createMixedRepresentableObjectFor(axis1_type, axis1_spin, axis2_type, axis2_spin, axis3_type, axis3_spin, axis4_type, axis4_spin, T);
    }
};

}  // namespace GQCP
