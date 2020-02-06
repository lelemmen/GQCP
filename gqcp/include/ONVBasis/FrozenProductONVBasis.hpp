// This file is part of GQCG-gqcp.
// 
// Copyright (C) 2017-2019  the GQCG developers
// 
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
// 
#pragma once


#include "ONVBasis/BaseONVBasis.hpp"
#include "ONVBasis/BaseFrozenCoreONVBasis.hpp"
#include "ONVBasis/FrozenONVBasis.hpp"
#include "ONVBasis/ProductONVBasis.hpp"
#include "Operator/SecondQuantized/USQHamiltonian.hpp"


namespace GQCP {


/**
 *  A class that represents the product of two frozen ONV basiss (alpha and beta).
 */
class FrozenProductONVBasis: public BaseFrozenCoreONVBasis {
private:
    size_t X;  // number of frozen orbitals/electrons

    FrozenONVBasis frozen_fock_space_alpha;
    FrozenONVBasis frozen_fock_space_beta;

    ProductONVBasis active_product_fock_space;  // active (non-frozen) product ONV basis containing only the active electrons (N_alpha-X, N_beta-X) and orbitals (K-X)

public:
    // CONSTRUCTORS
    /**
     *  @param K            the total number of orbitals (equal for alpha and beta)
     *  @param N_alpha      the total number of alpha electrons
     *  @param N_beta       the total number of beta electrons
     *  @param X            the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    FrozenProductONVBasis(size_t K, size_t N_alpha, size_t N_beta, size_t X);

    /**
     *  @param fock_space       (to be frozen) full product ONV basis
     *  @param X                the number of frozen orbitals and electrons (equal for alpha and beta)
     */
    FrozenProductONVBasis(const ProductONVBasis& fock_space, size_t X);


    // GETTERS
    size_t get_N_alpha() const { return this->frozen_fock_space_alpha.get_N(); }
    size_t get_N_beta() const { return this->frozen_fock_space_beta.get_N(); }
    size_t get_number_of_frozen_orbitals() const { return this->X;}

    const FrozenONVBasis& get_frozen_fock_space_alpha() const { return this->frozen_fock_space_alpha; }
    const FrozenONVBasis& get_frozen_fock_space_beta() const { return this->frozen_fock_space_beta; }

    const ProductONVBasis& get_active_product_fock_space() const { return this->active_product_fock_space; }
    ONVBasisType get_type() const override { return ONVBasisType::FrozenProductONVBasis; }

    using BaseFrozenCoreONVBasis::evaluateOperatorDense;
    using BaseFrozenCoreONVBasis::evaluateOperatorDiagonal;


    // UNRESTRICTED
    /**
     *  Evaluate the Hamiltonian in a dense matrix
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param diagonal_values                bool to indicate if diagonal values will be calculated
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    SquareMatrix<double> evaluateOperatorDense(const USQHamiltonian<double>& usq_hamiltonian, bool diagonal_values) const;

    /**
     *  Evaluate the diagonal of the Hamiltonian
     *
     *  @param usq_hamiltonian          the Hamiltonian expressed in an unrestricted orthonormal basis 
     *
     *  @return the Hamiltonian's diagonal evaluation in a vector with the dimension of the ONV basis
     */
    VectorX<double> evaluateOperatorDiagonal(const USQHamiltonian<double>& usq_hamiltonian) const;

    /**
     *  Evaluate the unrestricted Hamiltonian in a matrix vector product
     *
     *  @param usq_hamiltonian                the Hamiltonian expressed in an unrestricted orthonormal basis 
     *  @param x                              the vector upon which the evaluation acts 
     *  @param diagonal                       the diagonal evaluated in the ONV basis
     *
     *  @return the Hamiltonian's evaluation in a dense matrix with the dimensions of the ONV basis
     */
    VectorX<double> evaluateOperatorMatrixVectorProduct(const USQHamiltonian<double>& usq_hamiltonian, const VectorX<double>& x, const VectorX<double>& diagonal) const;
};


}  // namespace GQCP
