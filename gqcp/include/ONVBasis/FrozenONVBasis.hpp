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
#include "ONVBasis/ONVManipulator.hpp"
#include "ONVBasis/ONVBasis.hpp"


namespace GQCP {


/**
 *  A class that represents a frozen ONV basis: this is a subspace of the N-electron ONV basis in which the first X orbitals are always occupied
 */
class FrozenONVBasis: public BaseFrozenCoreONVBasis, public ONVManipulator<FrozenONVBasis> {
protected:
    size_t X;  // number of frozen orbitals/electrons
    ONVBasis active_fock_space;  // active (non-frozen) ONV basis containing only the active electrons (N-X) and orbitals (K-X)

public:
    // CONSTRUCTORS
    /**
     *  @param K        the total number of orbitals
     *  @param N        the total number of electrons
     *  @param X        the number of frozen orbitals and electrons
     */
    FrozenONVBasis(size_t K, size_t N, size_t X);

    /**
     *  @param fock_space       (to be frozen) full ONV basis
     *  @param X                the number of frozen orbitals and electrons
     */
    FrozenONVBasis(const ONVBasis& fock_space, size_t X);


    // DESTRUCTOR
    ~FrozenONVBasis() override = default;


    // GETTERS
    size_t get_number_of_frozen_orbitals() const { return this->X; }
    const ONVBasis& get_active_fock_space() const { return this->active_fock_space; }
    ONVBasisType get_type() const override { return ONVBasisType::FrozenONVBasis; }


    // OVERRIDEN PUBLIC METHODS
    /**
     *  @param representation       a representation of an ONV
     *
     *  @return the next bitstring permutation in the frozen ONV basis
     *
     *      Examples (first orbital is frozen):
     *          0101 -> 1001
     *         01101 -> 10011
     */
    size_t ulongNextPermutation(size_t representation) const override;

    /**
     *  @param representation      a representation of an ONV
     *
     *  @return the address (i.e. the ordering number) of the given ONV
     */
    size_t getAddress(size_t representation) const override;

    /**
      *  Calculate unsigned representation for a given address
      *
      *  @param address                 the address of the representation is calculated
      *
      *  @return unsigned representation of the address
      */
    size_t calculateRepresentation(size_t address) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a one electron operator
     */
    size_t countOneElectronCouplings(const ONV& onv) const override;

    /**
     *  @param onv       the ONV
     *
     *  @return the amount of ONVs (with a larger address) this ONV would couple with given a two electron operator
     */
    size_t countTwoElectronCouplings(const ONV& onv) const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a one electron coupling scheme in the ONV basis
     */
    size_t countTotalOneElectronCouplings() const override;

    /**
     *  @return the amount non-zero (non-diagonal) couplings of a two electron coupling scheme in the ONV basis
     */
    size_t countTotalTwoElectronCouplings() const override;


    // PUBLIC METHODS
    /**
     *  If we have
     *      FrozenONVBasis fock_space;
     *
     *  This makes sure that we can call
     *      fock_space.getAddress(onv);
     *  instead of the syntax
     *      fock_space.ONVManipulator<FrozenONVBasis>::getAddress(onv);
     */
    using ONVManipulator<FrozenONVBasis>::getAddress;
};


}  // namespace GQCP
