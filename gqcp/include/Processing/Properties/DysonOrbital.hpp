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

#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A Dyson orbital. It encapsulates a set of Dyson amplitudes that are its expansion coefficients in a spinor basis.
 * 
 *  @tparam Scalar         The scalar type of the Dyson amplitudes.
 */
template <typename _Scalar>
class DysonOrbital {
public:
    // The scalar type used for a single parameter: real or complex.
    using Scalar = _Scalar;


private:
    // The Dyson amplitudes for the Dyson orbital. They are indicated as <N-1|a_p|N>, where 'p' is the index of spinor p.
    VectorX<Scalar> m_amplitudes;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  The default constructor.
     */
    DysonOrbital() = default;

    /**
     * Construct a Dyson orbital with all Dyson amplitudes set to zero.
     * 
     * @param K     The dimension of the amplitude vector, i.e. the dimension of a related orbital basis.
     */
    DysonOrbital(const size_t K) :
        m_amplitudes {VectorX<Scalar>::Zero(K)} {}


    /**
     * Construct a Dyson orbital with corresponding Dyson amplitudes.
     * 
     * @param amplitudes     The Dyson amplitudes for the Dyson orbital. They are indicated as <N-1|a_p|N>, where 'p' is the index of spinor p.
     */
    DysonOrbital(const VectorX<Scalar>& amplitudes) :
        m_amplitudes {amplitudes} {}

    /*
     *  MARK: Named constructors
     */

    /**
     *  Create a Dyson orbital from the formula for its amplitudes `<N_1|a_p|N>`.
     * 
     *  @param linear_expansion_J        The N-electron wave function in a spin-resolved ONV basis.
     *  @param linear_expansion_I        The (N-1)-electron wave function in a spin-resolved ONV basis. It should be expressed in the same orbital basis as the N-electron wave function.
     *
     *  @return A Dyson orbital incorporating Dyson amplitudes.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, DysonOrbital<double>> TransitionAmplitudes(const LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion_J, const LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion_I) {

        const auto onv_basis_J = linear_expansion_J.onvBasis();
        const auto onv_basis_I = linear_expansion_I.onvBasis();

        if ((onv_basis_J.alpha().numberOfElectrons() - onv_basis_I.alpha().numberOfElectrons() != 0) && (onv_basis_J.beta().numberOfElectrons() - onv_basis_I.beta().numberOfElectrons() != 1)) {
            if ((onv_basis_J.alpha().numberOfElectrons() - onv_basis_I.alpha().numberOfElectrons() != 1) && (onv_basis_J.beta().numberOfElectrons() - onv_basis_I.beta().numberOfElectrons() != 0)) {
                throw std::runtime_error("DysonOrbital::TransitionAmplitudes(LinearExpansion, LinearExpansion): linear_expansion_I is not expressed in a spin-resolved ONV basis with one fewer electron than linear_expansion_J.");
            }
        }

        // Initialize environment variables.

        // The 'passive' ONV basis is the ONV basis that is equal for both wave functions.
        // The 'target' ONV basis has an electron difference of one.
        // We initialize the variables for the case in which they differ in one beta electron, if this isn't the case, we will update it later.
        auto passive_onv_basis_J = onv_basis_J.alpha();
        auto passive_onv_basis_I = onv_basis_I.alpha();
        auto target_onv_basis_J = onv_basis_J.beta();
        auto target_onv_basis_I = onv_basis_I.beta();

        // Mod variables relate to the modification of the address jump in coefficient index according to the ordering of the spin ONVs.
        size_t passive_mod_J = target_onv_basis_J.dimension();
        size_t passive_mod_I = target_onv_basis_I.dimension();
        size_t target_mod = 1;

        // If instead the ONV bases differ by one alpha electron we re-assign the variables to match the algorithm.
        if ((onv_basis_J.alpha().numberOfElectrons() - onv_basis_I.alpha().numberOfElectrons() == 1) && (onv_basis_J.beta().numberOfElectrons() - onv_basis_I.beta().numberOfElectrons() == 0)) {
            passive_onv_basis_J = target_onv_basis_J;
            passive_onv_basis_I = target_onv_basis_I;
            target_onv_basis_J = onv_basis_J.alpha();
            target_onv_basis_I = onv_basis_I.alpha();

            passive_mod_J = 1;
            passive_mod_I = 1;
            target_mod = passive_onv_basis_J.dimension();
        }

        const auto& ci_coeffs_J = linear_expansion_J.coefficients();
        const auto& ci_coeffs_I = linear_expansion_I.coefficients();

        VectorX<double> dyson_coeffs = VectorX<double>::Zero(onv_basis_J.alpha().numberOfOrbitals());

        // The actual algorithm to determine the Dyson amplitudes.

        // Since we want to calculate the overlap between two wave functions, the ONVs should have an equal number of electrons.
        // We therefore iterate over the ONVs of the 'target' ONV basis, which all have an electron more, and annihilate in one of the orbitals (let the index of that orbital be called 'p').
        // By calculating the overlap in the (N-1)-electron ONV basis, we can calculate the contributions to the  'p'-th coefficient (i.e. the Dyson amplitude) of the Dyson orbital.
        SpinUnresolvedONV onv = target_onv_basis_J.constructONVFromAddress(0);

        for (size_t Jt = 0; Jt < target_onv_basis_J.dimension(); Jt++) {           // Jt loops over addresses of the target ONV basis.
            int sign = -1;                                                         // Total phase factor of all the annihilations that have occurred.
            for (size_t e = 0; e < target_onv_basis_J.numberOfElectrons(); e++) {  // Loop over electrons in the ONV.

                // Annihilate on the corresponding orbital, to make sure we can calculate overlaps in the (N-1)-'target' ONV basis.
                sign *= -1;
                size_t p = onv.occupationIndexOf(e);
                onv.annihilate(p);

                // Now, we calculate the overlap in the (N-1)-electron 'target' ONV basis.
                // In order to access the correct coefficients for the, we need the address of the resulting (annihilated) ONV inside the 'target' (N-1)-electron ONV basis.
                size_t address = target_onv_basis_I.addressOf(onv.unsignedRepresentation());

                double coeff = 0;
                for (size_t Jp = 0; Jp < passive_onv_basis_J.dimension(); Jp++) {                                                                // Jp loops over the addresses of the passive ONV basis.
                    coeff += sign * ci_coeffs_J(Jt * target_mod + Jp * passive_mod_J) * ci_coeffs_I(address * target_mod + Jp * passive_mod_I);  // Access the indices of the coefficient vectors.
                }
                dyson_coeffs(p) += coeff;
                onv.create(p);  // Allow the iteration to continue with the original ONV.
            }

            if (Jt < target_onv_basis_J.dimension() - 1) {  // Prevent the last permutation from occurring.
                target_onv_basis_J.transformONVToNextPermutation(onv);
            }
        }  // Target address (Jt) loop.

        return DysonOrbital<Scalar>(dyson_coeffs);
    }

    /**
     *  Create a Dyson orbital from the formula for its amplitudes `<N_1|a_p|N>`.
     * 
     *  @param linear_expansion_J        The N-electron wave function in a spin-unresolved ONV basis.
     *  @param linear_expansion_I        The (N-1)-electron wave function in a spin-unresolved ONV basis. It should be expressed in the same orbital basis as the N-electron wave function.
     *
     *  @return A Dyson orbital incorporating Dyson amplitudes.
     */
    template <typename Z = Scalar>
    static enable_if_t<std::is_same<Z, double>::value, DysonOrbital<double>> TransitionAmplitudes(const LinearExpansion<double, SpinUnresolvedONVBasis>& linear_expansion_J, const LinearExpansion<double, SpinUnresolvedONVBasis>& linear_expansion_I) {

        // Initialize environment variables.

        const auto onv_basis_J = linear_expansion_J.onvBasis();
        const auto onv_basis_I = linear_expansion_I.onvBasis();

        if (onv_basis_J.numberOfElectrons() - onv_basis_I.numberOfElectrons() != 1) {
            throw std::runtime_error("DysonOrbital::TransitionAmplitudes(LinearExpansion, LinearExpansion): linear_expansion_I is not expressed in a spin-unresolved ONV basis with one fewer electron than linear_expansion_J.");
        }

        const auto& ci_coeffs_J = linear_expansion_J.coefficients();
        const auto& ci_coeffs_I = linear_expansion_I.coefficients();

        VectorX<double> dyson_coeffs = VectorX<double>::Zero(onv_basis_J.numberOfOrbitals());

        // The actual algorithm to determine the Dyson amplitudes.

        // Since we want to calculate the overlap between two wave functions, the ONVs should have an equal number of electrons.
        // We therefore iterate over the ONVs of the N-electron ONV basis, which all have an electron more, and annihilate in one of the orbitals (let the index of that orbital be called 'p').
        // By calculating the overlap in the (N-1)-electron ONV basis, we can calculate the contributions to the  'p'-th coefficient (i.e. the Dyson amplitude) of the Dyson orbital.
        SpinUnresolvedONV onv = onv_basis_J.constructONVFromAddress(0);
        double coeff = 0;

        for (size_t J = 0; J < onv_basis_J.dimension(); J++) {              // J loops over addresses of the N-electron spin-unresolved ONV basis.
            int sign = -1;                                                  // Total phase factor of all the annihilations that have occurred.
            for (size_t e = 0; e < onv_basis_J.numberOfElectrons(); e++) {  // Loop over electrons in the ONV.

                // Annihilate on the corresponding orbital, to make sure we can calculate overlaps in the (N-1)-electron ONV basis.
                sign *= -1;
                size_t p = onv.occupationIndexOf(e);
                onv.annihilate(p);

                // Now, we calculate the overlap in the (N-1)-electron ONV basis.
                // In order to access the correct coefficients for the, we need the address of the resulting (annihilated) ONV inside the (N-1)-electron ONV basis.
                size_t address = onv_basis_I.addressOf(onv.unsignedRepresentation());

                coeff = sign * ci_coeffs_J(J) * ci_coeffs_I(address);  // Access the indices of the coefficient vectors.

                dyson_coeffs(p) += coeff;
                onv.create(p);  // Allow the iteration to continue with the original ONV.
            }

            if (J < onv_basis_J.dimension() - 1) {  // Prevent the last permutation from occurring.
                onv_basis_J.transformONVToNextPermutation(onv);
            }
        }  // Loop over J ONVs.

        return DysonOrbital<Scalar>(dyson_coeffs);
    }

    /*
     *  MARK: Access
     */

    /**
     *  @return The Dyson amplitudes for the Dyson orbital. They are indicated as <N-1|a_p|N>, where 'p' is the index of spinor p.
     */
    const VectorX<Scalar>& amplitudes() const { return this->m_amplitudes; }
};

}  // namespace GQCP
