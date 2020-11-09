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


#include "Basis/Integrals/IntegralCalculator.hpp"
#include "Basis/MullikenPartitioning.hpp"
#include "Basis/SpinorBasis/SimpleSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/Spinor.hpp"
#include "Basis/Transformations/JacobiRotation.hpp"
#include "Basis/Transformations/RTransformationMatrix.hpp"
#include "Mathematical/Representation/SquareMatrix.hpp"
#include "Operator/FirstQuantized/Operator.hpp"
#include "Operator/SecondQuantized/EvaluatableScalarRSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQTwoElectronOperator.hpp"
#include "Utilities/aliases.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A restricted spin-orbital basis, i.e. a spin-orbital basis where the alpha- and beta-spinors are equal.
 * 
 *  @tparam _ExpansionScalar        The scalar type used to represent an expansion coefficient of the spin-orbitals in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                  The type of shell the underlying scalar basis contains.
 */
template <typename _ExpansionScalar, typename _Shell>
class RSpinOrbitalBasis:
    public SimpleSpinOrbitalBasis<_ExpansionScalar, _Shell, RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type of the base spinor basis.
    using BaseSpinorBasis = SimpleSpinorBasis<_ExpansionScalar, RSpinOrbitalBasis<_ExpansionScalar, _Shell>>;

    // The type of transformation matrix that is naturally related to a GSpinorBasis.
    using TM = RTransformationMatrix<ExpansionScalar>;  // TODO: Rename to TransformationMatrix once the class is gone

    // The type that is used for representing the primitive for a basis function of this spin-orbital basis' underlying AO basis.
    using Primitive = typename Shell::Primitive;

    // The type that is used for representing the underlying basis functions of this spin-orbital basis.
    using BasisFunction = typename Shell::BasisFunction;

    // The type that is used to represent a spatial orbital for this spin-orbital basis.
    using SpatialOrbital = LinearCombination<product_t<ExpansionScalar, typename BasisFunction::CoefficientScalar>, BasisFunction>;


public:
    /*
     *  MARK: Constructors
     */

    // Inherit `SimpleSpinOrbitalBasis`'s constructors.
    using SimpleSpinOrbitalBasis<_ExpansionScalar, _Shell, RSpinOrbitalBasis<_ExpansionScalar, _Shell>>::SimpleSpinOrbitalBasis;


    /*
     *  MARK: General information
     */

    /**
     *  @return The number of different spatial orbitals that are used in this restricted spin-orbital basis.
     */
    size_t numberOfSpatialOrbitals() const { return this->scalar_basis.numberOfBasisFunctions(); }

    /**
     *  @return The number of spinors that are described by this restricted spin-orbital basis.
     */
    size_t numberOfSpinors() const { return 2 * this->numberOfSpatialOrbitals(); /* alpha and beta spinors are equal*/ }

    /**
     *  @return The number of spin-orbitals that are described by this restricted spin-orbital basis.
     */
    size_t numberOfSpinOrbitals() const { return this->numberOfSpinors(); }


    /*
     *  MARK: Quantization of first-quantized operators
     */

    /**
     *  Quantize a first-quantized one-electron operator in this restricted spin-orbital basis.
     * 
     *  @param fq_one_op                            The first-quantized one-electron operator.
     * 
     *  @tparam FQOneElectronOperator               The type of the first-quantized one-electron operator.
     * 
     *  @return The second-quantized operator corresponding to the given first-quantized operator, i.e. expressed in/projected onto this spin-orbital basis.
     */
    template <typename FQOneElectronOperator>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> RSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBasis());  // in AO/scalar basis

        ResultOperator op {one_op_par};           // op for 'operator'
        op.transform(this->coefficientMatrix());  // now in the spatial/spin-orbital basis
        return op;
    }


    /**
     *  @param fq_op        the first-quantized Coulomb operator
     * 
     *  @return the second-quantized operator corresponding to the Coulomb operator
     */
    auto quantize(const CoulombRepulsionOperator& fq_op) const -> RSQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Vectorizer> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = RSQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Vectorizer>;

        const auto one_op_par = IntegralCalculator::calculateLibintIntegrals(fq_op, this->scalarBasis());  // in AO/scalar basis

        ResultOperator op {one_op_par};           // op for 'operator'
        op.transform(this->coefficientMatrix());  // now in spatial/spin-orbital basis
        return op;
    }


    /**
     *  Quantize the (one-electron) electronic density operator.
     * 
     *  @param fq_density_op                    the first-quantized density operator
     * 
     *  @return the second-quantized density operator
     */
    EvaluatableScalarRSQOneElectronOperator<ScalarFunctionProduct<LinearCombination<double, BasisFunction>>> quantize(const ElectronicDensityOperator& fq_density_op) const {

        using Evaluatable = ScalarFunctionProduct<LinearCombination<double, BasisFunction>>;  // the evaluatable type for the density operator
        using ResultOperator = EvaluatableScalarRSQOneElectronOperator<Evaluatable>;

        // There aren't any 'integrals' to be calculated for the density operator: we can just multiply every pair of spatial orbitals.
        const auto phi = this->spatialOrbitals();
        const auto K = this->numberOfSpatialOrbitals();

        SquareMatrix<Evaluatable> rho_par {K};  // the matrix representation ('par' for 'parameters') of the second-quantized (one-electron) density operator
        for (size_t p = 0; p < K; p++) {
            for (size_t q = 0; q < K; q++) {
                rho_par(p, q) = phi[p] * phi[q];
            }
        }

        return ResultOperator(rho_par);
    }


    /*
     *  MARK: Orbitals
     */

    /**
     *  @return the set of spatial orbitals that is associated to this spin-orbital basis
     */
    std::vector<SpatialOrbital> spatialOrbitals() const {

        // The spatial orbitals are a linear combination of the basis functions, where every column of the coefficient matrix describes one expansion of a spatial orbital in terms of the basis functions.
        const auto basis_functions = this->scalar_basis.basisFunctions();
        const auto& C = this->C;


        // For all spatial orbitals, proceed to calculate the contraction between the associated coefficient matrix column and the basis functions.
        std::vector<SpatialOrbital> spatial_orbitals;
        spatial_orbitals.reserve(this->numberOfSpatialOrbitals());
        for (size_t p = 0; p < this->numberOfSpatialOrbitals(); p++) {

            // Calculate the spatial orbitals as a contraction between a column of the coefficient matrix and the basis functions.
            SpatialOrbital spatial_orbital {};
            for (size_t mu = 0; mu < basis_functions.size(); mu++) {
                const auto coefficient = this->C.col(p)(mu);
                const auto& function = basis_functions[mu];
                spatial_orbital.append({coefficient}, {function});
            }

            spatial_orbitals.push_back(spatial_orbital);
        }

        return spatial_orbitals;
    }


    /**
     *  @return the set of spin-orbitals that is associated to this spin-orbital basis
     */
    std::vector<Spinor<ExpansionScalar, BasisFunction>> spinOrbitals() const {

        // The spin-orbitals for a restricted spin-orbital basis can be easily constructed from the spatial orbitals, by assigning a zero component once for the beta component of the spin-orbital and once for the alpha component of the spin-orbital.
        const auto spatial_orbitals = this->spatialOrbitals();

        std::vector<Spinor<ExpansionScalar, BasisFunction>> spin_orbitals;
        spin_orbitals.reserve(this->numberOfSpinors());
        for (const auto& spatial_orbital : spatial_orbitals) {

            // Add the alpha- and beta-spin-orbitals accordingly.
            const Spinor<ExpansionScalar, BasisFunction> alpha_spin_orbital {spatial_orbital, 0};  // the '0' int literal can be converted to a zero LinearCombination
            spin_orbitals.push_back(alpha_spin_orbital);

            const Spinor<ExpansionScalar, BasisFunction> beta_spin_orbital {0, spatial_orbital};
            spin_orbitals.push_back(beta_spin_orbital);
        }

        return spin_orbitals;
    }


    /*
     *  MARK: Mulliken partitioning
     */

    /**
     *  Partition this set of restricted spin-orbitals according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken partitioning.
     * 
     *  @return A `MullikenPartitioning` for the AOs selected by the supplied selector function.
     */
    MullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const BasisFunction&)>& selector) const {

        const auto basis_functions = this->scalarBasis().basisFunctions();

        // Find the indices of those basis functions for which the selector returns true.
        std::vector<size_t> ao_indices;
        for (size_t i = 0; i < this->numberOfSpatialOrbitals(); i++) {
            const auto& basis_function = basis_functions[i];

            if (selector(basis_function)) {
                ao_indices.push_back(i);
            }
        }

        return MullikenPartitioning<ExpansionScalar> {ao_indices, this->coefficientMatrix()};
    }


    /**
     *  Partition this set of restricted spin-orbitals according to the Mulliken partitioning scheme.
     * 
     *  @param selector             A function that returns true for shells that should be included the Mulliken partitioning.
     * 
     *  @return A `MullikenPartitioning` for the AOs selected by the supplied selector function.
     */
    MullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const Shell&)>& selector) const {

        const auto shells = this->scalarBasis().shellSet().asVector();

        // Find the indices of those basis functions for which the shell selector returns true.
        std::vector<size_t> ao_indices;
        size_t bf_index = 0;
        for (size_t shell_index = 0; shell_index < shells.size(); shell_index++) {
            const auto& shell = shells[shell_index];

            // If a shell has to be included, include all indices of the basis functions in it.
            const auto number_of_bf_in_shell = shell.numberOfBasisFunctions();
            if (selector(shell)) {
                for (size_t i = 0; i < number_of_bf_in_shell; i++) {
                    ao_indices.push_back(bf_index);
                    bf_index++;
                }
            } else {
                // Increase the current BF index to accommodate to the next shell.
                bf_index += number_of_bf_in_shell;
            }
        }

        return MullikenPartitioning<ExpansionScalar> {ao_indices, this->coefficientMatrix()};
    }


    /**
     *  @param ao_list     indices of the AOs used for the Mulliken populations
     *
     *  @return the Mulliken operator for a set of given AO indices
     *
     *  @note this method is only available for real matrix representations
     */
    template <typename S = ExpansionScalar, typename = IsReal<S>>
    ScalarRSQOneElectronOperator<double> calculateMullikenOperator(const std::vector<size_t>& ao_list) const {

        const auto K = this->numberOfSpatialOrbitals();
        if (ao_list.size() > K) {
            throw std::invalid_argument("RSpinOrbitalBasis::calculateMullikenOperator(std::vector<size_t>): Too many AOs are selected in the given ao_list");
        }

        const SquareMatrix<double> p_a = SquareMatrix<double>::PartitionMatrix(ao_list, K);                       // the partitioning matrix
        const auto S_AO = IntegralCalculator::calculateLibintIntegrals(Operator::Overlap(), this->scalar_basis);  // the overlap matrix expressed in the AO basis

        ScalarRSQOneElectronOperator<double> mulliken_op {0.5 * (this->C.adjoint() * p_a * S_AO * this->C + this->C.adjoint() * S_AO * p_a * this->C)};
        return mulliken_op;
    }
};


/*
 *  MARK: Convenience aliases
 */
template <typename ExpansionScalar, typename Shell>
using RSpinorBasis = RSpinOrbitalBasis<ExpansionScalar, Shell>;


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename _ExpansionScalar, typename _Shell>
struct SpinorBasisTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of transformation matrix that is naturally related to an RSpinOrbitalBasis.
    using TM = RTransformationMatrix<ExpansionScalar>;  // TODO: Rename to TransformationMatrix once the class is gone

    // The second-quantized representation of the overlap operator related to an RSpinOrbitalBasis.
    using SQOverlapOperator = ScalarRSQOneElectronOperator<ExpansionScalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of transformation matrix that is naturally related to an `RSpinOrbitalBasis`.
    using TM = RTransformationMatrix<_ExpansionScalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<RSpinOrbitalBasis<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation that is naturally related to an `RSpinOrbitalBasis`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
