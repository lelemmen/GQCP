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


#include "Basis/MullikenPartitioning/GMullikenPartitioning.hpp"
#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/ScalarBasis/LondonGTOShell.hpp"
#include "Basis/ScalarBasis/ScalarBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/SimpleSpinorBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/GTransformation.hpp"
#include "Molecule/Molecule.hpp"
#include "Operator/FirstQuantized/CoulombRepulsionOperator.hpp"
#include "Operator/FirstQuantized/DiamagneticOperator.hpp"
#include "Operator/FirstQuantized/ElectronicDipoleOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpinSquaredOperator.hpp"
#include "Operator/FirstQuantized/ElectronicSpin_zOperator.hpp"
#include "Operator/FirstQuantized/FQMolecularHamiltonian.hpp"
#include "Operator/FirstQuantized/FQMolecularPauliHamiltonian.hpp"
#include "Operator/FirstQuantized/KineticOperator.hpp"
#include "Operator/FirstQuantized/NuclearAttractionOperator.hpp"
#include "Operator/FirstQuantized/OrbitalZeemanOperator.hpp"
#include "Operator/FirstQuantized/OverlapOperator.hpp"
#include "Operator/FirstQuantized/SpinZeemanOperator.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/GSQTwoElectronOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "Utilities/type_traits.hpp"


namespace GQCP {


/**
 *  A general spinor basis, i.e. a spinor basis without any restrictions on the expansion of the alpha and beta components of the spinors in terms of the underlying (possibly different) scalar bases.
 *
 *  @tparam _ExpansionScalar        The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
 *  @tparam _Shell                  The type of shell the underlying scalar bases contain.
 *
 *  @note The individual columns of the coefficient matrix represent the spinors of this basis; they are not ordered by increasing single-particle energy.
 */
template <typename _ExpansionScalar, typename _Shell>
class GSpinorBasis:
    public SimpleSpinorBasis<_ExpansionScalar, GSpinorBasis<_ExpansionScalar, _Shell>> {
public:
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of shell the underlying scalar bases contain.
    using Shell = _Shell;

    // The type of the base spinor basis.
    using Base = SimpleSpinorBasis<_ExpansionScalar, GSpinorBasis<_ExpansionScalar, _Shell>>;

    // The type of transformation that is naturally related to a `GSpinorBasis`.
    using Transformation = GTransformation<ExpansionScalar>;

    // The type that is used for representing the primitive for a basis function of this spin-orbital basis' underlying AO basis.
    using Primitive = typename Shell::Primitive;

    // The type that is used for representing the underlying basis functions of this spin-orbital basis.
    using BasisFunction = typename Shell::BasisFunction;


private:
    // The scalar bases for the alpha and beta components of the spinors.
    SpinResolved<ScalarBasis<Shell>> scalar_bases;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  Create a `GSpinorBasis` from scalar bases that are not necessarily equal.
     *
     *  @param scalar_bases                 The scalar bases for the alpha and beta components of the spinors.
     *  @param C                            The transformation that relates the current set of spinors with the atomic spinors.
     */
    GSpinorBasis(const SpinResolved<ScalarBasis<Shell>>& scalar_bases, const Transformation& C) :
        Base(C),
        scalar_bases {scalar_bases} {

        // Check if the dimensions of the given objects are compatible
        const auto K_alpha = scalar_bases.alpha().numberOfBasisFunctions();
        const auto K_beta = scalar_bases.beta().numberOfBasisFunctions();

        if (C.numberOfOrbitals() != K_alpha + K_beta) {
            throw std::invalid_argument("GSpinorBasis(const SpinResolved<ScalarBasis<Shell>>&, const Transformation&): The given dimensions of the scalar bases and coefficient matrix are incompatible.");
        }
    }

    /**
     *  Create a `GSpinorBasis` from scalar bases that are not necessarily equal.
     *
     *  @param alpha_scalar_basis           The scalar basis in which the alpha components of the spinors are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta components of the spinors are expanded.
     *  @param C                            The transformation that relates the current set of spinors with the atomic spinors.
     */
    GSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis, const Transformation& C) :
        GSpinorBasis(SpinResolved<ScalarBasis<Shell>> {alpha_scalar_basis, beta_scalar_basis}, C) {}


    /**
     *  Construct a generalized spinor basis in which both underlying scalar bases are equal.
     * ^
     *  @param scalar_basis             The scalar basis in which both the alpha and beta components are expanded.
     *  @param C                        The transformation that relates the current set of spinors with the atomic spinors.
     */
    GSpinorBasis(const ScalarBasis<Shell>& scalar_basis, const Transformation& C) :
        GSpinorBasis(scalar_basis, scalar_basis, C) {}


    /**
     *  Construct a generalized spinor basis with two different underlying scalar basis, and a coefficient matrix being the identity. The resulting spinor basis corresponds to the atomic spinors (AOs).
     *
     *  @param alpha_scalar_basis           The scalar basis in which the alpha components of the spinors are expanded.
     *  @param beta_scalar_basis            The scalar basis in which the beta components of the spinors are expanded.
     */
    GSpinorBasis(const ScalarBasis<Shell>& alpha_scalar_basis, const ScalarBasis<Shell>& beta_scalar_basis) :
        GSpinorBasis(alpha_scalar_basis, beta_scalar_basis,
                     Transformation::Identity(alpha_scalar_basis.numberOfBasisFunctions(), beta_scalar_basis.numberOfBasisFunctions())) {}


    /**
     *  Construct a generalized spinor basis in which both underlying scalar bases are equal, and a coefficient matrix being the identity. The resulting spinor basis corresponds to the atomic spinors (AOs).
     *
     *  @param scalar_basis             The scalar basis in which both the alpha and beta components are expanded.
     */
    GSpinorBasis(const ScalarBasis<Shell>& scalar_basis) :
        GSpinorBasis(scalar_basis, scalar_basis) {}


    /**
     *  Construct a generalized spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the nuclear framework. The resulting spinor basis corresponds to the non-orthogonal atomic spinors (AOs).
     *
     *  @param nuclear_framework        The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    GSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name) :
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name)) {}


    /**
     *  Construct a generalized spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spinor basis corresponds to the non-orthogonal atomic spinors (AOs).
     *
     *  @param molecule                 The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name)) {}


    /**
     *  Construct a generalized spinor basis with an underlying scalar basis (equal for both the alpha and beta components) that is made by placing shells corresponding to the basisset specification on every nucleus of the molecule. The resulting spinor basis corresponds to the non-orthogonal London atomic spinors (AOs).
     *
     *  @param molecule                 The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name            The name of the basisset, e.g. "STO-3G".
     *  @param B                        The homogeneous magnetic field.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    template <typename Z = Shell>
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name, const HomogeneousMagneticField& B,
                 typename std::enable_if<std::is_same<Z, LondonGTOShell>::value>::type* = 0) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name, B)) {}


    /**
     *  Construct a generalized spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the nuclear framework. The resulting spinor basis corresponds to the non-orthogonal atomic spinors (AOs).
     *
     *  @param nuclear_framework            The nuclear framework containing the nuclei on which the shells should be centered.
     *  @param basisset_name_alpha          The name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component.
     *  @param basisset_name_beta           The name of the basisset, e.g. "STO-3G", used for the expansion of the beta component.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    GSpinorBasis(const NuclearFramework& nuclear_framework, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        GSpinorBasis(ScalarBasis<Shell>(nuclear_framework, basisset_name_alpha),
                     ScalarBasis<Shell>(nuclear_framework, basisset_name_beta)) {}


    /**
     *  Construct a generalized spinor basis with a underlying scalar bases made by placing shells corresponding to the basisset specifications on every nucleus of the molecule. The resulting spinor basis corresponds to the non-orthogonal atomic spinors (AOs).
     *
     *  @param molecule                     The molecule containing the nuclei on which the shells should be centered.
     *  @param basisset_name_alpha          The name of the basisset, e.g. "STO-3G", used for the expansion of the alpha component.
     *  @param basisset_name_beta           The name of the basisset, e.g. "STO-3G", used for the expansion of the beta component.
     *
     *  @note The normalization factors of the spherical (or axis-aligned Cartesian) GTO primitives are embedded in the contraction coefficients of the underlying shells.
     */
    GSpinorBasis(const Molecule& molecule, const std::string& basisset_name_alpha, const std::string& basisset_name_beta) :
        GSpinorBasis(ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_alpha),
                     ScalarBasis<Shell>(molecule.nuclearFramework(), basisset_name_beta)) {}


    /*
     *  MARK: Named constructors
     */

    /**
     *  Convert a restricted spin-orbital basis into a generalized framework, yielding a generalized coefficient matrix that is spin-blocked out.
     *
     *  @param r_spinor_basis           The restricted spinor basis.
     *
     *  @return The restricted spinor basis as a generalized one.
     */
    static GSpinorBasis<ExpansionScalar, Shell> FromRestricted(const RSpinOrbitalBasis<ExpansionScalar, Shell>& r_spinor_basis) {

        // Create an USpinOrbitalBasis from the restricted one and use ::FromUnrestricted.
        const auto u_spinor_basis = USpinOrbitalBasis<ExpansionScalar, Shell>::FromRestricted(r_spinor_basis);

        return GSpinorBasis<ExpansionScalar, Shell>::FromUnrestricted(u_spinor_basis);
    }


    /**
     *  Convert an unrestricted spin-orbital basis into a generalized framework, yielding a generalized coefficient matrix that is spin-blocked out.
     *
     *  @param u_spinor_basis           The unrestricted spinor basis.
     *
     *  @return The generalized spinor basis corresponding to the unrestricted spin-orbital basis.
     *
     *  @note We assume that the unrestricted spin-orbital basis has equal underlying scalar bases for the alpha- and beta-spin-orbitals.
     */
    static GSpinorBasis<ExpansionScalar, Shell> FromUnrestricted(const USpinOrbitalBasis<ExpansionScalar, Shell>& u_spinor_basis) {

        const auto C_general = GTransformation<ExpansionScalar>::FromUnrestricted(u_spinor_basis.expansion());

        return GSpinorBasis<ExpansionScalar, Shell>(u_spinor_basis.alpha().scalarBasis(), C_general);  // Assume the alpha- and beta- scalar bases are equal.
    }


    /*
     *  MARK: Coefficients
     */

    /**
     *  @param sigma        Alpha or beta.
     *
    *  @return The number of coefficients that are used for the expansion of the requested spin-component of a spinor.
     */
    size_t numberOfCoefficients(const Spin& sigma) const { return this->scalarBases().component(sigma).numberOfBasisFunctions(); }


    /*
     *  MARK: General info
     */

    /**
     *  @return The number of spinors that are described by this generalized spinor basis.
     */
    size_t numberOfSpinors() const {

        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);

        return K_alpha + K_beta;
    }


    /*
     *  MARK: Quantization of first-quantized operators (GTOShell)
     */

    /**
     *  Quantize a scalar spin-independent one-electron operator in this general spinor basis. Spin-independent one-electron operators are those whose two-component matrix operator form contains the same scalar operator in the top-left and bottom-right corner.
     *
     *  @param fq_one_op            A spin-independent first-quantized operator.
     *
     *  @return The second-quantized representation of the given operator.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, GSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;


        // The strategy for calculating the matrix representation of the one-electron operator in this spinor basis is to:
        //  1. Express the operator in the underlying scalar bases; and
        //  2. Afterwards transform them using the current coefficient matrix.
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();
        SquareMatrix<ResultScalar> f = SquareMatrix<ResultScalar>::Zero(M);  // the total result

        // 1. Express the operator in the underlying scalar bases: spin-independent operators only have alpha-alpha and beta-beta blocks.
        const auto F_alpha = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBases().alpha());
        const auto F_beta = IntegralCalculator::calculateLibintIntegrals(fq_one_op, this->scalarBases().beta());

        f.topLeftCorner(K_alpha, K_alpha) = F_alpha;
        f.bottomRightCorner(K_beta, K_beta) = F_beta;

        // 2. Transform using the current coefficient matrix.
        ResultOperator op {{f}};  // op for 'operator'
        op.transform(this->expansion());
        return op;
    }


    /**
     *  Quantize the electronic spin operator in this general spinor basis.
     *
     *  @param fq_one_op        The (first-quantized) electronic spin operator.
     *
     *  @return The electronic spin operator expressed in this spinor basis.
     */
    template <typename Z = Shell>
    auto quantize(const ElectronicSpinOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, GSQOneElectronOperator<product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>, ElectronicSpinOperator::Vectorizer>> {

        using ResultScalar = product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQOneElectronOperator<ResultScalar, ElectronicSpinOperator::Vectorizer>;

        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();

        // The strategy to quantize the spin operator is as follows.
        //  1. First, calculate the necessary overlap integrals over the scalar bases.
        //  2. Then, construct the scalar basis representations of the components of the spin operator by placing the overlaps into the correct blocks.
        //  3. Transform the components (in scalar basis) with the current coefficient matrix to yield the components in spinor basis.

        SquareMatrix<ResultScalar> S_x = SquareMatrix<ResultScalar>::Zero(M);
        SquareMatrix<ResultScalar> S_y = SquareMatrix<ResultScalar>::Zero(M);
        SquareMatrix<ResultScalar> S_z = SquareMatrix<ResultScalar>::Zero(M);


        // 1. Calculate the necessary overlap integrals over the scalar bases.
        const auto S_aa = IntegralCalculator::calculateLibintIntegrals(OverlapOperator(), this->scalarBases().alpha());
        const auto S_ab = IntegralCalculator::calculateLibintIntegrals(OverlapOperator(), this->scalarBases().alpha(), this->scalarBases().beta());
        const auto S_ba = IntegralCalculator::calculateLibintIntegrals(OverlapOperator(), this->scalarBases().beta(), this->scalarBases().alpha());
        const auto S_bb = IntegralCalculator::calculateLibintIntegrals(OverlapOperator(), this->scalarBases().beta());


        // 2. Place the overlaps into the correct blocks.
        S_x.block(0, K_alpha, K_alpha, K_beta) = 0.5 * S_ab;
        S_x.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * S_ba;

        using namespace GQCP::literals;
        S_y.block(0, K_alpha, K_alpha, K_beta) = -0.5 * 1_ii * S_ab;
        S_y.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * 1_ii * S_ba;

        S_z.topLeftCorner(K_alpha, K_alpha) = 0.5 * S_aa;
        S_z.bottomRightCorner(K_beta, K_beta) = -0.5 * S_bb;


        // 3. Transform using the coefficient matrix
        ResultOperator spin_op {std::vector<SquareMatrix<ResultScalar>> {S_x, S_y, S_z}};  // 'op' for operator
        spin_op.transform(this->expansion());
        return spin_op;
    }


    /**
     *  Quantize the Coulomb operator in this general spinor basis.
     *
     *  @param coulomb_op           The first-quantized Coulomb operator.
     *
     *  @return The second-quantized operator corresponding to the Coulomb operator.
     */
    template <typename Z = Shell>
    auto quantize(const CoulombRepulsionOperator& coulomb_op) const -> enable_if_t<std::is_same<Z, GTOShell>::value, GSQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Vectorizer>> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Vectorizer>;

        // The strategy for calculating the matrix representation of the two-electron operator in this spinor basis is to:
        //  1. Calculate the Coulomb integrals in the underlying scalar bases;
        //  2. Place the calculated integrals as 'blocks' in the larger representation, so that we can;
        //  3. Transform the operator using the current coefficient matrix.

        // 1. Calculate the Coulomb integrals in the underlying scalar bases.
        const auto g_aaaa = IntegralCalculator::calculateLibintIntegrals(CoulombRepulsionOperator(), this->scalarBases().alpha());
        const auto g_aabb = IntegralCalculator::calculateLibintIntegrals(CoulombRepulsionOperator(), this->scalarBases().alpha(), this->scalarBases().beta());
        const auto g_bbaa = IntegralCalculator::calculateLibintIntegrals(CoulombRepulsionOperator(), this->scalarBases().beta(), this->scalarBases().alpha());
        const auto g_bbbb = IntegralCalculator::calculateLibintIntegrals(CoulombRepulsionOperator(), this->scalarBases().beta());


        // 2. Place the calculated integrals as 'blocks' in the larger representation
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);

        const auto M = this->numberOfSpinors();
        auto g_par = SquareRankFourTensor<ResultScalar>::Zero(M);  // 'par' for 'parameters'

        // Primed indices are indices in the larger representation, normal ones are those in the smaller tensors.
        for (size_t mu_ = 0; mu_ < M; mu_++) {  // mu 'prime'
            const size_t mu = mu_ % K_alpha;

            for (size_t nu_ = 0; nu_ < M; nu_++) {  // nu 'prime'
                const size_t nu = nu_ % K_alpha;

                for (size_t rho_ = 0; rho_ < M; rho_++) {  // rho 'prime'
                    const size_t rho = rho_ % K_alpha;

                    for (size_t lambda_ = 0; lambda_ < M; lambda_++) {  // lambda 'prime'
                        const size_t lambda = lambda_ % K_alpha;

                        if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_aaaa(mu, nu, rho, lambda);
                        } else if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_aabb(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_bbaa(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g_bbbb(mu, nu, rho, lambda);
                        }
                    }
                }
            }
        }


        // 3. Transform the operator using the current coefficient matrix.
        ResultOperator g_op {g_par};  // 'op' for 'operator'
        g_op.transform(this->expansion());
        return g_op;
    }


    /*
     *  MARK: Quantization of first-quantized operators (LondonGTOShell)
     */

    /**
     *  Quantize a spin-independent one-electron operator in this general spinor basis. Spin-independent one-electron operators are those whose two-component matrix operator form contains the same scalar operator in the top-left and bottom-right corner.
     *
     *  @param fq_one_op            A spin-independent first-quantized operator.
     *
     *  @return The second-quantized representation of the given operator.
     */
    template <typename FQOneElectronOperator, typename Z = Shell>
    auto quantize(const FQOneElectronOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQOneElectronOperator<product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>, typename FQOneElectronOperator::Vectorizer>> {

        using ResultScalar = product_t<typename FQOneElectronOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQOneElectronOperator<ResultScalar, typename FQOneElectronOperator::Vectorizer>;
        using Vectorizer = typename FQOneElectronOperator::Vectorizer;

        const auto N = FQOneElectronOperator::NumberOfComponents;
        const auto& vectorizer = FQOneElectronOperator::vectorizer;


        // The strategy for calculating the matrix representation of the one-electron operator in this spinor basis is to:
        //  1. Express the operator in the underlying scalar bases and;
        //  2. Afterwards transform them using the current coefficient matrix.
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();


        // 1. Express the operator in the underlying scalar bases: spin-independent operators only have alpha-alpha and beta-beta blocks.
        auto engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(fq_one_op);
        const auto F_aa = GQCP::IntegralCalculator::calculate(engine, this->scalarBases().alpha().shellSet(), this->scalarBases().alpha().shellSet());
        const auto F_bb = GQCP::IntegralCalculator::calculate(engine, this->scalarBases().beta().shellSet(), this->scalarBases().beta().shellSet());

        // For each of the components of the operator, place the scalar basis representations into the spinor basis representation.
        std::array<SquareMatrix<ResultScalar>, N> fs;
        for (size_t i = 0; i < N; i++) {
            SquareMatrix<ResultScalar> f = SquareMatrix<ResultScalar>::Zero(M);

            f.topLeftCorner(K_alpha, K_alpha) = F_aa[i];
            f.bottomRightCorner(K_beta, K_beta) = F_bb[i];

            fs[i] = f;
        }


        // 2. Transform using the current coefficient matrix.
        StorageArray<SquareMatrix<ResultScalar>, Vectorizer> array {fs, vectorizer};
        ResultOperator op {array};  // 'op' for 'operator'.
        op.transform(this->expansion());
        return op;
    }


    /**
     *  Quantize the electronic spin operator in this general spinor basis.
     *
     *  @param fq_one_op        The (first-quantized) electronic spin operator.
     *
     *  @return The electronic spin operator expressed in this spinor basis.
     */
    template <typename Z = Shell>
    auto quantize(const ElectronicSpinOperator& fq_one_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQOneElectronOperator<product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>, ElectronicSpinOperator::Vectorizer>> {

        using ResultScalar = product_t<ElectronicSpinOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQOneElectronOperator<ResultScalar, ElectronicSpinOperator::Vectorizer>;

        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);
        const auto M = this->numberOfSpinors();

        // The strategy to quantize the spin operator is as follows.
        //  1. First, calculate the necessary overlap integrals over the scalar bases.
        //  2. Then, construct the scalar basis representations of the components of the spin operator by placing the overlaps into the correct blocks.
        //  3. Transform the components (in scalar basis) with the current coefficient matrix to yield the components in spinor basis.

        SquareMatrix<ResultScalar> S_x = SquareMatrix<ResultScalar>::Zero(M);
        SquareMatrix<ResultScalar> S_y = SquareMatrix<ResultScalar>::Zero(M);
        SquareMatrix<ResultScalar> S_z = SquareMatrix<ResultScalar>::Zero(M);


        // 1. Calculate the necessary overlap integrals over the scalar bases.
        auto overlap_engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(OverlapOperator());
        const auto S_aa = GQCP::IntegralCalculator::calculate(overlap_engine, this->scalarBases().alpha().shellSet(), this->scalarBases().alpha().shellSet())[0];
        const auto S_ab = GQCP::IntegralCalculator::calculate(overlap_engine, this->scalarBases().alpha().shellSet(), this->scalarBases().beta().shellSet())[0];
        const auto S_ba = GQCP::IntegralCalculator::calculate(overlap_engine, this->scalarBases().beta().shellSet(), this->scalarBases().alpha().shellSet())[0];
        const auto S_bb = GQCP::IntegralCalculator::calculate(overlap_engine, this->scalarBases().beta().shellSet(), this->scalarBases().beta().shellSet())[0];


        // 2. Place the overlaps into the correct blocks.
        S_x.block(0, K_alpha, K_alpha, K_beta) = 0.5 * S_ab;
        S_x.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * S_ba;

        using namespace GQCP::literals;
        S_y.block(0, K_alpha, K_alpha, K_beta) = -0.5 * 1_ii * S_ab;
        S_y.block(K_alpha, 0, K_beta, K_alpha) = 0.5 * 1_ii * S_ba;

        S_z.topLeftCorner(K_alpha, K_alpha) = 0.5 * S_aa;
        S_z.bottomRightCorner(K_beta, K_beta) = -0.5 * S_bb;


        // 3. Transform using the coefficient matrix.
        ResultOperator spin_op {std::vector<SquareMatrix<ResultScalar>> {S_x, S_y, S_z}};  // 'op' for operator.
        spin_op.transform(this->expansion());
        return spin_op;
    }


    /**
     *  Quantize the orbital Zeeman operator in this general spinor basis.
     *
     *  @param op               The (first-quantized) orbital Zeeman operator.
     *
     *  @return The orbital Zeeman operator expressed in this spinor basis.
     */
    template <typename Z = Shell>
    auto quantize(const OrbitalZeemanOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQOneElectronOperator<product_t<OrbitalZeemanOperator::Scalar, ExpansionScalar>, OrbitalZeemanOperator::Vectorizer>> {

        // Return the orbital Zeeman operator as a contraction beween the magnetic field and the angular momentum operator.
        const auto L = this->quantize(op.angularMomentum());
        const auto& B = op.magneticField().strength();
        return 0.5 * L.dot(B);
    }


    /**
     *  Quantize the diamagnetic operator in this general spinor basis.
     *
     *  @param op               The (first-quantized) diamagnetic operator.
     *
     *  @return The diamagnetic operator expressed in this spinor basis.
     */
    template <typename Z = Shell>
    auto quantize(const DiamagneticOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQOneElectronOperator<product_t<DiamagneticOperator::Scalar, ExpansionScalar>, DiamagneticOperator::Vectorizer>> {

        using ResultScalar = product_t<DiamagneticOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQOneElectronOperator<ResultScalar, DiamagneticOperator::Vectorizer>;


        // Return the diamagnetic operator as a contraction beween the magnetic field and the electronic quadrupole operator.
        const auto Q = this->quantize(op.electronicQuadrupole()).allParameters();
        const auto& B = op.magneticField().strength();

        // Prepare some variables.
        const auto& B_x = B(CartesianDirection::x);
        const auto& B_y = B(CartesianDirection::y);
        const auto& B_z = B(CartesianDirection::z);

        const auto& Q_xx = Q[DyadicCartesianDirection::xx];
        const auto& Q_xy = Q[DyadicCartesianDirection::xy];
        const auto& Q_xz = Q[DyadicCartesianDirection::xz];
        const auto& Q_yy = Q[DyadicCartesianDirection::yy];
        const auto& Q_yz = Q[DyadicCartesianDirection::yz];
        const auto& Q_zz = Q[DyadicCartesianDirection::zz];


        SquareMatrix<ResultScalar> D_par = 0.125 * ((std::pow(B_y, 2) + std::pow(B_z, 2)) * Q_xx +
                                                    (std::pow(B_x, 2) + std::pow(B_z, 2)) * Q_yy +
                                                    (std::pow(B_x, 2) + std::pow(B_y, 2)) * Q_zz -
                                                    2 * B_x * B_y * Q_xy -
                                                    2 * B_x * B_z * Q_xz -
                                                    2 * B_y * B_z * Q_yz);

        return ResultOperator {D_par};
    }


    /**
     *  Quantize the spin Zeeman operator in this general spinor basis.
     *
     *  @param op               The (first-quantized) spin Zeeman  operator.
     *
     *  @return The diamagnetic operator expressed in this spinor basis.
     */
    template <typename Z = Shell>
    auto quantize(const SpinZeemanOperator& op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQOneElectronOperator<product_t<SpinZeemanOperator::Scalar, ExpansionScalar>, SpinZeemanOperator::Vectorizer>> {

        // Return the spin Zeeman operator as a contraction beween the magnetic field and the spin operator.
        const auto S = this->quantize(ElectronicSpinOperator());
        const auto& B = op.magneticField().strength();
        return S.dot(B);
    }


    /**
     *  Quantize the Coulomb operator in this general spinor basis.
     *
     *  @param coulomb_op           The first-quantized Coulomb operator.
     *
     *  @return The second-quantized operator corresponding to the Coulomb operator.
     *
     *  @note For efficiency reasons (two-electron integrals over London orbitals are particularly slow to calculate), we assume that the alpha and beta scalar bases are equal.
     */
    template <typename Z = Shell>
    auto quantize(const CoulombRepulsionOperator& coulomb_op) const -> enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQTwoElectronOperator<product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>, CoulombRepulsionOperator::Vectorizer>> {

        using ResultScalar = product_t<CoulombRepulsionOperator::Scalar, ExpansionScalar>;
        using ResultOperator = GSQTwoElectronOperator<ResultScalar, CoulombRepulsionOperator::Vectorizer>;

        // The strategy for calculating the matrix representation of the two-electron operator in this spinor basis is to:
        //  1. Calculate the Coulomb integrals in the underlying scalar bases;
        //  2. Place the calculated integrals as 'blocks' in the larger representation, so that we can;
        //  3. Transform the operator using the current coefficient matrix.

        // 1. Calculate the Coulomb integrals in the underlying scalar bases.
        auto coulomb_engine = GQCP::IntegralEngine::InHouse<GQCP::LondonGTOShell>(CoulombRepulsionOperator());
        const auto g = GQCP::IntegralCalculator::calculate(coulomb_engine, this->scalarBases().alpha().shellSet(), this->scalarBases().alpha().shellSet())[0];


        // 2. Place the calculated integrals as 'blocks' in the larger representation.
        const auto K_alpha = this->numberOfCoefficients(Spin::alpha);
        const auto K_beta = this->numberOfCoefficients(Spin::beta);

        const auto M = this->numberOfSpinors();
        auto g_par = SquareRankFourTensor<ResultScalar>::Zero(M);  // 'par' for 'parameters'.

        // Primed indices are indices in the larger representation, normal ones are those in the smaller tensors.
        for (size_t mu_ = 0; mu_ < M; mu_++) {  // Mu 'prime'.
            const size_t mu = mu_ % K_alpha;

            for (size_t nu_ = 0; nu_ < M; nu_++) {  // Nu 'prime'.
                const size_t nu = nu_ % K_alpha;

                for (size_t rho_ = 0; rho_ < M; rho_++) {  // Rho 'prime'.
                    const size_t rho = rho_ % K_alpha;

                    for (size_t lambda_ = 0; lambda_ < M; lambda_++) {  // Lambda 'prime'.
                        const size_t lambda = lambda_ % K_alpha;

                        if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g(mu, nu, rho, lambda);
                        } else if ((mu_ < K_alpha) && (nu_ < K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ < K_alpha) && (lambda_ < K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g(mu, nu, rho, lambda);
                        } else if ((mu_ >= K_alpha) && (nu_ >= K_alpha) && (rho_ >= K_alpha) && (lambda_ >= K_alpha)) {
                            g_par(mu_, nu_, rho_, lambda_) = g(mu, nu, rho, lambda);
                        }
                    }
                }
            }
        }


        // 3. Transform the operator using the current coefficient matrix.
        ResultOperator g_op {g_par};  // 'op' for 'operator'.
        g_op.transform(this->expansion());
        return g_op;
    }


    /**
     *  Quantize the molecular magnetic Hamiltonian.
     *
     *  @param fq_hamiltonian           The molecular magnetic Hamiltonian.
     *
     *  @return The second-quantized molecular magnetic Hamiltonian.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQHamiltonian<ExpansionScalar>> quantize(const FQMolecularMagneticHamiltonian& fq_hamiltonian) const {

        const auto T = this->quantize(fq_hamiltonian.kinetic());
        const auto OZ = this->quantize(fq_hamiltonian.orbitalZeeman());
        const auto D = this->quantize(fq_hamiltonian.diamagnetic());

        const auto V = this->quantize(fq_hamiltonian.nuclearAttraction());

        const auto g = this->quantize(fq_hamiltonian.coulombRepulsion());

        return GSQHamiltonian<ExpansionScalar> {{T, OZ, D, V}, {g}};
    }


    /**
     *  Quantize the molecular Pauli Hamiltonian.
     *
     *  @param fq_hamiltonian           The molecular Pauli Hamiltonian.
     *
     *  @return The second-quantized molecular Pauli Hamiltonian.
     */
    template <typename Z = Shell>
    enable_if_t<std::is_same<Z, LondonGTOShell>::value, GSQHamiltonian<ExpansionScalar>> quantize(const FQMolecularPauliHamiltonian& fq_hamiltonian) const {

        const auto T = this->quantize(fq_hamiltonian.kinetic());
        const auto OZ = this->quantize(fq_hamiltonian.orbitalZeeman());
        const auto D = this->quantize(fq_hamiltonian.diamagnetic());

        const auto SZ = this->quantize(fq_hamiltonian.spinZeeman());

        const auto V = this->quantize(fq_hamiltonian.nuclearAttraction());

        const auto g = this->quantize(fq_hamiltonian.coulombRepulsion());

        return GSQHamiltonian<ExpansionScalar> {{T, OZ, D, SZ, V}, {g}};
    }


    /*
     *  MARK: Quantization of first-quantized operators
     */

    /**
     *  Quantize the molecular Hamiltonian.
     *
     *  @param fq_hamiltonian           The molecular Hamiltonian.
     *
     *  @return The second-quantized molecular Hamiltonian.
     */
    GSQHamiltonian<ExpansionScalar> quantize(const FQMolecularHamiltonian& fq_hamiltonian) const {

        const auto T = this->quantize(fq_hamiltonian.kinetic());
        const auto V = this->quantize(fq_hamiltonian.nuclearAttraction());

        const auto g = this->quantize(fq_hamiltonian.coulombRepulsion());

        return GSQHamiltonian<ExpansionScalar> {T + V, g};
    }


    /**
     *  Quantize the electronic S^2 operator in this spinor basis.
     *
     *  @param fq_S2_op             The first-quantized electronic S^2 operator.
     *
     *  @return The second-quantized representation of the electronic S^2 operator in this spinor basis.
     *
     *  @note Due to the way one-electron products are implemented, this quantization will only provide the correct result if the spinor basis is orthonormal.
     */
    ScalarGSQOneElectronOperatorProduct<ExpansionScalar> quantize(const ElectronicSpinSquaredOperator& fq_S2_op) const {

        // Prepare the 3 vector components of the electronic spin operator S = (S_x, S_y, S_z).
        const auto S_op = this->quantize(ElectronicSpinOperator());
        const auto& S_x_op = S_op(CartesianDirection::x);
        const auto& S_y_op = S_op(CartesianDirection::y);
        const auto& S_z_op = S_op(CartesianDirection::z);


        // Calculate S^2 = S_x^2 + S_y^2 + S_z^2.
        const auto S_x_2 = S_x_op * S_x_op;
        const auto S_y_2 = S_y_op * S_y_op;
        const auto S_z_2 = S_z_op * S_z_op;

        return S_x_2 + S_y_2 + S_z_2;
    }


    /*
     *  MARK: Scalar basis
     */

    /**
     *  @return The scalar bases for the alpha and beta components of the spinors.
     */
    const SpinResolved<ScalarBasis<Shell>>& scalarBases() const { return this->scalar_bases; }


    /**
     *  MARK: Mulliken partitioning
     */

    /**
     *  Partition this set of generalized spinors according to the Mulliken partitioning scheme.
     *
     *  @param selector             A function that returns true for basis functions that should be included the Mulliken partitioning.
     *
     *  @return A `GMullikenPartitioning` for the AOs selected by the supplied selector function.
     *
     *  @note The underlying scalar bases are assumed to be equal.
     */
    GMullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const BasisFunction&)>& selector) const {

        // Assume the underlying scalar bases are equal, and proceed to work with the one for the alpha component.
        const auto ao_indices = this->scalarBases().alpha().basisFunctionIndices(selector);
        return GMullikenPartitioning<ExpansionScalar> {ao_indices, this->expansion()};
    }


    /**
     *  Partition this set of generalized spinors according to the Mulliken partitioning scheme.
     *
     *  @param selector             A function that returns true for shells that should be included the Mulliken partitioning.
     *
     *  @return A `GMullikenPartitioning` for the AOs selected by the supplied selector function.
     *
     *  @note The underlying scalar bases are assumed to be equal.
     */
    GMullikenPartitioning<ExpansionScalar> mullikenPartitioning(const std::function<bool(const Shell&)>& selector) const {

        // Assume the underlying scalar bases are equal, and proceed to work with the one for the alpha component.
        const auto ao_indices = this->scalarBases().alpha().basisFunctionIndices(selector);
        return GMullikenPartitioning<ExpansionScalar> {ao_indices, this->expansion()};
    }
};


/*
 *  MARK: SpinorBasisTraits
 */

/**
 *  A type that provides compile-time information on spinor bases that is otherwise not accessible through a public class alias.
 */
template <typename _ExpansionScalar, typename _Shell>
struct SpinorBasisTraits<GSpinorBasis<_ExpansionScalar, _Shell>> {
    // The scalar type used to represent an expansion coefficient of the spinors in the underlying scalar orbitals: real or complex.
    using ExpansionScalar = _ExpansionScalar;

    // The type of transformation that is naturally related to a `GSpinorBasis`.
    using Transformation = GTransformation<ExpansionScalar>;

    // The second-quantized representation of the overlap operator related to the derived spinor basis.
    using SQOverlapOperator = ScalarGSQOneElectronOperator<ExpansionScalar>;
};


/*
 *  MARK: BasisTransformableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `BasisTransformable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct BasisTransformableTraits<GSpinorBasis<_ExpansionScalar, _Shell>> {

    // The type of transformation that is naturally related to a `GSpinorBasis`.
    using Transformation = GTransformation<_ExpansionScalar>;
};


/*
 *  MARK: JacobiRotatableTraits
 */

/**
 *  A type that provides compile-time information related to the abstract interface `JacobiRotatable`.
 */
template <typename _ExpansionScalar, typename _Shell>
struct JacobiRotatableTraits<GSpinorBasis<_ExpansionScalar, _Shell>> {

    // The type of Jacobi rotation that is naturally related to a `GSpinorBasis`.
    using JacobiRotationType = JacobiRotation;
};


}  // namespace GQCP
