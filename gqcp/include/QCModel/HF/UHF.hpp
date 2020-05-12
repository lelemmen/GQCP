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


#include "Basis/SpinorBasis/Spin.hpp"
#include "Basis/TransformationMatrix.hpp"
#include "Mathematical/Representation/Matrix.hpp"
#include "Processing/RDM/OneRDM.hpp"
#include "QCModel/HF/RHF.hpp"


namespace GQCP {
namespace QCModel {


/**
 *  The unrestricted Hartree-Fock wave function model.
 * 
 *  @tparam _Scalar         the scalar representation of the coefficients in the coefficient matrix
 */
template <typename _Scalar>
class UHF {
public:
    using Scalar = _Scalar;


private:
    size_t N_alpha;  // the number of alpha electrons
    size_t N_beta;   // the number of beta electrons

    VectorX<double> orbital_energies_alpha;  // sorted by ascending energy
    VectorX<double> orbital_energies_beta;   // sorted by ascending energy

    TransformationMatrix<Scalar> C_alpha;  // the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
    TransformationMatrix<Scalar> C_beta;   // the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  A basic constructor that sets all the properties.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C_alpha                                  the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
     *  @param C_beta                                   the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis
     */
    UHF(const size_t N_alpha, const size_t N_beta, const VectorX<double>& orbital_energies_alpha, const VectorX<double>& orbital_energies_beta, const TransformationMatrix<Scalar>& C_alpha, const TransformationMatrix<Scalar>& C_beta) :
        N_alpha {N_alpha},
        N_beta {N_beta},
        orbital_energies_alpha {orbital_energies_alpha},
        orbital_energies_beta {orbital_energies_beta},
        C_alpha {C_alpha},
        C_beta {C_beta} {

        // Check for valid arguments.
        const auto K_alpha = C_alpha.dimension();  // number of alpha spatial orbitals
        const auto K_beta = C_beta.dimension();    // number of beta spatial orbitals

        if (N_alpha > K_alpha) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given alpha electrons cannot be larger than the number of alpha spatial orbitals.");
        }

        if (N_beta > K_beta) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given beta electrons cannot be larger than the number of beta spatial orbitals.");
        }

        if (K_alpha != orbital_energies_alpha.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given alpha-spin-orbital energies does not match the number of alpha spin-orbitals.");
        }

        if (K_beta != orbital_energies_beta.size()) {
            throw std::invalid_argument("UHF(const size_t, const size_t, const VectorX<double>&, const VectorX<double>&, const TransformationMatrix<Scalar>&, const TransformationMatrix<Scalar>&): The number of given beta-spin-orbital energies does not match the beta of beta spin-orbitals.");
        }
    }


    /**
     *  A constructor that initializes the orbital energies to zeros.
     * 
     *  @param N_alpha                                  the number of alpha electrons
     *  @param N_beta                                   the number of beta electrons
     *  @param orbital_energies_alpha                   the orbital energies for the alpha-spin-orbitals, sorted by ascending energy
     *  @param orbital_energies_beta                    the orbital energies for the beta-spin-orbitals, sorted by ascending energy
     *  @param C_alpha                                  the coefficient matrix that expresses every alpha spatial orbital (as a column) in its underlying scalar basis
     *  @param C_beta                                   the coefficient matrix that expresses every beta spatial orbital (as a column) in its underlying scalar basis
     */
    UHF(const size_t N_alpha, const size_t N_beta, const TransformationMatrix<Scalar>& C_alpha, const TransformationMatrix<Scalar>& C_beta) :
        UHF(N_alpha, N_beta,
            GQCP::VectorX<double>::Zero(C_alpha.dimension()),
            GQCP::VectorX<double>::Zero(C_beta.dimension()),
            C_alpha, C_beta) {
    }


    /**
     *  Convert an RHF wave function model to an UHF wave function model.
     * 
     *  @param rhf_model            an RHF wave function model
     */
    UHF(const GQCP::QCModel::RHF<Scalar>& rhf_model) :
        UHF(rhf_model.numberOfElectrons(Spin::alpha), rhf_model.numberOfElectrons(Spin::beta),
            rhf_model.orbitalEnergies(), rhf_model.orbitalEnergies(),
            rhf_model.coefficientMatrix(), rhf_model.coefficientMatrix()) {}


    /*
     *  PUBLIC STATIC METHODS
     */

    /**
     *  @param P_sigma              the spin-sigma density matrix in a scalar basis
     *  @param H_core_sigma         the spin-sigma core Hamiltonian expressed in the same scalar basis
     *  @param F_sigma              the spin-sigma Fock matrix in the same scalar basis
     *
     *  @return the UHF electronic energy for the sigma electrons
     */
    static double calculateElectronicEnergy(const OneRDM<Scalar>& P_sigma, const ScalarSQOneElectronOperator<Scalar>& H_core_sigma, const ScalarSQOneElectronOperator<Scalar>& F_sigma) {

        // First, calculate the sum of H_core and F (this saves a contraction).
        const ScalarSQOneElectronOperator<Scalar> Z_sigma = H_core_sigma + F_sigma;

        // Convert the matrices Z and P to a tensor, as contractions are only implemented for tensors.
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_sigma_tensor {P_sigma.data(), P_sigma.rows(), P_sigma.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> Z_sigma_tensor {Z_sigma.parameters().data(), Z_sigma.parameters().rows(), Z_sigma.parameters().cols()};

        // Specify the contraction pair
        // To calculate the electronic energy, we must perform a double contraction.
        //      0.5 P_sigma(mu nu) P_sigma(mu nu)
        Eigen::array<Eigen::IndexPair<int>, 2> contraction_pair = {Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1)};

        // Calculate the double contraction (with prefactor 0.5).
        Tensor<Scalar, 0> contraction = 0.5 * P_sigma_tensor.contract(Z_sigma_tensor, contraction_pair);

        // As the double contraction of two rank-2 tensors is a scalar (a tensor of rank 0), we should access the value as (0).
        return contraction(0);
    }


    /**
     *  @param K_sigma          the number of spatial orbitals for one of the spin components
     *  @param N_sigma          the number of sigma electrons, i.e. the number of occupied sigma spin-orbitals
     *
     *  @return the sigma-spin UHF 1-RDM expressed in an orthonormal sigma spin-orbital basis
     */
    static OneRDM<Scalar> calculateOrthonormalBasis1RDM(const size_t K_sigma, const size_t N_sigma) {

        // The 1-RDM for UHF looks like (for K=5, N=3)
        //    1  0  0  0  0
        //    0  1  0  0  0
        //    0  0  1  0  0
        //    0  0  0  0  0
        //    0  0  0  0  0

        OneRDM<Scalar> D_MO = OneRDM<Scalar>::Zero(K_sigma, K_sigma);
        D_MO.topLeftCorner(N_sigma, N_sigma) = SquareMatrix<Scalar>::Identity(N_sigma, N_sigma);

        return D_MO;
    }


    /**
     *  @param C_sigma          the coefficient matrix that expresses the sigma spin-orbitals (as a column) in its underlying scalar basis
     *  @param N_sigma          the number of sigma electrons, i.e. the number of occupied sigma spin-orbitals
     *
     *  @return the sigma-spin UHF 1-RDM expressed in the underlying scalar basis
     */
    static OneRDM<Scalar> calculateScalarBasis1RDM(const TransformationMatrix<double>& C_sigma, const size_t N_sigma) {

        const auto K_sigma = C_sigma.dimension();
        const auto D_orthonormal = UHF<Scalar>::calculateOrthonormalBasis1RDM(K_sigma, N_sigma);

        // Transform the 1-RDM in an orthonormal basis to the underlying scalar basis
        return C_sigma.conjugate() * D_orthonormal * C_sigma.transpose();
    }


    /**
     *  Calculate the UHF direct (Coulomb) matrix for spin sigma.
     * 
     *  @param P_alpha              the UHF alpha density matrix expressed in the underlying scalar orbital basis
     *  @param P_beta               the UHF beta density matrix expressed in the (same) underlying scalar orbital basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the (same) underlying scalar orbital basis
     * 
     *  @return the UHF direct (Coulomb) matrix for spin sigma
     */
    static ScalarSQOneElectronOperator<Scalar> calculateScalarBasisDirectMatrix(const OneRDM<Scalar>& P_alpha, const OneRDM<Scalar>& P_beta, const SQHamiltonian<Scalar>& sq_hamiltonian) {

        // To perform the contraction, we will first have to convert the density matrices into tensors (since contractions are only implemented for tensors).
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_alpha_tensor {P_alpha.data(), P_alpha.rows(), P_alpha.cols()};
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_beta_tensor {P_beta.data(), P_beta.rows(), P_beta.cols()};

        // Specify the contraction pairs for the direct contractions:
        //      (mu nu|rho lambda) P(rho lambda)
        Eigen::array<Eigen::IndexPair<int>, 2> direct_contraction_pair = {Eigen::IndexPair<int>(2, 0), Eigen::IndexPair<int>(3, 1)};

        // Do the actual contractions, and convert the given tensor back to a matrix.
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        Tensor<Scalar, 2> J_alpha_tensor = g.contract(P_alpha_tensor, direct_contraction_pair);
        Tensor<Scalar, 2> J_beta_tensor = g.contract(P_beta_tensor, direct_contraction_pair);

        Tensor<Scalar, 2> J_tensor = J_alpha_tensor.Eigen() + J_beta_tensor.Eigen();

        Eigen::Map<Eigen::MatrixXd> J {J_tensor.data(), J_tensor.dimension(0), J_tensor.dimension(1)};
        return ScalarSQOneElectronOperator<Scalar>(J);
    }


    /**
     *  Calculate the UHF exchange matrix for spin sigma.
     * 
     *  @param P_sigma              the UHF sigma density matrix expressed in the underlying scalar orbital basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the (same) underlying scalar orbital basis
     * 
     *  @return the UHF direct (Coulomb) matrix for spin sigma
     */
    static ScalarSQOneElectronOperator<Scalar> calculateScalarBasisExchangeMatrix(const OneRDM<Scalar>& P_sigma, const SQHamiltonian<Scalar>& sq_hamiltonian) {

        // To perform the contraction, we will first have to convert the density matrix into a tensor (since contractions are only implemented for tensors).
        Eigen::TensorMap<Eigen::Tensor<const Scalar, 2>> P_sigma_tensor {P_sigma.data(), P_sigma.rows(), P_sigma.cols()};

        // Specify the contraction pairs for the exchange contraction:
        //      (mu rho|lambda nu) P(lambda rho)
        Eigen::array<Eigen::IndexPair<int>, 2> exchange_contraction_pair = {Eigen::IndexPair<int>(1, 1), Eigen::IndexPair<int>(2, 0)};

        // Do the actual contraction, and convert the given tensor back to a matrix.
        const auto& g = sq_hamiltonian.twoElectron().parameters();
        Tensor<Scalar, 2> K_sigma_tensor = g.contract(P_sigma_tensor, exchange_contraction_pair);

        Eigen::Map<Eigen::MatrixXd> K_sigma {K_sigma_tensor.data(), K_sigma_tensor.dimension(0), K_sigma_tensor.dimension(1)};
        return ScalarSQOneElectronOperator<Scalar>(K_sigma);
    }


    /**
     *  Calculate the UHF Fock matrix F = H_core + G, in which G is a contraction of the density matrix and the two-electron integrals
     *
     *  @param P                    the RHF density matrix in a scalar basis
     *  @param sq_hamiltonian       the Hamiltonian expressed in the same scalar basis
     *
     *  @return the RHF Fock matrix expressed in the scalar basis
     */
    static ScalarSQOneElectronOperator<Scalar> calculateScalarBasisFockMatrix(const Spin sigma, const OneRDM<Scalar>& P_alpha, const OneRDM<Scalar>& P_beta, const SQHamiltonian<Scalar>& sq_hamiltonian) {

        // F_sigma = H_core + (J_alpha + J_beta) - K_sigma
        const auto& H_core = sq_hamiltonian.core();
        const auto J = UHF<Scalar>::calculateScalarBasisDirectMatrix(P_alpha, P_beta, sq_hamiltonian);

        switch (sigma) {
        case Spin::alpha: {
            const auto K_alpha = UHF<Scalar>::calculateScalarBasisExchangeMatrix(P_alpha, sq_hamiltonian);
            return H_core + J - K_alpha;
            break;
        }

        case Spin::beta: {
            const auto K_beta = UHF<Scalar>::calculateScalarBasisExchangeMatrix(P_beta, sq_hamiltonian);
            return H_core + J - K_beta;
            break;
        }
        }
    }


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the sigma-spin UHF 1-RDM expressed in an orthonormal sigma spin-orbital basis for these UHF model parameters
     */
    OneRDM<Scalar> calculateOrthonormalBasis1RDM(const Spin sigma) const {

        const auto K_sigma = this->numberOfSpinOrbitals(sigma);
        const auto N_sigma = this->numberOfElectrons(sigma);

        return UHF<Scalar>::calculateOrthonormalBasis1RDM(K_sigma, N_sigma);
    }


    /**
     *  @param sigma            alpha or beta
     *
     *  @return the sigma-spin UHF 1-RDM expressed in the underlying scalar basis for these UHF model parameters
     */
    OneRDM<Scalar> calculateScalarBasis1RDM(const Spin sigma) const {

        const auto C_sigma = this->coefficientMatrix(sigma);
        const auto N_sigma = this->numberOfElectrons(sigma);

        return UHF<Scalar>::calculateScalarBasis1RDM(C_sigma, N_sigma);
    }


    /**
     *  @param sigma            alpha or beta
     *
     *  @return the coefficient matrix that expresses the sigma spin-orbitals (as a column) in its underlying scalar basis
     */
    const TransformationMatrix<Scalar>& coefficientMatrix(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->C_alpha;
            break;
        }

        case Spin::beta: {
            return this->C_beta;
            break;
        }
        }
    }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of sigma electrons that these UHF model parameters describe, i.e. the number of occupied sigma-spin-orbitals
     */
    size_t numberOfElectrons(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->N_alpha;
            break;
        }

        case Spin::beta: {
            return this->N_beta;
            break;
        }
        }
    }

    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the number of sigma spin-orbitals that these UHF model parameters describe
     */
    size_t numberOfSpinOrbitals(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->C_alpha.dimension();
            break;
        }

        case Spin::beta: {
            return this->C_beta.dimension();
            break;
        }
        }
    }


    /**
     *  @param sigma            alpha or beta
     * 
     *  @return the orbital energies of the sigma-spin-orbitals
     */
    const VectorX<double>& orbitalEnergies(const Spin sigma) const {

        switch (sigma) {
        case Spin::alpha: {
            return this->orbital_energies_alpha;
            break;
        }

        case Spin::beta: {
            return this->orbital_energies_beta;
            break;
        }
        }
    }
};


}  // namespace QCModel
}  // namespace GQCP
