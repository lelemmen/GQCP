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

#define BOOST_TEST_MODULE "SpinResolvedONV"

#include <boost/test/unit_test.hpp>

#include "ONVBasis/SpinResolvedONV.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "QCModel/HF/UHF.hpp"


/**
 *  Check if the creation of the 'RHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(RHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_P=3 electron electron pairs, the 'RHF' ONV should be ("00111" = 7) ("00111" = 7)
    const size_t K = 5;
    const size_t N_P = 3;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_P, 7};

    const GQCP::SpinResolvedONV reference {alpha_onv, alpha_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::RHF(K, N_P));
}


/**
 *  Check if the creation of the 'UHF' ONV works as expected.
 */
BOOST_AUTO_TEST_CASE(UHF) {

    // For K=5 alpha-spinorbitals, K=5 beta-spinorbitals, and N_alpha=3 alpha-electrons and N_beta=2 beta-electrons, the 'UHF' ONV should be ("00111" = 7) ("00011" = 3)
    const size_t K = 5;
    const size_t N_alpha = 3;
    const size_t N_beta = 2;
    const GQCP::SpinUnresolvedONV alpha_onv {K, N_alpha, 7};
    const GQCP::SpinUnresolvedONV beta_onv {K, N_beta, 3};

    const GQCP::SpinResolvedONV reference {alpha_onv, beta_onv};
    BOOST_CHECK(reference == GQCP::SpinResolvedONV::UHF(K, N_alpha, N_beta));
}


/**
 *  Check if the overlap between an RHF-related ONV and an UHF-related ONV works as expected.
 * 
 *  We don't really have a reference implementation, but we can check if the overlaps are equal to 1 if we use RHF orbitals and UHF orbitals that have the same alpha- and beta-part.
 * 
 *  The system under consideration is H2 with a 6-31G** basisset.
 */
BOOST_AUTO_TEST_CASE(RHF_UHF_overlap) {

    // Find the canonical RHF parameters.
    const auto h2 = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N_P = h2.numberOfElectrons() / 2;

    const GQCP::RSpinorBasis<double, GQCP::GTOShell> r_spinor_basis {h2, "6-31G**"};
    const auto K = r_spinor_basis.numberOfSpatialOrbitals();

    auto sq_hamiltonian = GQCP::SQHamiltonian<double>::Molecular(r_spinor_basis, h2);  // in an AO basis

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(h2.numberOfElectrons(), sq_hamiltonian, r_spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {sq_hamiltonian};

    const auto rhf_parameters = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment).groundStateParameters();


    // Convert the RHF parameters into UHF parameters.
    const GQCP::QCModel::UHF<double> uhf_parameters {rhf_parameters};
    const auto u_spinor_basis = GQCP::USpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);


    // Check if the RHF determinant has overlap 1 with the corresponding UHF determinant.
    const auto rhf_determinant = GQCP::SpinResolvedONV::RHF(K, N_P);
    const auto uhf_determinant = GQCP::SpinResolvedONV::UHF(K, N_P, N_P);
    BOOST_CHECK(std::abs(rhf_determinant.calculateOverlap(uhf_determinant, r_spinor_basis, u_spinor_basis) - 1.0) < 1.0e-07);
}
