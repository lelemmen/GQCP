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

#define BOOST_TEST_MODULE "CCSD"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Basis/SpinorBasis/USpinOrbitalBasis.hpp"
#include "Basis/Transformations/transform.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemSolver.hpp"
#include "ONVBasis/SpinUnresolvedONV.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"
#include "Operator/FirstQuantized/NuclearRepulsionOperator.hpp"
#include "Operator/SecondQuantized/SQHamiltonian.hpp"
#include "QCMethod/CC/CCSD.hpp"
#include "QCMethod/CC/CCSDEnvironment.hpp"
#include "QCMethod/CC/CCSDSolver.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "QCMethod/HF/GHF/GHF.hpp"
#include "QCMethod/HF/GHF/GHFSCFEnvironment.hpp"
#include "QCMethod/HF/GHF/GHFSCFSolver.hpp"
#include "QCMethod/HF/RHF/DiagonalRHFFockMatrixObjective.hpp"
#include "QCMethod/HF/RHF/RHF.hpp"
#include "QCMethod/HF/RHF/RHFSCFSolver.hpp"
#include "Utilities/complex.hpp"


/**
 *  Check if the implementation of spinor-CCSD is correct, by comparing with a reference by crawdad (https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2305).
 *
 *  The system under consideration is H2O in an STO-3G basisset.
 */
BOOST_AUTO_TEST_CASE(h2o_crawdad) {

    // Prepare the canonical RHF spin-orbital basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2o_crawdad.xyz");
    const auto N = molecule.numberOfElectrons();

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> r_spinor_basis {molecule, "STO-3G"};
    const auto r_sq_hamiltonian = r_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // in an AO basis.
    const auto K = r_spinor_basis.numberOfSpatialOrbitals();

    auto rhf_environment = GQCP::RHFSCFEnvironment<double>::WithCoreGuess(N, r_sq_hamiltonian, r_spinor_basis.overlap().parameters());
    auto plain_rhf_scf_solver = GQCP::RHFSCFSolver<double>::Plain();
    const GQCP::DiagonalRHFFockMatrixObjective<double> objective {r_sq_hamiltonian};
    const auto rhf_qc_structure = GQCP::QCMethod::RHF<double>().optimize(objective, plain_rhf_scf_solver, rhf_environment);
    const auto rhf_parameters = rhf_qc_structure.groundStateParameters();

    r_spinor_basis.transform(rhf_parameters.expansion());


    // Check if the intermediate RHF results are correct. We can't continue if this isn't the case.
    const auto rhf_energy = rhf_qc_structure.groundStateEnergy() + GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();
    const double ref_rhf_energy = -74.942079928192;
    BOOST_REQUIRE(std::abs(rhf_energy - ref_rhf_energy) < 1.0e-09);


    // Create a GSpinorBasis since we have implement spinor-CCSD, and quantize the molecular Hamiltonian in it.
    const auto g_spinor_basis = GQCP::GSpinorBasis<double, GQCP::GTOShell>::FromRestricted(r_spinor_basis);
    const auto M = g_spinor_basis.numberOfSpinors();
    const auto g_sq_hamiltonian = g_spinor_basis.quantize(GQCP::FQMolecularHamiltonian(molecule));  // In the canonical restricted spin-orbitals.

    // Create the GHF ONV (which is actually just the RHF ONV, since we're using the canonical RHF orbitals) and the corresponding orbital space.
    const auto reference_onv = GQCP::SpinUnresolvedONV::GHF(2 * K, N, rhf_parameters.spinOrbitalEnergiesBlocked());
    const auto orbital_space = reference_onv.orbitalSpace();


    // Initialize an environment suitable for CCSD.
    auto environment = GQCP::CCSDEnvironment<double>::PerturbativeCCSD(g_sq_hamiltonian, orbital_space);

    // Since we're working with a Hartree-Fock reference, the perturbative amplitudes actually correspond to the MP2 amplitudes. This means that the initial CCSD energy correction is the MP2 energy correction.
    const double ref_mp2_correction_energy = -0.049149636120;
    const auto& t1 = environment.t1_amplitudes.back();

    BOOST_REQUIRE(t1.asImplicitMatrixSlice().asMatrix().isZero(1.0e-08));  // for a HF reference, the perturbative T1 amplitudes are zero

    const auto initial_ccsd_correction_energy = environment.correlation_energies.back();
    BOOST_REQUIRE(std::abs(initial_ccsd_correction_energy - ref_mp2_correction_energy) < 1.0e-10);


    // Prepare the CCSD solver and optimize the CCSD model parameters.
    auto solver = GQCP::CCSDSolver<double>::Plain();
    const auto ccsd_qc_structure = GQCP::QCMethod::CCSD<double>().optimize(solver, environment);

    const auto ccsd_correlation_energy = ccsd_qc_structure.groundStateEnergy();

    const double ref_ccsd_correlation_energy = -0.070680088376;
    BOOST_CHECK(std::abs(ccsd_correlation_energy - ref_ccsd_correlation_energy) < 1.0e-08);
}


/**
 *  Validate the correctness of complex-valued CCSD by checking it against the complex FCI results for H2.
 * 
 *  The calculations are performed in a 6-31G basis set, using London orbitals to describe the magnetic Hamiltonian.
 */
BOOST_AUTO_TEST_CASE(CCSD_vs_CI_complex) {

    // Set up the magnetic Hamiltonian in the GHF canonical orbital basis.
    const auto molecule = GQCP::Molecule::ReadXYZ("data/h2.xyz");
    const auto N = molecule.numberOfElectrons();
    const auto nuclear_repulsion = GQCP::NuclearRepulsionOperator(molecule.nuclearFramework()).value();

    const GQCP::HomogeneousMagneticField B {{0, 0, 0.1}};


    GQCP::GSpinorBasis<GQCP::complex, GQCP::LondonGTOShell> spinor_basis {molecule, "6-31G", B};
    const auto M = spinor_basis.numberOfSpinors();
    const auto S = spinor_basis.overlap();
    auto hamiltonian = spinor_basis.quantize(GQCP::FQMolecularPauliHamiltonian(molecule, B));  // Here, the magnetic Hamiltonian is expressed in the AO basis.


    // Solve the GHF SCF equations.
    auto ghf_environment = GQCP::GHFSCFEnvironment<GQCP::complex>::WithComplexlyTransformedCoreGuess(N, hamiltonian, S);
    auto scf_solver = GQCP::GHFSCFSolver<GQCP::complex>::Plain();
    const auto ghf_qc_structure = GQCP::QCMethod::GHF<GQCP::complex>().optimize(scf_solver, ghf_environment);
    const auto& ghf_parameters = ghf_qc_structure.groundStateParameters();
    const auto rhf_energy = ghf_qc_structure.groundStateEnergy();

    const auto& C = ghf_parameters.expansion();
    spinor_basis.transform(C);
    hamiltonian.transform(C);


    // Calculate the CCSD electronic energy.
    const auto orbital_space = ghf_parameters.orbitalSpace();
    auto ccsd_environment = GQCP::CCSDEnvironment<GQCP::complex>::PerturbativeCCSD(hamiltonian, orbital_space);


    // Prepare the CCSD solver and optimize the CCSD model parameters.
    auto ccsd_solver = GQCP::CCSDSolver<GQCP::complex>::Plain();
    const auto ccsd_qc_structure = GQCP::QCMethod::CCSD<GQCP::complex>().optimize(ccsd_solver, ccsd_environment);
    const auto ccsd_correlation_energy = ccsd_qc_structure.groundStateEnergy();
    const auto ccsd_electronic_energy = ccsd_correlation_energy + rhf_energy;


    // Calculate the FCI electronic energy and check if it matches the CCSD result.
    const GQCP::SpinUnresolvedONVBasis full_onv_basis {M, N};
    const GQCP::SpinUnresolvedSelectedONVBasis onv_basis {full_onv_basis};

    auto fci_environment = GQCP::CIEnvironment::Dense(hamiltonian, onv_basis);
    auto fci_solver = GQCP::EigenproblemSolver::Dense<GQCP::complex>();
    auto fci_qc_structure = GQCP::QCMethod::CI<GQCP::complex, GQCP::SpinUnresolvedSelectedONVBasis>(onv_basis).optimize(fci_solver, fci_environment);

    const auto fci_electronic_energy = fci_qc_structure.groundStateEnergy();

    BOOST_CHECK(std::abs(ccsd_electronic_energy - fci_electronic_energy) < 1.0e-08);
}
