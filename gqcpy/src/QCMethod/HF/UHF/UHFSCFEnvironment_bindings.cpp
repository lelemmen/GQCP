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

#include "QCMethod/HF/UHF/UHFSCFEnvironment.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


void bindUHFSCFEnvironment(py::module& module) {
    py::class_<GQCP::UHFSCFEnvironment<double>>(module, "UHFSCFEnvironment", "An algorithm environment that can be used with standard UHF SCF solvers.")

        // CONSTRUCTORS

        .def(py::init([](const size_t N_alpha, const size_t N_beta, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S, const Eigen::MatrixXd& C_alpha_initial, const Eigen::MatrixXd& C_beta_initial) {
                 return GQCP::UHFSCFEnvironment<double>(N_alpha, N_beta, sq_hamiltonian, GQCP::SquareMatrix<double>(S), GQCP::TransformationMatrix<double>(C_alpha_initial), GQCP::TransformationMatrix<double>(C_beta_initial));
             }),
             py::arg("N_alpha"),
             py::arg("N_beta"),
             py::arg("sq_hamiltonian"),
             py::arg("S"),
             py::arg("C_alpha_initial"),
             py::arg("C_beta_initial"),
             "A constructor that initializes the environment with initial guesses for the alpha and beta coefficient matrices.")

        .def(py::init([](const GQCP::QCModel::RHF<double>& rhf_parameters, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                 return GQCP::UHFSCFEnvironment<double>(rhf_parameters, sq_hamiltonian, GQCP::SquareMatrix<double>(S));
             }),
             py::arg("rhf_parameters"),
             py::arg("sq_hamiltonian"),
             py::arg("S"),
             "A constructor that initializes the environment from converged RHF model parameters.")

        .def_static(
            "WithCoreGuess",
            [](const size_t N_alpha, const size_t N_beta, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const Eigen::MatrixXd& S) {  // use an itermediary Eigen matrix for the Python binding, since Pybind11 doesn't accept our types that are derived from Eigen::Matrix
                return GQCP::UHFSCFEnvironment<double>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, GQCP::SquareMatrix<double>(S));
            },
            "Initialize an UHF SCF environment with initial coefficient matrices (equal for alpha and beta) that is obtained by diagonalizing the core Hamiltonian matrix.")


        // Bind read-write members/properties, exposing intermediary environment variables to the Python interface.
        .def_readwrite("N_alpha", &GQCP::UHFSCFEnvironment<double>::N_alpha)
        .def_readwrite("N_beta", &GQCP::UHFSCFEnvironment<double>::N_beta)

        .def_readwrite("electronic_energies", &GQCP::UHFSCFEnvironment<double>::electronic_energies)

        .def_readwrite("orbital_energies_alpha", &GQCP::UHFSCFEnvironment<double>::orbital_energies_alpha)
        .def_readwrite("orbital_energies_beta", &GQCP::UHFSCFEnvironment<double>::orbital_energies_beta)

        .def_property(
            "S",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                return environment.S;
            },
            [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& S) {
                environment.S = GQCP::SquareMatrix<double>(S);
            })


        // Define read-only 'getters'
        .def_readonly(
            "coefficient_matrices_alpha",
            &GQCP::UHFSCFEnvironment<double>::coefficient_matrices_alpha)

        .def_readonly(
            "coefficient_matrices_beta",
            &GQCP::UHFSCFEnvironment<double>::coefficient_matrices_beta)

        .def_readonly(
            "density_matrices",
            &GQCP::UHFSCFEnvironment<double>::density_matrices)

        .def_readonly(
            "error_vectors_alpha",
            &GQCP::UHFSCFEnvironment<double>::error_vectors_alpha)

        .def_readonly(
            "error_vectors_beta",
            &GQCP::UHFSCFEnvironment<double>::error_vectors_beta)

        // Define getters for non-native components
        .def(
            "density_matrices_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                std::vector<GQCP::MatrixX<double>> alpha_density_matrices;
                std::transform(environment.density_matrices.begin(), environment.density_matrices.end(), std::back_inserter(alpha_density_matrices), [](const GQCP::SpinResolvedOneDM<double>& D) { return D.alpha(); });

                return alpha_density_matrices;
            })

        .def(
            "density_matrices_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                std::vector<GQCP::MatrixX<double>> beta_density_matrices;
                std::transform(environment.density_matrices.begin(), environment.density_matrices.end(), std::back_inserter(beta_density_matrices), [](const GQCP::SpinResolvedOneDM<double>& D) { return D.beta(); });

                return beta_density_matrices;
            })

        .def(
            "fock_matrices_alpha",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                std::vector<GQCP::MatrixX<double>> alpha_fock_matrices;
                std::transform(environment.fock_matrices.begin(), environment.fock_matrices.end(), std::back_inserter(alpha_fock_matrices), [](const GQCP::ScalarUSQOneElectronOperator<double>& F) { return F.parameters(GQCP::Spin::alpha); });

                return alpha_fock_matrices;
            })

        .def(
            "fock_matrices_beta",
            [](const GQCP::UHFSCFEnvironment<double>& environment) {
                std::vector<GQCP::MatrixX<double>> beta_fock_matrices;
                std::transform(environment.fock_matrices.begin(), environment.fock_matrices.end(), std::back_inserter(beta_fock_matrices), [](const GQCP::ScalarUSQOneElectronOperator<double>& F) { return F.parameters(GQCP::Spin::beta); });

                return beta_fock_matrices;
            })

        // Bind methods for the replacement of the most current iterates.
        .def("replace_current_coefficient_matrix_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_coefficient_matrix_alpha) {
                 environment.coefficient_matrices_alpha.pop_back();
                 environment.coefficient_matrices_alpha.push_back(GQCP::SquareMatrix<double>(new_coefficient_matrix_alpha));
             })

        .def("replace_current_coefficient_matrix_beta",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_coefficient_matrix_beta) {
                 environment.coefficient_matrices_beta.pop_back();
                 environment.coefficient_matrices_beta.push_back(GQCP::SquareMatrix<double>(new_coefficient_matrix_beta));
             })

        .def("replace_current_density_matrix_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_density_matrix_alpha) {
                 const auto last_spin_resolved_density = environment.density_matrices.back();
                 const auto last_density_matrix_beta = last_spin_resolved_density.beta();
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(GQCP::SpinResolvedOneDM<double>(new_density_matrix_alpha, last_density_matrix_beta));
             })

        .def("replace_current_density_matrix_beta",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_density_matrix_beta) {
                 const auto last_spin_resolved_density = environment.density_matrices.back();
                 const auto last_density_matrix_alpha = last_spin_resolved_density.beta();
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(GQCP::SpinResolvedOneDM<double>(last_density_matrix_alpha, new_density_matrix_beta));
             })

        .def("replace_current_fock_matrix_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_fock_matrix_alpha) {
                 const auto last_USQ_Fock = environment.fock_matrices.back();
                 const auto last_fock_matrix_beta = last_USQ_Fock.parameters(GQCP::Spin::beta);
                 GQCP::ScalarUSQOneElectronOperator<double> new_fock {new_fock_matrix_alpha, last_fock_matrix_beta};

                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(new_fock);
             })

        .def("replace_current_fock_matrix_beta",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_fock_matrix_beta) {
                 const auto last_USQ_Fock = environment.fock_matrices.back();
                 const auto last_fock_matrix_alpha = last_USQ_Fock.parameters(GQCP::Spin::alpha);
                 GQCP::ScalarUSQOneElectronOperator<double> new_fock {last_fock_matrix_alpha, new_fock_matrix_beta};

                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(new_fock);
             })

        .def("replace_current_error_vectors_alpha",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_error_vectors_alpha) {
                 environment.error_vectors_alpha.pop_back();
                 environment.error_vectors_alpha.push_back(GQCP::SquareMatrix<double>(new_error_vectors_alpha));
             })

        .def("replace_current_error_vectors_beta",
             [](GQCP::UHFSCFEnvironment<double>& environment, const Eigen::MatrixXd& new_error_vectors_beta) {
                 environment.error_vectors_beta.pop_back();
                 environment.error_vectors_beta.push_back(GQCP::SquareMatrix<double>(new_error_vectors_beta));
             });
}


}  // namespace gqcpy
