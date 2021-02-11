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
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to UHF SCF environments.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `UHFSCFEnvironment`.
 */
template <typename Class>
void bindQCMethodUHFSCFEnvironmentInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::Scalar;


    py_class

        /*
         *  MARK: Named constructors specific to UHF
         */

        .def(py::init([](const size_t N_alpha, const size_t N_beta, const USQHamiltonian<Scalar>& sq_hamiltonian, const ScalarUSQOneElectronOperator<Scalar>& S, const UTransformation<Scalar>& C_initial) {
                 return UHFSCFEnvironment<Scalar>(N_alpha, N_beta, sq_hamiltonian, S, C_initial);
             }),
             py::arg("N_alpha"),
             py::arg("N_beta"),
             py::arg("sq_hamiltonian"),
             py::arg("S"),
             py::arg("C_initial"),
             "A constructor that initializes the environment with initial guesses for the alpha and beta coefficient matrices.")

        .def(py::init([](const QCModel::RHF<Scalar>& rhf_parameters, const USQHamiltonian<Scalar>& sq_hamiltonian, const ScalarUSQOneElectronOperator<Scalar>& S) {
                 return UHFSCFEnvironment<Scalar>(rhf_parameters, sq_hamiltonian, S);
             }),
             py::arg("rhf_parameters"),
             py::arg("sq_hamiltonian"),
             py::arg("S"),
             "A constructor that initializes the environment from converged RHF model parameters.")

        .def_static(
            "WithCoreGuess",
            [](const size_t N_alpha, const size_t N_beta, const USQHamiltonian<Scalar>& sq_hamiltonian, const ScalarUSQOneElectronOperator<Scalar>& S) {
                return UHFSCFEnvironment<Scalar>::WithCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
            },
            "Initialize an UHF SCF environment with initial coefficient matrices (equal for alpha and beta) that is obtained by diagonalizing the core Hamiltonian matrix.")

        .def_static(
            "WithTransformedCoreGuess",
            [](const size_t N_alpha, const size_t N_beta, const USQHamiltonian<Scalar>& sq_hamiltonian, const ScalarUSQOneElectronOperator<Scalar>& S, const std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>&)>& transformation_function) {
                return UHFSCFEnvironment<Scalar>::WithTransformedCoreGuess(N_alpha, N_beta, sq_hamiltonian, S, transformation_function);
            },
            "Initialize an UHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently applying the given unary transformation function.")

        /*
         *  MARK: Read-write members & properties specific to UHF
         */

        .def_property(
            "S",
            [](const UHFSCFEnvironment<Scalar>& environment) {
                return environment.S;
            },
            [](UHFSCFEnvironment<Scalar>& environment, const ScalarUSQOneElectronOperator<Scalar>& S) {
                environment.S = S;
            })

        /*
         *  MARK: Replacements specific to UHF
         */

        .def("replace_current_coefficient_matrices",
             [](UHFSCFEnvironment<Scalar>& environment, const Eigen::MatrixXd& new_coefficient_matrix_alpha, const Eigen::MatrixXd& new_coefficient_matrix_beta) {
                 environment.coefficient_matrices.pop_back();
                 environment.coefficient_matrices.push_back(UTransformation<Scalar>(UTransformationComponent<Scalar>(new_coefficient_matrix_alpha), UTransformationComponent<Scalar>(new_coefficient_matrix_beta)));
             })

        .def("replace_current_density_matrices",
             [](UHFSCFEnvironment<Scalar>& environment, const Eigen::MatrixXd& new_density_matrix_alpha, const Eigen::MatrixXd& new_density_matrix_beta) {
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(SpinResolved1DM<Scalar>(SpinResolved1DMComponent<Scalar>(new_density_matrix_alpha), SpinResolved1DMComponent<Scalar>(new_density_matrix_beta)));
             })

        .def("replace_current_density_matrices",
             [](UHFSCFEnvironment<Scalar>& environment, const SpinResolved1DM<Scalar>& D) {
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(D);
             })

        .def("replace_current_fock_matrices",
             [](UHFSCFEnvironment<Scalar>& environment, const Eigen::MatrixXd& new_fock_matrix_alpha, const Eigen::MatrixXd& new_fock_matrix_beta) {
                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(ScalarUSQOneElectronOperator<Scalar>(new_fock_matrix_alpha, new_fock_matrix_beta));
             })

        .def("replace_current_error_vectors",
             [](UHFSCFEnvironment<Scalar>& environment, const Eigen::MatrixXd& new_error_vector_alpha, const Eigen::MatrixXd& new_error_vector_beta) {
                 environment.error_vectors.pop_back();
                 environment.error_vectors.push_back(SpinResolved<VectorX<Scalar>>(new_error_vector_alpha, new_error_vector_beta));
             });
}


/**
 *  Add Python bindings for UHF SCF environments.
 */
void bindUHFSCFEnvironments(py::module& module) {

    // Provide bindings for real-valued UHF SCF environments.
    py::class_<UHFSCFEnvironment<double>> py_UHFSCFEnvironment_d {module, "UHFSCFEnvironment_d", "An algorithm environment that can be used with real-valued UHF SCF solvers."};

    bindQCMethodUHFSCFEnvironmentInterface(py_UHFSCFEnvironment_d);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_UHFSCFEnvironment_d);


    // Provide bindings for complex-valued UHF SCF environments.
    py::class_<UHFSCFEnvironment<complex>> py_UHFSCFEnvironment_cd {module, "UHFSCFEnvironment_cd", "An algorithm environment that can be used with complex-valued UHF SCF solvers."};

    py_UHFSCFEnvironment_cd
        .def_static(
            "WithComplexlyTransformedCoreGuess",
            [](const size_t N_alpha, const size_t N_beta, const USQHamiltonian<complex>& sq_hamiltonian, const ScalarUSQOneElectronOperator<complex>& S) {
                return UHFSCFEnvironment<complex>::WithComplexlyTransformedCoreGuess(N_alpha, N_beta, sq_hamiltonian, S);
            },
            "Initialize an UHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently adding/subtracting a small complex value from certain elements.");

    bindQCMethodUHFSCFEnvironmentInterface(py_UHFSCFEnvironment_cd);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_UHFSCFEnvironment_cd);
}


}  // namespace gqcpy
