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

#include "QCMethod/HF/GHF/GHFSCFEnvironment.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to GHF SCF environments.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `GHFSCFEnvironment`.
 */
template <typename Class>
void bindQCMethodGHFSCFEnvironmentInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::Scalar;


    py_class

        /*
         *  MARK: Named constructors specific to GHF
         */

        .def_static(
            "WithCoreGuess",
            [](const size_t N, const GSQHamiltonian<Scalar>& hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S) {
                return GHFSCFEnvironment<Scalar>::WithCoreGuess(N, hamiltonian, S);
            },
            "Initialize an GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.")

        .def_static(
            "WithTransformedCoreGuess",
            [](const size_t N, const GSQHamiltonian<Scalar>& hamiltonian, const ScalarGSQOneElectronOperator<Scalar>& S, const std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>&)>& transformation_function) {
                return GHFSCFEnvironment<Scalar>::WithTransformedCoreGuess(N, hamiltonian, S, transformation_function);
            },
            "Initialize an GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently applying the given unary transformation function.")


        /*
         *  MARK: Read-write members & properties specific to GHF
         */

        .def_property(
            "S",
            [](const GHFSCFEnvironment<Scalar>& environment) {
                return environment.S;
            },
            [](GHFSCFEnvironment<Scalar>& environment, ScalarGSQOneElectronOperator<Scalar>& S) {
                environment.S = S;
            })

        /*
         *  MARK: Replacements specific to GHF
         */

        .def("replace_current_coefficient_matrix",
             [](GHFSCFEnvironment<Scalar>& environment, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& new_coefficient_matrix) {
                 environment.coefficient_matrices.pop_back();
                 environment.coefficient_matrices.push_back(SquareMatrix<Scalar>(new_coefficient_matrix));
             })

        .def("replace_current_density_matrix",
             [](GHFSCFEnvironment<Scalar>& environment, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& new_density_matrix) {
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(SquareMatrix<Scalar>(new_density_matrix));
             })

        .def("replace_current_fock_matrix",
             [](GHFSCFEnvironment<Scalar>& environment, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& new_fock_matrix) {
                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(SquareMatrix<Scalar>(new_fock_matrix));
             })

        .def("replace_current_error_vectors",
             [](GHFSCFEnvironment<Scalar>& environment, const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& new_error_vectors) {
                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(SquareMatrix<Scalar>(new_error_vectors));
             });
}


/**
 *  Add Python bindings for GHF SCF environments.
 */
void bindGHFSCFEnvironments(py::module& module) {

    // Provide bindings for real-valued GHF SCF environments.
    py::class_<GHFSCFEnvironment<double>> py_GHFSCFEnvironment_d {module, "GHFSCFEnvironment_d", "An algorithm environment that can be used with real-valued GHF SCF solvers."};

    bindQCMethodGHFSCFEnvironmentInterface(py_GHFSCFEnvironment_d);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_GHFSCFEnvironment_d);


    // Provide bindings for complex-valued GHF SCF environments.
    py::class_<GHFSCFEnvironment<complex>> py_GHFSCFEnvironment_cd {module, "GHFSCFEnvironment_cd", "An algorithm environment that can be used with complex-valued GHF SCF solvers."};

    py_GHFSCFEnvironment_cd
        .def_static(
            "WithComplexlyTransformedCoreGuess",
            [](const size_t N, const GSQHamiltonian<complex>& hamiltonian, const ScalarGSQOneElectronOperator<complex>& S) {
                return GHFSCFEnvironment<complex>::WithComplexlyTransformedCoreGuess(N, hamiltonian, S);
            },
            "Initialize an GHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently adding/subtracting a small complex value from certain elements.");

    bindQCMethodGHFSCFEnvironmentInterface(py_GHFSCFEnvironment_cd);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_GHFSCFEnvironment_cd);
}


}  // namespace gqcpy
