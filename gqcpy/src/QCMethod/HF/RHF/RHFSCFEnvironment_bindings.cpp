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

#include "QCMethod/HF/RHF/RHFSCFEnvironment.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to RHF SCF environments.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `RHFSCFEnvironment`.
 */
template <typename Class>
void bindQCMethodRHFSCFEnvironmentInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::Scalar;


    py_class

        /*
         *  MARK: Named constructors specific to RHF
         */

        .def_static(
            "WithCoreGuess",
            [](const size_t N, const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQOneElectronOperator<Scalar>& S) {
                return RHFSCFEnvironment<Scalar>::WithCoreGuess(N, sq_hamiltonian, S);
            },
            "Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix.")

        .def_static(
            "WithTransformedCoreGuess",
            [](const size_t N, const RSQHamiltonian<Scalar>& sq_hamiltonian, const ScalarRSQOneElectronOperator<Scalar>& S, const std::function<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>(const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>&)>& transformation_function) {
                return RHFSCFEnvironment<Scalar>::WithTransformedCoreGuess(N, sq_hamiltonian, S, transformation_function);
            },
            "Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently applying the given unary transformation function.")

        /*
         *  MARK: Read-write members & properties specific to RHF
         */

        .def_property(
            "S",
            [](const RHFSCFEnvironment<Scalar>& environment) {
                return environment.S;
            },
            [](RHFSCFEnvironment<Scalar>& environment, const ScalarRSQOneElectronOperator<Scalar>& S) {
                environment.S = S;
            })

        /*
         *  MARK: Replacements specific to RHF
         */

        .def("replace_current_coefficient_matrix",
             [](RHFSCFEnvironment<Scalar>& environment, const RTransformation<Scalar>& new_coefficient_matrix) {
                 environment.coefficient_matrices.pop_back();
                 environment.coefficient_matrices.push_back(new_coefficient_matrix);
             })

        .def("replace_current_density_matrix",
             [](RHFSCFEnvironment<Scalar>& environment, const Orbital1DM<Scalar>& new_density_matrix) {
                 environment.density_matrices.pop_back();
                 environment.density_matrices.push_back(new_density_matrix);
             })

        .def("replace_current_fock_matrix",
             [](RHFSCFEnvironment<Scalar>& environment, const ScalarRSQOneElectronOperator<Scalar>& new_fock_matrix) {
                 environment.fock_matrices.pop_back();
                 environment.fock_matrices.push_back(new_fock_matrix);
             })

        .def("replace_current_error_vectors",
             [](RHFSCFEnvironment<Scalar>& environment, const Eigen::MatrixXd& new_error_vectors) {
                 environment.error_vectors.pop_back();
                 environment.error_vectors.push_back(VectorX<Scalar>(new_error_vectors));
             });
}


/**
 *  Add Python bindings for RHF SCF environments.
 */
void bindRHFSCFEnvironments(py::module& module) {

    // Provide bindings for real-valued RHF SCF environments.
    py::class_<RHFSCFEnvironment<double>> py_RHFSCFEnvironment_d {module, "RHFSCFEnvironment_d", "An algorithm environment that can be used with real-valued RHF SCF solvers."};

    bindQCMethodRHFSCFEnvironmentInterface(py_RHFSCFEnvironment_d);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_RHFSCFEnvironment_d);


    // Provide bindings for complex-valued RHF SCF environments.
    py::class_<RHFSCFEnvironment<complex>> py_RHFSCFEnvironment_cd {module, "RHFSCFEnvironment_cd", "An algorithm environment that can be used with complex-valued RHF SCF solvers."};

    py_RHFSCFEnvironment_cd
        .def_static(
            "WithComplexlyTransformedCoreGuess",
            [](const size_t N, const RSQHamiltonian<complex>& sq_hamiltonian, const ScalarRSQOneElectronOperator<complex>& S) {
                return RHFSCFEnvironment<complex>::WithComplexlyTransformedCoreGuess(N, sq_hamiltonian, S);
            },
            "Initialize an RHF SCF environment with an initial coefficient matrix that is obtained by diagonalizing the core Hamiltonian matrix and subsequently adding/subtracting a small complex value from certain elements.");

    bindQCMethodRHFSCFEnvironmentInterface(py_RHFSCFEnvironment_cd);
    bindQCMethodHartreeFockSCFEnvironmentInterface(py_RHFSCFEnvironment_cd);
}


}  // namespace gqcpy
