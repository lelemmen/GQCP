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

#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "ONVBasis/SpinUnresolvedSelectedONVBasis.hpp"
#include "QCMethod/CI/CI.hpp"
#include "QCMethod/CI/CIEnvironment.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `QCMethod::CI`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `RSpinOrbitalBasis`.
 */
template <typename Class>
void bindQCMethodCIInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;

    using Scalar = typename Type::Scalar;
    using ONVBasis = typename Type::ONVBasis;


    py_class

        /*
         *  MARK: Constructors
         */

        .def(py::init<const ONVBasis, const size_t>(),
             py::arg("onv_basis"),
             py::arg("number_of_states") = 1)


        /*
         *  MARK: Optimization
         */

        .def(
            "optimize",
            [](const QCMethod::CI<Scalar, ONVBasis>& qc_method, Algorithm<EigenproblemEnvironment<Scalar>>& solver, EigenproblemEnvironment<Scalar>& environment) {
                return qc_method.optimize(solver, environment);
            },
            "Optimize the CI wave function model: find the linear expansion coefficients.")

        .def(
            "optimize",
            [](const QCMethod::CI<Scalar, ONVBasis>& qc_method, IterativeAlgorithm<EigenproblemEnvironment<Scalar>>& solver, EigenproblemEnvironment<Scalar>& environment) {
                return qc_method.optimize(solver, environment);
            },
            "Optimize the CI wave function model: find the linear expansion coefficients.");
}


/**
 *  Register multiple templated `QCMethod::CI` to the gqcpy module and expose parts of their C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindQCMethodCI(py::module& module) {

    py::class_<QCMethod::CI<double, SpinResolvedONVBasis>> py_CI_d_SpinResolvedONVBasis {module, "CI_d_SpinResolvedONVBasis", "Real-valued configuration interaction in a spin-resolved ONV basis."};
    bindQCMethodCIInterface(py_CI_d_SpinResolvedONVBasis);

    py::class_<QCMethod::CI<double, SeniorityZeroONVBasis>> py_CI_d_SeniorityZeroONVBasis {module, "CI_d_SeniorityZeroONVBasis", "Real-valued configuration interaction in a seniority zero ONV basis."};
    bindQCMethodCIInterface(py_CI_d_SeniorityZeroONVBasis);

    py::class_<QCMethod::CI<double, SpinResolvedSelectedONVBasis>> py_CI_d_SpinResolvedSelectedONVBasis {module, "CI_d_SpinResolvedSelectedONVBasis", "Real-valued configuration interaction in a spin-resolved selected ONV basis."};
    bindQCMethodCIInterface(py_CI_d_SpinResolvedSelectedONVBasis);


    py::class_<QCMethod::CI<double, SpinUnresolvedONVBasis>> py_CI_d_SpinUnresolvedONVBasis {module, "CI_d_SpinUnresolvedONVBasis", "Real-valued configuration interaction in a spin-unresolved ONV basis."};
    bindQCMethodCIInterface(py_CI_d_SpinUnresolvedONVBasis);

    py::class_<QCMethod::CI<double, SpinUnresolvedSelectedONVBasis>> py_CI_d_SpinUnesolvedSelectedONVBasis {module, "CI_d_SpinUnresolvedSelectedONVBasis", "Real-valued configuration interaction in a spin-unresolved selected ONV basis."};
    bindQCMethodCIInterface(py_CI_d_SpinUnesolvedSelectedONVBasis);

    py::class_<QCMethod::CI<complex, SpinUnresolvedSelectedONVBasis>> py_CI_cd_SpinUnesolvedSelectedONVBasis {module, "CI_cd_SpinUnresolvedSelectedONVBasis", "Complex-valued configuration interaction in a spin-unresolved selected ONV basis."};
    bindQCMethodCIInterface(py_CI_cd_SpinUnesolvedSelectedONVBasis);
}


}  // namespace gqcpy
