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

#include "QCMethod/CC/CCDSolver.hpp"
#include "Utilities/complex.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `CCDSolver`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `CCDSolver`.
 */
template <typename Class>
void bindCCDSolverInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;

    using Scalar = typename Type::Scalar;


    py_class
        .def_static(
            "Plain",
            [](const double threshold, const size_t maximum_number_of_iterations) {
                return CCDSolver<Scalar>::Plain(threshold, maximum_number_of_iterations);
            },
            py::arg("threshold") = 1.0e-08,
            py::arg("maximum_number_of_iterations") = 128,
            "Return a plain CCD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion.")

        .def_static(
            "DIIS",
            [](const size_t minimum_subspace_dimension, const size_t maximum_subspace_dimension, const double threshold, const size_t maximum_number_of_iterations) {
                return CCDSolver<Scalar>::DIIS(minimum_subspace_dimension, maximum_subspace_dimension, threshold, maximum_number_of_iterations);
            },
            py::arg("minimum_subspace_dimension") = 6,
            py::arg("maximum_subspace_dimension") = 6,
            py::arg("threshold") = 1.0e-08,
            py::arg("maximum_number_of_iterations") = 128,
            "Return a DIIS CCD solver that uses the norm of the difference of consecutive amplitudes as a convergence criterion.");
}


/**
 *  Register multiple templated `CCDSolver` to the gqcpy module and expose parts of their C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindCCDSolver(py::module& module) {

    // Real-valued Python bindings.
    py::class_<CCDSolver<double>> py_CCDSolver_d {module, "CCDSolver_d", "A factory class that can construct CCD solvers for real-valued calculations."};
    bindCCDSolverInterface(py_CCDSolver_d);

    // Complex-valued Python bindings.
    py::class_<CCDSolver<complex>> py_CCDSolver_cd {module, "CCDSolver_cd", "A factory class that can construct CCD solvers for complex-valued calculations."};
    bindCCDSolverInterface(py_CCDSolver_cd);
}


}  // namespace gqcpy
