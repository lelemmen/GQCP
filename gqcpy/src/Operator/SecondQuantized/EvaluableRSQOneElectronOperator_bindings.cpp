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

#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Operator/SecondQuantized/EvaluableRSQOneElectronOperator.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `EvaluableRSQOneElectronOperatorInterface`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `EvaluableRSQOneElectronOperatorInterface`.
 */
template <typename Class>
void bindEvaluableRSQOneElectronOperatorInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::ExpansionScalar;


    /*
     *  MARK: Evaluations
     */

    py_class
        .def(
            "evaluate",
            [](const Type& sq_op, const Vector<double, 3>& r) {
                return sq_op.evaluate(r);
            },
            py::arg("r"),
            "Evaluate this one-electron operator for the given argument.");
}


/**
 *  Register multiple variants of `RSQDensityOperator` and `RSQCurrentDensityOperator` to the gqcpy module and expose a part of its C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which `RSQDensityOperator` and `RSQCurrentDensityOperator` should be registered.
 */
void bindEvaluableRSQOneElectronOperator(py::module& module) {

    // The type that represents a density distribution for a restricted spin-orbital basis.
    using DensityDistribution = RSpinOrbitalBasis<double, GTOShell>::DensityDistribution;

    // The type that represents a current density distribution for a restricted spin-orbital basis.
    using CurrentDensityDistribution = RSpinOrbitalBasis<double, GTOShell>::CurrentDensityDistribution;


    // Define a Python class `RSQDensityOperator` for the representation of the second-quantized density operator.
    py::class_<ScalarEvaluableRSQOneElectronOperator<DensityDistribution>> py_RSQDensityOperator {module, "RSQDensityOperator", "The second-quantized density operator."};
    bindEvaluableRSQOneElectronOperatorInterface(py_RSQDensityOperator);


    // Define a Python class `RSQCurrentDensityOperator` for the representation of the second-quantized current density operator.
    py::class_<VectorEvaluableRSQOneElectronOperator<CurrentDensityMatrixElement<complex, CartesianGTO>>> py_RSQCurrentDensityOperator {module, "RSQCurrentDensityOperator", "The second-quantized current density operator."};
    bindEvaluableRSQOneElectronOperatorInterface(py_RSQCurrentDensityOperator);
}


}  // namespace gqcpy
