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

#include "Basis/ScalarBasis/GTOShell.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Utilities/aliases.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `RSpinOrbitalBasis`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `RSpinOrbitalBasis`.
 */
template <typename Class>
void bindRSpinOrbitalBasisInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::ExpansionScalar;


    /*
     *  MARK: Constructors
     */

    py_class
        .def(py::init<const Molecule&, const std::string&>(),
             py::arg("molecule"),
             py::arg("basisset_name"));


    // Expose the `SimpleSpinorBasis` API to the Python class.
    bindSpinorBasisInterface(py_class);

    /*
     *  MARK: General information
     */

    py_class
        .def(
            "numberOfSpatialOrbitals",
            [](const RSpinOrbitalBasis<Scalar, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.numberOfSpatialOrbitals();
            },
            "Return the number of different spatial orbitals that are used in this spinor basis.")


        .def(
            "scalarBasis",
            [](const RSpinOrbitalBasis<Scalar, GTOShell>& spin_orbital_basis) {
                return spin_orbital_basis.scalarBasis();
            },
            "Return the underlying scalar basis with respect to which the basis coefficients are expressed.");


    /*
     *  MARK: Quantization of first-quantized operators
     */

    // Expose some quantization API to the Python class;
    bindSpinorBasisQuantizationInterface(py_class);

    // Expose some Mulliken API to the Python class;
    bindSpinorBasisMullikenInterface(py_class);
}


/**
 *  Register `RSpinOrbitalBasis_d` and `RSpinOrbitalBasis_cd` to the gqcpy module and expose parts of their C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindRSpinOrbitalBases(py::module& module) {

    // Define the Python class for `RSpinOrbitalBasis_d`.
    py::class_<RSpinOrbitalBasis<double, GTOShell>> py_RSpinOrbitalBasis_d {module, "RSpinOrbitalBasis_d", "A class that represents a real, restricted spinor basis with underlying GTO shells."};
    bindRSpinOrbitalBasisInterface(py_RSpinOrbitalBasis_d);

    py_RSpinOrbitalBasis_d
        .def(
            "quantize",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const ElectronicDipoleOperator& op) {
                return spin_orbital_basis.quantize(op);
            },
            "Return the electronic dipole operator expressed in this spinor basis.")

        .def(
            "quantize",
            [](const RSpinOrbitalBasis<double, GTOShell>& spin_orbital_basis, const ElectronicDensityOperator& op) {
                return spin_orbital_basis.quantize(op);
            },
            "Return the electronic density operator expressed in this spinor basis.");


    // Define the Python class for `RSpinOrbitalBasis_cd`.
    py::class_<RSpinOrbitalBasis<complex, GTOShell>> py_RSpinOrbitalBasis_cd {module, "RSpinOrbitalBasis_cd", "A class that represents a complex, restricted spinor basis with underlying GTO shells."};
    bindRSpinOrbitalBasisInterface(py_RSpinOrbitalBasis_cd);

    py_RSpinOrbitalBasis_cd
        .def(
            "quantize",
            [](const RSpinOrbitalBasis<complex, GTOShell>& spin_orbital_basis, const AngularMomentumOperator& op) {
                return spin_orbital_basis.quantize(op);
            },
            "Return the angular momentum operator expressed in this spinor basis.")

        .def(
            "quantize",
            [](const RSpinOrbitalBasis<complex, GTOShell>& spin_orbital_basis, const LinearMomentumOperator& op) {
                return spin_orbital_basis.quantize(op);
            },
            "Return the linear momentum operator expressed in this spinor basis.")

        .def(
            "quantize",
            [](const RSpinOrbitalBasis<complex, GTOShell>& spin_orbital_basis, const CurrentDensityOperator& op) {
                return spin_orbital_basis.quantize(op);
            },
            "Return the current density operator expressed in this spinor basis.");
}


}  // namespace gqcpy
