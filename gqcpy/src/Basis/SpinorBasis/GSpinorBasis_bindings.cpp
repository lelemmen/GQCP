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

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "gqcpy/include/interfaces.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


namespace py = pybind11;


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `GSpinorBasis` over `GTOShell`s.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `GSpinorBasis` over `GTOShell`s.
 */
template <typename Class>
void bindGTOGSpinorBasisInterface(Class& py_class) {

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


    /*
     *  MARK: Named constructors
     */

    py_class
        .def_static("FromRestricted",
                    [](const RSpinOrbitalBasis<Scalar, GTOShell>& spin_orbital_basis) {
                        return Type::FromRestricted(spin_orbital_basis);
                    })

        .def_static("FromUnrestricted",
                    [](const USpinOrbitalBasis<Scalar, GTOShell>& spin_orbital_basis) {
                        return Type::FromUnrestricted(spin_orbital_basis);
                    });


    // Expose some Mulliken API to the Python class;
    bindSpinorBasisMullikenInterface(py_class);
}


/**
 *  Add Python bindings for some APIs related to `GSpinorBasis`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `GSpinorBasis`.
 */
template <typename Class>
void bindGSpinorBasisInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using Scalar = typename Type::ExpansionScalar;


    /*
     *  MARK: Coefficients
     */

    py_class
        .def(
            "numberOfCoefficients",
            &Type::numberOfCoefficients,
            py::arg("sigma"),
            "Return the number of coefficients that are used for the expansion of the requested spin-component of a spinor");


    // Expose the `SimpleSpinorBasis` API to the Python class.
    bindSpinorBasisInterface(py_class);


    /*
     *  MARK: General info
     */

    py_class
        .def(
            "numberOfSpinors",
            &Type::numberOfSpinors,
            "The number of spinors that are described by this generalized spinor basis.");


    // Expose some quantization API to the Python class;
    bindSpinorBasisQuantizationInterface(py_class);
}


/**
 *  Add Python bindings for some APIs related to complex-valued `GSpinorBasis`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to complex-valued `GSpinorBasis`.
 */
template <typename Class>
void bindComplexGSpinorBasisInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;

    py_class
        .def(
            "quantize",
            [](const Type& spinor_basis, const ElectronicSpinOperator& op) {
                return spinor_basis.quantize(op);
            },
            "Return the electronic spin operator expressed in this spinor basis.")


        .def(
            "quantize",
            [](const Type& spinor_basis, const ElectronicSpinSquaredOperator& op) {
                return spinor_basis.quantize(op);
            },
            "Return the second-quantized representation of the electronic S^2 operator in this spinor basis.");
}


/**
 *  Register `GSpinorBasis_d`, `GSpinorBasis_cd` and `London_GSpinorBasis` to the gqcpy module and expose parts of their C++ interfaces to Python.
 * 
 *  @param module           The Pybind11 module in which the classes should be registered.
 */
void bindGSpinorBases(py::module& module) {

    // Define the Python class for `GSpinorBasis_d`.
    py::class_<GSpinorBasis<double, GTOShell>> py_GSpinorBasis_d {module, "GSpinorBasis_d", "A class that represents a real, (generalized) spinor basis with underlying GTO shells."};

    bindGSpinorBasisInterface(py_GSpinorBasis_d);
    bindGTOGSpinorBasisInterface(py_GSpinorBasis_d);


    // Define the Python class for `GSpinorBasis_cd`.
    py::class_<GSpinorBasis<complex, GTOShell>> py_GSpinorBasis_cd {module, "GSpinorBasis_cd", "A class that represents a complex, (generalized) spinor basis with underlying GTO shells."};

    bindGSpinorBasisInterface(py_GSpinorBasis_cd);
    bindGTOGSpinorBasisInterface(py_GSpinorBasis_cd);
    bindComplexGSpinorBasisInterface(py_GSpinorBasis_cd);

    // Define the Python class for `LondonGSpinorBasis`.
    py::class_<GSpinorBasis<complex, LondonGTOShell>> py_LondonGSpinorBasis {module, "LondonGSpinorBasis", "A class that represents a complex, (generalized) spinor basis with underlying London GTO shells."};

    py_LondonGSpinorBasis
        .def(
            "quantize",
            [](const GSpinorBasis<complex, LondonGTOShell>& spinor_basis, const FQMolecularMagneticHamiltonian& op) {
                return spinor_basis.quantize(op);
            },
            "Return the molecular magnetic Hamiltonian expressed in this spinor basis.")

        .def(
            "quantize",
            [](const GSpinorBasis<complex, LondonGTOShell>& spinor_basis, const FQMolecularPauliHamiltonian& op) {
                return spinor_basis.quantize(op);
            },
            "Return the molecular Pauli Hamiltonian expressed in this spinor basis.");

    bindGSpinorBasisInterface(py_LondonGSpinorBasis);
    bindComplexGSpinorBasisInterface(py_LondonGSpinorBasis);
    bindLondonSpinorBasisQuantizationInterface(py_LondonGSpinorBasis);
}


}  // namespace gqcpy
