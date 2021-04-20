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

#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "ONVBasis/SpinResolvedONVBasis.hpp"
#include "ONVBasis/SpinResolvedSelectedONVBasis.hpp"
#include "QCModel/CI/LinearExpansion.hpp"
#include "Utilities/complex.hpp"
#include "gqcpy/include/utilities.hpp"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


// Due to the different kinds of `LinearExpansion`s we can encounter, we opt for an interface approach in order to avoid duplicate code.

/**
 *  Add Python bindings for some APIs related to the linear expansion of a CI calculation.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `LinearExpansion`.
 */
template <typename Class>
void bindQCModelCILinearExpansionInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using ONVBasis = typename Type::ONVBasis;

    py_class

        /*
         *  MARK: Constructors
         */

        .def(py::init([](const ONVBasis& onv_basis, const Eigen::MatrixXd& coefficients) {
                 return Type(onv_basis, coefficients);
             }),
             py::arg("onv_basis"),
             py::arg("coefficients"),
             "Return a linear expansion inside the given ONV basis, with corresponding expansion coefficients.")

        .def_static(
            "Constant",
            [](const ONVBasis& onv_basis) {
                return Type::Constant(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a normalized coefficient vector (i.e. all the coefficients are equal).")

        .def_static(
            "HartreeFock",
            [](const ONVBasis& onv_basis) {
                return Type::HartreeFock(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion that represents the Hartree-Fock wave function.")

        .def_static(
            "Random",
            [](const ONVBasis& onv_basis) {
                return Type::Random(onv_basis);
            },
            py::arg("onv_basis"),
            "Return a linear expansion with a random, normalized coefficient vector, with coefficients uniformly distributed in [-1, +1] before any normalization.")

        .def_static(
            "Normalized",
            [](const ONVBasis& onv_basis, const Eigen::MatrixXd& coefficients) {
                return Type::Normalized(onv_basis, coefficients);
            },
            py::arg("onv_basis"),
            py::arg("coefficients"),
            "Return a normalized linear expansion inside a given ONV basis with possibly non-normalized coefficients.")


        /*
         * MARK: Acces
         */

        .def("coefficients",
             &Type::coefficients,
             "Return the expansion coefficients of this linear expansion wave function model.")

        .def("onvBasis",
             &Type::onvBasis,
             "Return the ONV basis that is related to this linear expansion wave function model.");
}


/**
 *  Add Python bindings for some APIs related to density matrix calculation in a linear expansion of a CI calculation.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `LinearExpansion` density matrix calculation.
 */
template <typename Class>
void bindQCModelCILinearExpansionDensityMatrixInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;
    using ONVBasis = typename Type::ONVBasis;

    py_class

        /*
         *  MARK: Spin unresolved
         */

        .def(
            "calculate1DM",
            [](const Type& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the one-electron density matrix (1-DM) for a wave function expansion of the specified type.")

        .def(
            "calculate2DM",
            [](const Type& linear_expansion) {
                return linear_expansion.calculate2DM();
            },
            "Return the two-electron density matrix (2-DM) for a wave function expansion of the specified type.")


        /*
         * MARK: Spin resolved
         */

        .def(
            "calculateSpinResolved1DM",
            [](const Type& linear_expansion) {
                return linear_expansion.calculateSpinResolved1DM();
            },
            "Return the spin-resolved 1-DM for a wave function expansion of the specified type.")

        .def(
            "calculateSpinResolved2DM",
            [](const Type& linear_expansion) {
                return linear_expansion.calculateSpinResolved2DM();
            },
            "Return the spin resolved 2-DM for a wave function expansion of the specified type.");
}


void bindLinearExpansions(py::module& module) {

    // Define the python class related to a linear expansion of a Seniority Zero ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<double, SeniorityZeroONVBasis>> py_LinearExpansion_SeniorityZero {module, "LinearExpansion_SeniorityZero", "The linear expansion (configuration interaction) wave function model in a seniority-zero ONV basis."};

    // Expose the bindings exclusive to this type of linear expansion.
    py_LinearExpansion_SeniorityZero

        /*
         * MARK: Entropy
         */
        .def(
            "calculateShannonEntropy",
            [](const LinearExpansion<double, SeniorityZeroONVBasis>& linear_expansion) {
                return linear_expansion.calculateShannonEntropy();
            },
            "Return the Shannon entropy (information content) of the wave function.")

        .def(
            "calculateSingleOrbitalEntropy",
            [](const LinearExpansion<double, SeniorityZeroONVBasis>& linear_expansion, const size_t orbital_index) {
                return linear_expansion.calculateSingleOrbitalEntropy(orbital_index);
            },
            py::arg("orbital_index"),
            "Return the single orbital entropy of the orbital at the specified index, according to the formula of Boguslawski (https://doi.org/10.1002/qua.24832).");

    // Expose the linear expansion interface.
    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SeniorityZero);

    // Expose the linear expansion density matrix interface.
    bindQCModelCILinearExpansionDensityMatrixInterface(py_LinearExpansion_SeniorityZero);

    // Define the python class related to a linear expansion of a Spin Resolved ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<double, SpinResolvedONVBasis>> py_LinearExpansion_SpinResolved {module, "LinearExpansion_SpinResolved", "The linear expansion (configuration interaction) wave function model in a spin-resolved ONV basis."};

    // Expose the bindings exclusive to this type of linear expansion.
    py_LinearExpansion_SpinResolved

        /*
         * MARK: Constructors
         */

        .def_static(
            "FromONVProjection",
            [](const SpinResolvedONV& onv, const RSpinOrbitalBasis<double, GTOShell>& r_spinor_basis, const USpinOrbitalBasis<double, GTOShell>& u_spinor_basis) {
                return LinearExpansion<double, SpinResolvedONVBasis>::FromONVProjection(onv, r_spinor_basis, u_spinor_basis);
            },
            py::arg("onv"),
            py::arg("r_spinor_basis"),
            py::arg("u_spinor_basis"),
            "Create the linear expansion of the given spin-resolved ONV that is expressed in the given USpinOrbitalBasis, by projection onto the spin-resolved ONVs expressed with respect to the given RSpinOrbitalBasis.")


        /*
         * MARK: Basis transformations
         */

        .def(
            "basisTransform",
            [](LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion, const Eigen::MatrixXd& T) {
                return linear_expansion.basisTransform(RTransformation<double>(T));
            },
            py::arg("T"),
            "Update the expansion coefficients of this linear expansion so that they correspond to the situation after a transformation of the underlying spinor basis with the given basis transformation.")


        /*
         * MARK: Iteration
         */

        .def(
            "forEach",
            [](const LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion, const std::function<void(const double, const SpinResolvedONV)>& callback) {
                return linear_expansion.forEach(callback);
            },
            py::arg("callback"),
            "Iterate over all expansion coefficients and corresponding ONVs, and apply the given callback function.")


        /*
         * MARK: Entropy
         */

        .def(
            "calculateShannonEntropy",
            [](const LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculateShannonEntropy();
            },
            "Return the Shannon entropy (information content) of the wave function.")

        .def(
            "calculateSingleOrbitalEntropy",
            [](const LinearExpansion<double, SpinResolvedONVBasis>& linear_expansion, const size_t orbital_index) {
                return linear_expansion.calculateSingleOrbitalEntropy(orbital_index);
            },
            py::arg("orbital_index"),
            "Return the single orbital entropy of the orbital at the specified index, according to the formula of Boguslawski (https://doi.org/10.1002/qua.24832).");

    // Expose the linear expansion interface.
    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SpinResolved);

    // Expose the linear expansion density matrix interface.
    bindQCModelCILinearExpansionDensityMatrixInterface(py_LinearExpansion_SpinResolved);

    // Define the python class related to a linear expansion of a Spin Resolved Selected ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<double, SpinResolvedSelectedONVBasis>> py_LinearExpansion_SpinResolvedSelected {module, "LinearExpansion_SpinResolvedSelected", "The linear expansion (configuration interaction) wave function model in a spin-resolved selected ONV basis."};

    // Expose the bindings exclusive to this type of linear expansion
    py_LinearExpansion_SpinResolvedSelected

        /*
         * MARK: Constructors
         */

        .def_static(
            "FromGAMESSUS",
            [](const std::string& GAMESSUS_filename) {
                return LinearExpansion<double, SpinResolvedSelectedONVBasis>::FromGAMESSUS(GAMESSUS_filename);
            },
            py::arg("GAMESSUS_filename"),
            "Return the corresponding spin-resolved selected linear expansion from a given GAMESS-US file.")

        /*
         *  MARK: Entropy
         */

        .def(
            "calculateShannonEntropy",
            [](const LinearExpansion<double, SpinResolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculateShannonEntropy();
            },
            "Return the Shannon entropy (information content) of the wave function.");


    // Expose the linear expansion interface.
    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SpinResolvedSelected);

    // Expose the linear expansion density matrix interface.
    bindQCModelCILinearExpansionDensityMatrixInterface(py_LinearExpansion_SpinResolvedSelected);


    // Define the python class related to a linear expansion of a Spin Unresolved ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<double, SpinUnresolvedONVBasis>> py_LinearExpansion_SpinUnresolved {module, "LinearExpansion_SpinUnresolved", "The linear expansion (configuration interaction) wave function model in a spin-unresolved ONV basis."};

    // Expose the bindings exclusive to this type of linear expansion
    py_LinearExpansion_SpinUnresolved

        /*
         * MARK: Constructors
         */

        .def_static(
            "FromONVProjection",
            [](const SpinUnresolvedONV& onv_of, const GSpinorBasis<double, GTOShell>& spinor_basis_on, const GSpinorBasis<double, GTOShell>& spinor_basis_of) {
                return LinearExpansion<double, SpinUnresolvedONVBasis>::FromONVProjection(onv_of, spinor_basis_on, spinor_basis_of);
            },
            py::arg("onv_of"),
            py::arg("spinor_basis_on"),
            py::arg("spinor_basis_of"),
            "Create the linear expansion of the given spin-unresolved ONV that is expressed in the given GSpinorBasis, by projection onto the spin-resolved ONVs expressed with respect to another given GSpinorBasis.")

        /*
         * MARK: Density matrices
         */

        .def(
            "calculate1DM",
            [](const LinearExpansion<double, SpinUnresolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the generalized one-electron density matrix.")

        .def(
            "calculateNDMElement",
            [](const LinearExpansion<double, SpinUnresolvedONVBasis>& linear_expansion, const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) {
                return linear_expansion.calculateNDMElement(bra_indices, ket_indices);
            },
            py::arg("bra_indices"),
            py::arg("ket_indices"),
            "Return an element of the N-DM, as specified by the given bra and ket indices. `calculateNDMElement({0, 1}, {2, 1})` would calculate an element of the 2-NDM d^{(2)} (0, 1, 1, 2) corresponding the operator string: `a^dagger_0 a^dagger_1 a_2 a_1`.")

        /*
         *  MARK: Entropy
         */

        .def(
            "calculateShannonEntropy",
            [](const LinearExpansion<double, SpinUnresolvedONVBasis>& linear_expansion) {
                return linear_expansion.calculateShannonEntropy();
            },
            "Return the Shannon entropy (information content) of the wave function.");

    // Expose the linear expansion interface.
    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SpinUnresolved);


    // Define the python class related to a linear expansion of a spin unresolved selected ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<double, SpinUnresolvedSelectedONVBasis>> py_LinearExpansion_SpinUnresolvedSelected_d {module, "LinearExpansion_SpinResolvedSelected_d", "The real-valued linear expansion (configuration interaction) wave function model in a spin-unresolved selected ONV basis."};

    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SpinUnresolvedSelected_d);

    py_LinearExpansion_SpinUnresolvedSelected_d

        /*
         *  MARK: Density matrices
         */

        .def(
            "calculate1DM",
            [](const LinearExpansion<double, SpinUnresolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the generalized (G) 1-DM.")

        .def(
            "calculate2DM",
            [](const LinearExpansion<double, SpinUnresolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculate2DM();
            },
            "Return the generalized (G) 2-DM.")

        .def(
            "calculateNDMElement",
            [](const LinearExpansion<double, SpinUnresolvedSelectedONVBasis>& linear_expansion, const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) {
                return linear_expansion.calculateNDMElement(bra_indices, ket_indices);
            },
            py::arg("bra_indices"),
            py::arg("ket_indices"),
            "Return an element of the N-DM, as specified by the given bra and ket indices. `calculateNDMElement({0, 1}, {2, 1})` would calculate an element of the 2-NDM d^{(2)} (0, 1, 1, 2) corresponding the operator string: `a^dagger_0 a^dagger_1 a_2 a_1`.")

        /*
         *  MARK: Entropy
         */

        .def(
            "calculateShannonEntropy",
            [](const LinearExpansion<double, SpinUnresolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculateShannonEntropy();
            },
            "Return the Shannon entropy (information content) of the wave function.");


    // Define the python class related to a linear expansion of a spin unresolved selected ONV basis and expose the necessary interfaces.
    py::class_<LinearExpansion<complex, SpinUnresolvedSelectedONVBasis>> py_LinearExpansion_SpinUnresolvedSelected_cd {module, "LinearExpansion_SpinResolvedSelected_cd", "The complex-valued linear expansion (configuration interaction) wave function model in a spin-unresolved selected ONV basis."};

    bindQCModelCILinearExpansionInterface(py_LinearExpansion_SpinUnresolvedSelected_cd);

    py_LinearExpansion_SpinUnresolvedSelected_cd

        /*
         *  MARK: Density matrices
         */

        .def(
            "calculate1DM",
            [](const LinearExpansion<complex, SpinUnresolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculate1DM();
            },
            "Return the generalized (G) 1-DM.")

        .def(
            "calculate2DM",
            [](const LinearExpansion<complex, SpinUnresolvedSelectedONVBasis>& linear_expansion) {
                return linear_expansion.calculate2DM();
            },
            "Return the generalized (G) 2-DM.")

        .def(
            "calculateNDMElement",
            [](const LinearExpansion<double, SpinUnresolvedSelectedONVBasis>& linear_expansion, const std::vector<size_t>& bra_indices, const std::vector<size_t>& ket_indices) {
                return linear_expansion.calculateNDMElement(bra_indices, ket_indices);
            },
            py::arg("bra_indices"),
            py::arg("ket_indices"),
            "Return an element of the N-DM, as specified by the given bra and ket indices. `calculateNDMElement({0, 1}, {2, 1})` would calculate an element of the 2-NDM d^{(2)} (0, 1, 1, 2) corresponding the operator string: `a^dagger_0 a^dagger_1 a_2 a_1`.");
}


}  // namespace gqcpy
