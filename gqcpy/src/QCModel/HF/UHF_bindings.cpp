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

#include "QCModel/HF/UHF.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


void bindQCModelUHF(py::module& module) {
    py::class_<GQCP::QCModel::UHF<double>>(module, "QCModel_UHF", "The unrestricted Hartree-Fock wave function model.")

        // PUBLIC METHODS

        .def(
            "calculateOrthonormalBasis1DM",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters) {
                return uhf_parameters.calculateOrthonormalBasis1DM();
            },
            "Return the sigma-spin UHF 1-DM expressed in an orthonormal sigma spin-orbital basis for these UHF model parameters.")

        .def(
            "calculateScalarBasis1DM",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters) {
                return uhf_parameters.calculateScalarBasis1DM();
            },
            "Return the sigma-spin UHF 1-DM expressed in the underlying scalar basis for these UHF model parameters.")

        .def(
            "coefficientMatrix",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters, const GQCP::Spin sigma) {
                return uhf_parameters.coefficientMatrix(sigma);
            },
            py::arg("sigma"),
            "Return the coefficient matrix that expresses the sigma spin-orbitals (as a column) in its underlying scalar basis.")

        .def(
            "numberOfElectrons",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters, const GQCP::Spin sigma) {
                return uhf_parameters.numberOfElectrons(sigma);
            },
            py::arg("sigma"),
            "Return the number of sigma electrons that these UHF model parameters describe, i.e. the number of occupied sigma-spin-orbitals.")

        .def(
            "numberOfSpinOrbitals",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters, const GQCP::Spin sigma) {
                return uhf_parameters.numberOfSpinOrbitals(sigma);
            },
            py::arg("sigma"),
            "Return the number of sigma spin-orbitals that these UHF model parameters describe.")

        .def(
            "numberOfSpinOrbitals",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters) {
                return uhf_parameters.numberOfSpinOrbitals();
            },
            "Return the number of spin-orbitals that these UHF model parameters describe.")

        .def(
            "orbitalEnergies",
            [](const GQCP::QCModel::UHF<double>& uhf_parameters, const GQCP::Spin sigma) {
                return uhf_parameters.orbitalEnergies(sigma);
            },
            py::arg("sigma"),
            "Return the orbital energies of the sigma-spin-orbitals.")

        .def(
            "spinOrbitalEnergiesBlocked",
            &GQCP::QCModel::UHF<double>::spinOrbitalEnergiesBlocked,
            "Return all the spin-orbital energies, with the alpha spin-orbital energies appearing before the beta spin-orbital energies");
}


}  // namespace gqcpy
