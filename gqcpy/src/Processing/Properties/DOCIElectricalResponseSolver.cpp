// This file is part of GQCG-gqcp.
//
// Copyright (C) 2017-2019  the GQCG developers
//
// GQCG-gqcp is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// GQCG-gqcp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with GQCG-gqcp.  If not, see <http://www.gnu.org/licenses/>.
//
#include "ONVBasis/SeniorityZeroONVBasis.hpp"
#include "Processing/Properties/CIElectricalResponseSolver.hpp"
#include "Processing/Properties/properties.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>


namespace py = pybind11;


namespace gqcpy {


// void bindDOCIElectricalResponseSolver(py::module& module) {

//     py::class_<GQCP::CIElectricalResponseSolver<GQCP::SeniorityZeroONVBasis>>(module, "vAP1roGElectricalResponseSolver", "A class whose instances can solve the response equations for DOCI.")

//         .def(py::init<const GQCP::LinearExpansion<GQCP::SeniorityZeroONVBasis>>(),
//             py::arg("linear_expansion")
//         )

//         .def("calculateWaveFunctionResponse",
//             [ ] (const GQCP::CIElectricalResponseSolver<GQCP::SeniorityZeroONVBasis>& response_solver, const GQCP::SQHamiltonian<double>& sq_hamiltonian, const GQCP::VectorSQOneElectronOperator<double>& dipole_op) {
//                 return response_solver.calculateWaveFunctionResponse(sq_hamiltonian, dipole_op);
//             },
//             "Solve the parameter-linear response equations and return the wave function response.",
//             py::arg("sq_hamiltonian"),
//             py::arg("dipole_op")
//         )

//         .def("calculateElectricPolarizability",
//             [ ] (const GQCP::CIElectricalResponseSolver<GQCP::SeniorityZeroONVBasis>& response_solver, const Eigen::Matrix<double, Eigen::Dynamic, 3>& x, const GQCP::VectorSQOneElectronOperator<double>& dipole_op) {
//                 const auto F_p = response_solver.calculateParameterResponseForce(dipole_op);

//                 return GQCP::calculateElectricPolarizability(F_p, x);
//             },
//             "Calculate the DOCI electric polarizability",
//             py::arg("x"),
//             py::arg("dipole_op")
//         )
//     ;
// }


}  // namespace gqcpy
