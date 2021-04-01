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

#include "Mathematical/Grid/Field.hpp"
#include "Utilities/aliases.hpp"

#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <functional>


namespace gqcpy {


// Provide some shortcuts for frequent namespaces.
namespace py = pybind11;
using namespace GQCP;


/**
 *  Add Python bindings for some APIs related to `Field`.
 * 
 *  @tparam Class               The type of the Pybind11 `class_` (generated by the compiler).
 *  
 *  @param py_class             The Pybind11 `class_` that should obtain APIs related to `Field`.
 */
template <typename Class>
void bindFieldInterface(Class& py_class) {

    // The C++ type corresponding to the Python class.
    using Type = typename Class::type;

    // The type of the evaluated function values.
    using T = typename Type::T;


    py_class

        /*
         *  MARK: Constructors
         */

        .def(py::init<const std::vector<T>&>(),
             py::arg("values"))


        /*
         *  MARK: Canonical mathematical operators
         */

        .def(py::self += py::self)
        .def(py::self + py::self)
        .def(-py::self)
        .def(py::self -= py::self)
        .def(py::self - py::self)


        /*
         *  MARK: General information
         */

        .def(
            "size",
            &Type::size,
            "Return the size of this field, i.e. the number of field values.")


        /*
         *  MARK: Mappings
         */

        .def(
            "map",
            [](Type& field, const std::function<T(const T&)>& function) {
                field.map(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values, in-place.")

        .def(
            "mapped",
            [](const Type& field, const std::function<T(const T&)>& function) {
                field.mapped(function);
            },
            py::arg("function"),
            "Apply a given function on each of this field's values.")


        /*
         *  MARK: Access
         */

        .def(
            "value",
            [](const Field<T>& field, const size_t index) {
                return field.value(index);
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "value",
            [](Field<T>& field, const size_t index) {
                return field.value(index);
            },
            py::arg("index"),
            "Return a read-only field value, corresponding to the given index.")

        .def(
            "values",
            &Field<T>::values,
            "Return the evaluated function values, in the order of the grid's loop.");
}


/**
 *  Register multiple variants of `Field` to the gqcpy module and expose a part of their C++ interface to Python.
 * 
 *  @param module           The Pybind11 module in which the variants of `Field` should be registered.
 */
void bindField(py::module& module) {

    // Register `py_Field_d` to the Python module.
    py::class_<Field<double>> py_Field_d {module, "Field_d", "A set of real-valued function values corresponding to points in space."};

    py_Field_d
        .def_static(
            "ReadCubeFile",
            [](const std::string& filename) {
                return Field<double>::ReadCubeFile(filename);
            },
            py::arg("filename"),
            "Parse a GAUSSIAN Cube file (http://paulbourke.net/dataformats/cube/) for scalar field values. The grid-associated information is discarded.");


    // Register `py_Field_cd` to the Python module.
    py::class_<Field<complex>> py_Field_cd {module, "Field_cd", "A set of complex-valued function values corresponding to points in space."};


    // Register `py_VectorField_d` to the Python module.
    py::class_<Field<Vector<double, 3>>> py_VectorField_d {module, "VectorField_d", "A set of real-valued 3D-vector function values corresponding to points in space."};

    py_VectorField_d
        .def_static(
            "ReadGridFile",
            [](const std::string& filename) {
                return Field<Vector<double, 3>>::ReadGridFile<3>(filename);
            },
            py::arg("filename"),
            "Parse an .igrid- or .rgrid-file and create the Field that is contained in it. The values for the grid are discarded.");


    // Register `py_VectorField_cd` to the Python module.
    py::class_<Field<Vector<complex, 3>>> py_VectorField_cd {module, "VectorField_cd", "A set of complex-valued 3D-vector function values corresponding to points in space."};


    // Register `py_MatrixField_d` to the Python module.
    py::class_<Field<Matrix<double, 3, 3>>> py_MatrixField_d {module, "MatrixField_d", "A set of real-valued 3D-matrix function values corresponding to points in space."};


    // Register `py_MatrixField_d` to the Python module.
    py::class_<Field<Matrix<complex, 3, 3>>> py_MatrixField_cd {module, "MatrixField_cd", "A set of complex-valued 3D-matrix function values corresponding to points in space."};
}


}  // namespace gqcpy
