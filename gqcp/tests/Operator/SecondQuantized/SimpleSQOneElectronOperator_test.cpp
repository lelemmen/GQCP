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

#define BOOST_TEST_MODULE "SimpleSQOneElectronOperator"

#include <boost/test/unit_test.hpp>

#include "Basis/SpinorBasis/GSpinorBasis.hpp"
#include "Basis/SpinorBasis/RSpinOrbitalBasis.hpp"
#include "Operator/SecondQuantized/GSQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/RSQOneElectronOperator.hpp"
#include "Physical/units.hpp"
#include "QCModel/HF/GHF.hpp"
#include "QCModel/HF/RHF.hpp"

#include <boost/math/constants/constants.hpp>


/**
 *  In this test suite, we'll test the behavior of `SimpleSQOneElectronOperator` through the derived classes `RSQOneElectronOperator` and `GSQOneElectronOperator`.
 */

/**
 *  Check the construction of one-electron operators from matrices.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_constructor) {

    // Set up some test matrices to be used as one-electron operator integrals.
    const auto square_matrix3 = GQCP::SquareMatrix<double>::Zero(3);
    const auto square_matrix4 = GQCP::SquareMatrix<double>::Zero(4);


    // Check a correct constructor: a square matrix may represent the integrals of a one-electron operator.
    BOOST_CHECK_NO_THROW(GQCP::ScalarRSQOneElectronOperator<double> operator_square {square_matrix3});

    // Check a throwing constructor: the dimensions of the square matrices must be equal.
    std::vector<GQCP::SquareMatrix<double>> matrices {square_matrix3, square_matrix3, square_matrix4};
    BOOST_CHECK_THROW(GQCP::VectorRSQOneElectronOperator<double> operator_rectangular {matrices}, std::invalid_argument);
}


/**
 *  Check if the `Zero` named constructor actually sets its parameters to zeros.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_zero) {

    const size_t dim = 2;

    // Check a zero constructor for scalar operators.
    const auto zero_scalar_operator = GQCP::ScalarRSQOneElectronOperator<double>::Zero(2);

    BOOST_CHECK_EQUAL(zero_scalar_operator.numberOfOrbitals(), dim);
    BOOST_CHECK(zero_scalar_operator.parameters().isZero(1.0e-08));


    // Check a zero constructor for vector operators.
    const auto zero_vector_operator = GQCP::VectorRSQOneElectronOperator<double>::Zero(2);

    BOOST_CHECK_EQUAL(zero_vector_operator.numberOfOrbitals(), dim);
    for (size_t i = 0; i < 3; i++) {
        BOOST_CHECK(zero_vector_operator.parameters(i).isZero(1.0e-08));
    }
}


/**
 *  Check if we can access a single component of a SimpleSQOneElectronOperator easily using the operator().
 */
BOOST_AUTO_TEST_CASE(operator_call) {

    const size_t dim = 4;

    // Set up three (random) matrix representations and construct the corresponding vector operator.
    const GQCP::SquareMatrix<double> M1 = GQCP::SquareMatrix<double>::Random(dim);
    const GQCP::SquareMatrix<double> M2 = GQCP::SquareMatrix<double>::Random(dim);
    const GQCP::SquareMatrix<double> M3 = GQCP::SquareMatrix<double>::Random(dim);

    const std::vector<GQCP::SquareMatrix<double>> matrices {M1, M2, M3};
    const GQCP::VectorRSQOneElectronOperator<double> vector_operator {matrices};

    // Check if we can access a single component of this operator.
    BOOST_CHECK(vector_operator(0).parameters().isApprox(M1, 1.0e-12));
    BOOST_CHECK(vector_operator(1).parameters().isApprox(M2, 1.0e-12));
    BOOST_CHECK(vector_operator(2).parameters().isApprox(M3, 1.0e-12));
}


/**
 *  Check if calculateExpectationValue throws when necessary.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_throw) {

    // Initialize a test operator and density matrices.
    const GQCP::ScalarRSQOneElectronOperator<double> f_operator {2};
    const GQCP::Orbital1DM<double> D_valid {GQCP::SquareMatrix<double>::Zero(2)};
    const GQCP::Orbital1DM<double> D_invalid {GQCP::SquareMatrix<double>::Zero(3)};

    BOOST_CHECK_NO_THROW(f_operator.calculateExpectationValue(D_valid));
    BOOST_CHECK_THROW(f_operator.calculateExpectationValue(D_invalid), std::invalid_argument);  // Invalid dimension for the density matrix.
}


/**
 *  Check if calculateExpectationValue shows the correct behaviour.
 */
BOOST_AUTO_TEST_CASE(calculateExpectationValue_behaviour) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op {M1};

    // Initialize a test density matrix.
    GQCP::SquareMatrix<double> D {dim};
    // clang-format off
    D << 0.0, 1.0,
         1.0, 0.0;
    // clang-format on
    const GQCP::Orbital1DM<double> DM {D};

    // Initialize a reference value and check the result.
    const double reference_expectation_value = 5.0;

    const double expectation_value = op.calculateExpectationValue(DM);  // A scalar-StorageArray can be implicitly converted to the underlying scalar.
    BOOST_CHECK(std::abs(expectation_value - reference_expectation_value) < 1.0e-08);
}


/**
 *  Check if the addition of operators works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_addition) {

    const size_t dim = 2;

    // Initialize two test matrices and convert them into operators.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};
    // clang-format on

    GQCP::SquareMatrix<double> M2 {dim};
    // clang-format off
    M2 << 5.0, 6.0,
          7.0, 8.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op2 {M2};


    // Initialize the reference and check the result of the operator addition.
    GQCP::SquareMatrix<double> M_sum_ref {dim};
    // clang-format off
    M_sum_ref <<  6.0,  8.0,
                 10.0, 12.0;
    // clang-format on

    const auto op_sum = op1 + op2;
    BOOST_CHECK(op_sum.parameters().isApprox(M_sum_ref, 1.0e-08));
}


/**
 *  Check if scalar multiplication of a one-electron operator works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_scalar_product) {

    const size_t dim = 2;
    const double scalar = 2.0;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};

    // Initialize the reference and check the result of the scalar multiplication.
    GQCP::SquareMatrix<double> M_ref {dim};
    // clang-format off
    M_ref <<  2.0,  4.0,
              6.0,  8.0;
    // clang-format on

    const auto op_result1 = scalar * op1;
    BOOST_CHECK(op_result1.parameters().isApprox(M_ref, 1.0e-08));

    const auto op_result2 = op1 * scalar;
    BOOST_CHECK(op_result2.parameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if the negation of a one-electron operator works as expected.
 */
BOOST_AUTO_TEST_CASE(SimpleSQOneElectronOperator_negation) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op1 {M1};

    // Initialize the reference and check the result of the negation.
    GQCP::SquareMatrix<double> M_ref {dim};
    // clang-format off
    M_ref <<  -1.0,  -2.0,
              -3.0,  -4.0;
    // clang-format on

    const auto op_result = -op1;
    BOOST_CHECK(op_result.parameters().isApprox(M_ref, 1.0e-08));
}


/**
 *  Check if dot product multiplication of a SimpleSQOneElectronOperator with a Vector with the same components is correctly implemented.
 */
BOOST_AUTO_TEST_CASE(dot) {

    // Set up a toy RSQOneElectronOperator with three components.
    GQCP::SquareMatrix<double> h_x {2};
    // clang-format off
    h_x << 1.0, 2.0,
           3.0, 4.0;
    // clang-format on

    GQCP::SquareMatrix<double> h_y {2};
    // clang-format off
    h_y << -1.0,  2.0,
            3.0, -4.0;
    // clang-format on

    GQCP::SquareMatrix<double> h_z {2};
    // clang-format off
    h_z <<  0.0,  3.0,
            3.0,  0.0;
    // clang-format on

    const std::vector<GQCP::SquareMatrix<double>> matrices {h_x, h_y, h_z};
    GQCP::VectorRSQOneElectronOperator<double> h_op {matrices};


    // Define a vector to do the dot product with, and check the result.
    GQCP::VectorX<double> a {3};
    a << 1.0, 2.0, 3.0;

    GQCP::SquareMatrix<double> dot_ref {2};
    // clang-format off
    dot_ref <<  -1.0,  15.0,
                18.0,  -4.0;
    // clang-format on

    BOOST_CHECK(h_op.dot(a).parameters().isApprox(dot_ref, 1.0e-12));
}


/**
 *  Check if the basis transformation works for a trivial case.
 */
BOOST_AUTO_TEST_CASE(transform_trivial) {

    const size_t dim = 3;

    // Initialize a test matrix and convert it into an operator.
    const GQCP::SquareMatrix<double> f = GQCP::SquareMatrix<double>::Random(dim);
    GQCP::ScalarRSQOneElectronOperator<double> op {f};

    // Initialize a trivial transformation, i.e. the identity transformation.
    const auto T = GQCP::RTransformation<double>::Identity(dim);

    // Check the in-place and the returning methods.
    const auto op_transformed = op.transformed(T);
    BOOST_CHECK(op_transformed.parameters().isApprox(f, 1.0e-12));

    op.transform(T);
    BOOST_CHECK(op.parameters().isApprox(f, 1.0e-12));
}


/**
 *  Check if the basis transformation method works as expected, for a non-trivial case.
 */
BOOST_AUTO_TEST_CASE(transform) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    GQCP::ScalarRSQOneElectronOperator<double> op {M1};

    // Initialize a corresponding transformation.
    GQCP::SquareMatrix<double> T_matrix {dim};
    // clang-format off
    T_matrix << 2.0, 3.0,
                4.0, 5.0;
    // clang-format on
    GQCP::RTransformation<double> T {T_matrix};

    // Initialize the reference matrix, corresponding to the result of the basis transformation.
    GQCP::SquareMatrix<double> ref {dim};
    // clang-format off
    ref << 108.0, 142.0,
           140.0, 184.0;
    // clang-format on

    // Check the in-place and 'returning' methods.
    const auto op_transformed = op.transformed(T);
    BOOST_CHECK(op_transformed.parameters().isApprox(ref, 1.0e-08));

    op.transform(T);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 *  Check if, if we transform with a transformation and its inverse, we get a zero operation.
 */
BOOST_AUTO_TEST_CASE(transform_and_inverse) {

    const size_t dim = 3;

    // Initialize a a transformation matrix and its inverse.
    GQCP::SquareMatrix<double> T_matrix {dim};
    // clang-format off
    T_matrix << 1,  0,  0,
                0, -2,  0,
                0,  0,  3;
    // clang-format on
    const GQCP::RTransformation<double> T {T_matrix};
    const auto T_inverse = T.inverse();

    // Initialize a random one-electron operator.
    const GQCP::SquareMatrix<double> f = GQCP::SquareMatrix<double>::Random(dim);
    GQCP::ScalarRSQOneElectronOperator<double> op {f};

    // Transform the one-electron operator with T and its inverse, and check if it is a zero operation.
    op.transformed(T).transformed(T_inverse);
    BOOST_CHECK(op.parameters().isApprox(f, 1.0e-12));
}


/**
 *  Check if the basis rotation method works as expected.
 */
BOOST_AUTO_TEST_CASE(rotate) {

    const size_t dim = 2;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> M1 {dim};
    // clang-format off
    M1 << 1.0, 2.0,
          3.0, 4.0;
    // clang-format on
    GQCP::ScalarRSQOneElectronOperator<double> op {M1};

    // Initialize a unitary transformation.
    GQCP::SquareMatrix<double> U_matrix {dim};
    // clang-format off
    U_matrix << 1.0,  0.0,
                0.0, -1.0;
    // clang-format on
    const GQCP::RTransformation<double> U {U_matrix};

    // Initialize the reference matrix, corresponding to the result of the basis transformation.
    GQCP::SquareMatrix<double> ref {dim};
    // clang-format off
    ref <<  1.0, -2.0,
           -3.0,  4.0;
    // clang-format on


    // Check the in-place and 'returning' methods.
    const auto op_rotated = op.rotated(U);
    BOOST_CHECK(op_rotated.parameters().isApprox(ref, 1.0e-08));

    op.rotate(U);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 *  Check a rotation test case.
 */
BOOST_AUTO_TEST_CASE(rotate_test_case) {

    const size_t dim = 3;


    // Initialize a test case one-electron operator.
    GQCP::SquareMatrix<double> f1 {dim};
    // clang-format off
    f1 << 1.0, 0.5, 0.0,
          0.5, 2.0, 0.0,
          0.0, 0.0, 1.0;
    // clang-format on
    const GQCP::ScalarRSQOneElectronOperator<double> op {f1};


    // Initialize a test unitary transformation.
    GQCP::SquareMatrix<double> U_matrix {dim};
    // clang-format off
    U_matrix << 1.0,  0.0,  0.0,
                0.0,  0.0, -1.0,
                0.0, -1.0,  0.0;
    // clang-format on
    const GQCP::RTransformation<double> U {U_matrix};


    // Initialize the reference, and check the result.
    GQCP::SquareMatrix<double> f_ref {dim};
    // clang-format off
    f_ref <<  1.0, 0.0, -0.5,
              0.0, 1.0,  0.0,
             -0.5, 0.0,  2.0;
    // clang-format on

    BOOST_CHECK(op.rotated(U).parameters().isApprox(f_ref, 1.0e-08));
}


/**
 *  Check that we can't use the `rotate` function with a non-unitary matrix.
 */
BOOST_AUTO_TEST_CASE(rotate_throws) {

    // Create a random one-electron operator.
    const size_t dim = 3;
    const GQCP::SquareMatrix<double> f = GQCP::SquareMatrix<double>::Random(dim);
    GQCP::ScalarRSQOneElectronOperator<double> op {f};


    // Check if rotating with a non-unitary transformation causes a throw.
    const GQCP::RTransformation<double> T = GQCP::RTransformation<double>::Random(dim);  // The probability of a random matrix being unitary approaches zero.
    BOOST_CHECK_THROW(op.rotate(T), std::invalid_argument);


    // Check if rotating with a unitary transformation is accepted.
    const GQCP::RTransformation<double> U = GQCP::RTransformation<double>::RandomUnitary(dim);
    BOOST_CHECK_NO_THROW(op.rotate(U));
}


/**
 *  Check if rotating with JacobiRotation is the same as with the corresponding Jacobi rotation matrix.
 */
BOOST_AUTO_TEST_CASE(rotate_jacobi_vs_matrix) {

    // Create a random one-electron operator.
    const size_t dim = 5;
    const GQCP::SquareMatrix<double> f = GQCP::SquareMatrix<double>::Random(dim);
    GQCP::ScalarRSQOneElectronOperator<double> op {f};

    // Create Jacobi rotation and the corresponding transformation.
    GQCP::JacobiRotation jacobi_rotation {4, 2, 56.81};
    const auto J = GQCP::RTransformation<double>::FromJacobi(jacobi_rotation, dim);


    // Rotate using both representations and check the result.
    BOOST_CHECK(op.rotated(jacobi_rotation).parameters().isApprox(op.rotated(J).parameters(), 1.0e-12));
}


/**
 *  Check if the Jacobi rotations are correctly applied, for a 3-dimensional case.
 */
BOOST_AUTO_TEST_CASE(jacobi_rotation_3) {

    const size_t dim = 3;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> f {dim};
    // clang-format off
    f << 1.0, 0.5, 0.0,
         0.5, 2.0, 0.0,
         0.0, 0.0, 1.0;
    // clang-format on
    GQCP::ScalarRSQOneElectronOperator<double> op {f};


    // Initialize the Jacobi rotation.
    const GQCP::JacobiRotation J {1, 0, boost::math::constants::half_pi<double>()};  // interchanges two orbitals and applies a sign change


    // Initialize the reference result.
    GQCP::SquareMatrix<double> ref {dim};
    // clang-format off
    ref <<  2.0, -0.5, 0.0,
           -0.5,  1.0, 0.0,
            0.0,  0.0, 1.0;
    // clang-format on


    // Check the in-place and 'returning' methods.
    const auto op_rotated = op.rotated(J);
    BOOST_CHECK(op_rotated.parameters().isApprox(ref, 1.0e-08));

    op.rotate(J);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 *  Check if the Jacobi rotations are correctly applied, for a 4-dimensional case.
 */
BOOST_AUTO_TEST_CASE(jacobi_rotation_4) {

    const size_t dim = 4;

    // Initialize a test matrix and convert it into an operator.
    GQCP::SquareMatrix<double> f {dim};
    // clang-format off
    f <<  1.0,  2.0,  3.0,  4.0,
          5.0,  6.0,  7.0,  8.0,
          9.0, 10.0, 11.0, 12.0,
         13.0, 14.0, 15.0, 16.0;
    // clang-format on
    GQCP::ScalarRSQOneElectronOperator<double> op {f};

    // Initialize the Jacobi rotation.
    const GQCP::JacobiRotation J {2, 1, boost::math::constants::half_pi<double>()};  // interchanges two orbitals and applies a sign change

    // Initialize the reference result.
    GQCP::SquareMatrix<double> ref {dim};
    // clang-format off
    ref <<   1.0,  3.0,  -2.0,  4.0,
             9.0, 11.0, -10.0, 12.0,
            -5.0, -7.0,   6.0, -8.0,
            13.0, 15.0, -14.0, 16.0;
    // clang-format on

    // Check the in-place and 'returning' methods.
    const auto op_rotated = op.rotated(J);
    BOOST_CHECK(op_rotated.parameters().isApprox(ref, 1.0e-08));

    op.rotate(J);
    BOOST_CHECK(op.parameters().isApprox(ref, 1.0e-08));
}


/**
 *  Check that the total RHF Mulliken population (i.e. Mulliken-partitioned number operator expectation value) of N2 is 14, for a restricted partitioning.
 */
BOOST_AUTO_TEST_CASE(mulliken_N2_STO_3G_restricted) {

    // Initialize the molecular Hamiltonian for N2 in the Löwdin-orthonormalized basis.
    const GQCP::Nucleus N1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{N1, N2}};

    GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell> spin_orbital_basis {molecule, "STO-3G"};
    spin_orbital_basis.lowdinOrthonormalize();


    // Set up the Mulliken partitioning scheme.
    using Shell = GQCP::RSpinOrbitalBasis<double, GQCP::GTOShell>::Shell;
    const auto mulliken_partitioning = spin_orbital_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "N";
        });

    // The Mulliken operator is the Mulliken-partitioned number (i.e. overlap) operator.
    const auto S = spin_orbital_basis.overlap();
    const auto mulliken_op = S.partitioned(mulliken_partitioning);


    // Create the RHF 1-DM for N2 and check the total Mulliken operator.
    const auto K = spin_orbital_basis.numberOfSpatialOrbitals();
    const auto N = molecule.numberOfElectrons();
    const auto D = GQCP::QCModel::RHF<double>::calculateOrthonormalBasis1DM(K, N);

    const double mulliken_population = mulliken_op.calculateExpectationValue(D);  // A scalar-StorageArray can be implicitly converted into the underlying scalar.
    BOOST_CHECK(std::abs(mulliken_population - N) < 1.0e-12);
}


/**
 *  Check that the total RHF Mulliken population (i.e. Mulliken-partitioned number operator expectation value) of N2 is 14, for a restricted partitioning.
 */
BOOST_AUTO_TEST_CASE(mulliken_N2_STO_3G_generalized) {

    // Initialize the molecular Hamiltonian for N2 in the Löwdin-orthonormalized basis.
    const GQCP::Nucleus N1 {7, 0.0, 0.0, 0.0};
    const GQCP::Nucleus N2 {7, 0.0, 0.0, GQCP::units::angstrom_to_bohr(1.134)};  // From CCCBDB, STO-3G geometry.
    const GQCP::Molecule molecule {{N1, N2}};

    GQCP::GSpinorBasis<double, GQCP::GTOShell> spinor_basis {molecule, "STO-3G"};
    spinor_basis.lowdinOrthonormalize();


    // Set up the Mulliken partitioning scheme.
    using Shell = GQCP::GSpinorBasis<double, GQCP::GTOShell>::Shell;
    const auto mulliken_partitioning = spinor_basis.mullikenPartitioning(
        [](const Shell& shell) {
            return shell.nucleus().element() == "N";
        });

    // The Mulliken operator is the Mulliken-partitioned number (i.e. overlap) operator.
    const auto S = spinor_basis.overlap();
    const auto mulliken_op = S.partitioned(mulliken_partitioning);


    // Create the GHF 1-DM for N2 and check the total Mulliken operator.
    const auto M = spinor_basis.numberOfSpinors();
    const auto N = molecule.numberOfElectrons();
    const auto D = GQCP::QCModel::GHF<double>::calculateOrthonormalBasis1DM(M, N);

    const double mulliken_population = mulliken_op.calculateExpectationValue(D);  // A scalar-StorageArray can be implicitly converted into the underlying scalar.
    BOOST_CHECK(std::abs(mulliken_population - N) < 1.0e-12);
}
