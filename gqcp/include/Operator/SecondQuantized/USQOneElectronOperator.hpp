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

#pragma once


#include "Basis/Transformations/JacobiRotationParameters.hpp"
#include "Basis/Transformations/TransformationMatrix.hpp"
#include "DensityMatrix/OneDM.hpp"
#include "DensityMatrix/SpinResolvedOneDM.hpp"
#include "DensityMatrix/TwoDM.hpp"
#include "Mathematical/Functions/ScalarFunction.hpp"
#include "Mathematical/Representation/QCMatrix.hpp"
#include "QuantumChemical/Spin.hpp"
#include "Utilities/type_traits.hpp"

#include <array>


namespace GQCP {


/**
 *  A class that represents an 'unrestricted second-quantized one-electron operator' suitable for the projection of the non-relativistic Hamiltonian onto an unrestricted spinor basis. It holds the matrix representation of its parameters for both spin components, which are (usually) integrals over first-quantized operators
 *
 *  @tparam _Scalar             the scalar type, i.e. the scalar representation of one of the parameters
 *  @tparam _Components         the number of components of the second-quantized operator
 */
template <typename _Scalar, size_t _Components>
class USQOneElectronOperator {
public:
    using Scalar = _Scalar;
    static constexpr auto Components = _Components;


private:
    std::array<QCMatrix<Scalar>, Components> fs_a;  // all the matrix representations of the alpha spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
    std::array<QCMatrix<Scalar>, Components> fs_b;  // all the matrix representations of the beta spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param fs_a    all the matrix representations of the alpha spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     *  @param fs_b    all the matrix representations of the beta spin component (hence the s) of the parameters (integrals) of the different components of this second-quantized operator
     */
    USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>& fs_a, const std::array<QCMatrix<Scalar>, Components>& fs_b) :
        fs_a {fs_a},
        fs_b {fs_b} {

        // Check if the given matrix representations have the same dimensions.
        const auto dimension_of_first_a = this->fs_a[0].numberOfOrbitals();
        const auto dimension_of_first_b = this->fs_b[0].numberOfOrbitals();

        for (size_t i = 1; i < Components; i++) {

            const auto dimension_of_ith_a = this->fs_a[i].numberOfOrbitals();
            const auto dimension_of_ith_b = this->fs_b[i].numberOfOrbitals();

            if ((dimension_of_first_a != dimension_of_ith_a) || (dimension_of_first_b != dimension_of_ith_b)) {
                throw std::invalid_argument("USQOneElectronOperator(const std::array<QCMatrix<Scalar>, Components>&, const std::array<QCMatrix<Scalar>, components>&): The given matrix representations do not have the same dimensions for either the alpha or beta component.");
            }
        }
    }


    /**
     *  A constructor for ScalarUSQOneElectronOperators that doesn't require the argument to be a vector of just one element.
     * 
     *  @param f_a            the matrix representation of the integrals of this scalar second-quantized operator, for the alpha spin component
     *  @param f_b            the matrix representation of the integrals of this scalar second-quantized operator, for the beta spin component
     * 
     *  @note This constructor is only available for ScalarUSQOneElectronOperators (for the std::enable_if, see https://stackoverflow.com/a/17842695/7930415)
     */
    template <size_t Z = Components>
    USQOneElectronOperator(const QCMatrix<Scalar>& f_a, const QCMatrix<Scalar>& f_b, typename std::enable_if<Z == 1>::type* = 0) :
        USQOneElectronOperator(std::array<QCMatrix<Scalar>, 1> {f_a}, std::array<QCMatrix<Scalar>, 1> {f_b}) {}


    /**
     *  Construct a one-electron operator with parameters that are zero, for both the alpha and beta spin components
     * 
     *  @param dim          the dimension of the matrix representation of the parameters, i.e. the number of orbitals/sites. This dimension is the same for the alpha and beta component.
     */
    USQOneElectronOperator(const size_t dim) {
        for (size_t i = 0; i < Components; i++) {
            this->fs_a[i] = QCMatrix<Scalar>::Zero(dim);
            this->fs_b[i] = QCMatrix<Scalar>::Zero(dim);
        }
    }


    /**
     *  Default constructor: construct a one-electron operator with parameters that are zero, for both spin components.
     */
    USQOneElectronOperator() :
        USQOneElectronOperator(0)  // dimensions of the representations are zero
    {}


    /*
     *  MARK: Accessing parameters
     */

    /**
     *  @param sigma                The requested spin component. This can be either alpha or beta.
     * 
     *  @return read-only matrix representations of all the parameters (integrals) of the different components of this second-quantized operator, for the requested spin component.
     */
    const std::array<QCMatrix<Scalar>, Components>& allParameters(const Spin sigma) const {

        if (sigma == Spin::alpha) {
            return this->fs_a;
        } else {
            return this->fs_b;
        };
    }


    /**
     *  @param sigma                The requested spin component. This can be either alpha or beta.
     * 
     *  @return writable matrix representations of all the parameters (integrals) of the different components of this second-quantized operator, for the requested spin component.
     */
    std::array<QCMatrix<Scalar>, Components>& allParameters(const Spin sigma) {

        if (sigma == Spin::alpha) {
            return this->fs_a;
        } else {
            return this->fs_b;
        };
    }


    /**
     *  @param i            the index
     * 
     *  @return the i-th component (i.e. the i-th alpha and i-th beta component) of this operator
     */
    USQOneElectronOperator<Scalar, 1> operator[](const size_t i) const {

        if (i >= Components) {
            throw std::invalid_argument("USQOneElectronOperator::operator[](const size_t): The given index is out of bounds.");
        }

        return USQOneElectronOperator<Scalar, 1> {this->fs_a[i], this->fs_b[i]};
    }


    /**
     *  @param sigma                    The requested spin component. This can be either alpha or beta.
     *  @param i                        The index of the component.
     * 
     *  @return a read-only the matrix representation of the parameters (integrals) of one of the the different components of this second-quantized operator, for the requested spin component.
     */
    const QCMatrix<Scalar>& parameters(const Spin sigma, const size_t i = 0) const {

        if (sigma == Spin::alpha) {
            return this->fs_a[i];
        } else {
            return this->fs_b[i];
        };
    }


    /**
     *  @param sigma                    The requested spin component. This can be either alpha or beta.
     *  @param i                        The index of the component.
     * 
     *  @return the writable matrix representation of the parameters (integrals) of one of the components of this second-quantized operator, for the requested spin component.
     */
    QCMatrix<Scalar>& parameters(const Spin sigma, const size_t i = 0) {

        if (sigma == Spin::alpha) {
            return this->fs_a[i];
        } else {
            return this->fs_b[i];
        };
    }


    /*
     *  MARK: Calculations
     */

    /**
     *  @param D                the spin-resolved 1-DM that represents the wave function
     *
     *  @return the expectation values of all components of the one-electron operator
     */
    Vector<Scalar, Components> calculateExpectationValue(const SpinResolvedOneDM<Scalar>& D) const {

        if (this->fs_a[0].numberOfOrbitals() != D.numberOfOrbitals(Spin::alpha) || this->fs_b[0].numberOfOrbitals() != D.numberOfOrbitals(Spin::beta)) {
            throw std::invalid_argument("USQOneElectronOperator::calculateExpectationValue(const OneDM<Scalar>&, const OneDM<Scalar>&): The given 1-DM is not compatible with the one-electron operator.");
        }

        std::array<Scalar, Components> expectation_values {};  // zero initialization

        for (size_t i = 0; i < Components; i++) {
            expectation_values[i] = (this->parameters(GQCP::Spin::alpha, i) * D.alpha()).trace() + (this->parameters(GQCP::Spin::beta, i) * D.beta()).trace();
        }

        return Eigen::Map<Eigen::Matrix<Scalar, Components, 1>>(expectation_values.data());  // convert std::array to Vector
    }


    /*
     *  MARK: General information
     */

    /**
     *  @return the sum of the alpha and beta dimensions
     */
    size_t numberOfOrbitals() const {
        return this->numberOfOrbitals(GQCP::Spin::alpha) + this->numberOfOrbitals(GQCP::Spin::beta);
    }


    /**
     *  @param sigma                The requested spin component. This can be either alpha or beta.
     * 
     *  @return the dimension of the matrices for the requested spin component.
     */
    size_t numberOfOrbitals(const Spin sigma) const {

        if (sigma == Spin::alpha) {
            return this->fs_a[0].numberOfOrbitals();
        } else {
            return this->fs_b[0].numberOfOrbitals();
        };
    }


    /*
     *  MARK: Basis transformations
     */

    /**
     *  In-place rotate the operator to another basis. The same matrix is used for the alpha and beta components.   
     * 
     *  @param U                            the (unitary) rotation matrix
     */
    void rotate(const TransformationMatrix<Scalar>& U) {

        // Transform the matrix representations of the components
        for (auto& f_a : this->allParameters(GQCP::Spin::alpha)) {
            f_a.basisRotate(U);
        }
        for (auto& f_b : this->allParameters(GQCP::Spin::beta)) {
            f_b.basisRotate(U);
        }
    }


    /**
     *  In-place rotate the operator using a unitary Jacobi rotation matrix constructed from the Jacobi rotation parameters. The same matrix is used to rotate the alpha and beta components.
     * 
     *  @param jacobi_rotation_parameters       the Jacobi rotation parameters (p, q, angle) that are used to specify a Jacobi rotation: we use the (cos, sin, -sin, cos) definition for the Jacobi rotation matrix
     */
    void rotate(const JacobiRotationParameters& jacobi_rotation_parameters) {

        // Transform the matrix representations of the components
        for (auto& f_a : this->allParameters(GQCP::Spin::alpha)) {
            f_a.basisRotate(jacobi_rotation_parameters);
        }
        for (auto& f_b : this->allParameters(GQCP::Spin::beta)) {
            f_b.basisRotate(jacobi_rotation_parameters);
        }
    }


    /**
     *  In-place transform the operator to another basis. The same transformation is applied on the alpha and beta components.
     * 
     *  @param T                            the transformation matrix
     */
    void transform(const TransformationMatrix<Scalar>& T) {

        // Transform the matrix representations of the components
        for (auto& f_a : this->allParameters(GQCP::Spin::alpha)) {
            f_a.basisTransform(T);
        }
        for (auto& f_b : this->allParameters(GQCP::Spin::beta)) {
            f_b.basisTransform(T);
        }
    }
};


/*
 *  MARK: Aliases
 */
template <typename Scalar>
using ScalarUSQOneElectronOperator = USQOneElectronOperator<Scalar, 1>;

template <typename Scalar>
using VectorUSQOneElectronOperator = USQOneElectronOperator<Scalar, 3>;


/*
 *  MARK: Vector-like arithmetic
 */

/**
 *  Add two one-electron operators by adding their parameters. The two alphas are added together and the two betas are added together. 
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator+(const USQOneElectronOperator<LHSScalar, Components>& lhs, const USQOneElectronOperator<RHSScalar, Components>& rhs) -> USQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto F_sum_a = lhs.allParameters(GQCP::Spin::alpha);
    auto F_sum_b = lhs.allParameters(GQCP::Spin::beta);
    for (size_t i = 0; i < Components; i++) {
        F_sum_a[i] += rhs.parameters(GQCP::Spin::alpha, i);
        F_sum_b[i] += rhs.parameters(GQCP::Spin::beta, i);
    }

    return USQOneElectronOperator<ResultScalar, Components>(F_sum_a, F_sum_b);
}


/**
 *  Multiply a one-electron operator with a scalar. The alpha and beta components are multiplied with the same scalar.
 * 
 *  @tparam Scalar              the scalar type of the scalar
 *  @tparam OperatorScalar      the scalar type of the operator
 * 
 *  @tparam scalar              the scalar of the scalar multiplication
 *  @tparam op                  the one-electron operator
 */
template <typename Scalar, typename OperatorScalar, size_t Components>
auto operator*(const Scalar& scalar, const USQOneElectronOperator<OperatorScalar, Components>& op) -> USQOneElectronOperator<product_t<Scalar, OperatorScalar>, Components> {

    using ResultScalar = product_t<Scalar, OperatorScalar>;

    auto fs_a = op.allParameters(GQCP::Spin::alpha);
    auto fs_b = op.allParameters(GQCP::Spin::beta);
    for (auto& f_a : fs_a) {
        f_a *= scalar;
    }
    for (auto& f_b : fs_b) {
        f_b *= scalar;
    }

    return USQOneElectronOperator<ResultScalar, Components>(fs_a, fs_b);
}


/**
 *  Negate a one-electron operator
 * 
 *  @tparam Scalar              the scalar type of the operator
 *  @tparam Components          the number of components of the one-electron operator
 * 
 *  @param op                   the operator
 */
template <typename Scalar, size_t Components>
USQOneElectronOperator<Scalar, Components> operator-(const USQOneElectronOperator<Scalar, Components>& op) {

    return (-1.0) * op;  // negation is scalar multiplication with (-1.0)
}


/**
 *  Subtract two one-electron operators by subtracting their parameters. The alphas are subtracted from the alphas and the betas from the betas.
 * 
 *  @tparam LHSScalar           the scalar type of the left-hand side
 *  @tparam RHSScalar           the scalar type of the right-hand side
 *  @tparam Components          the number of components of the one-electron operators
 * 
 *  @param lhs                  the left-hand side
 *  @param rhs                  the right-hand side
 */
template <typename LHSScalar, typename RHSScalar, size_t Components>
auto operator-(const USQOneElectronOperator<LHSScalar, Components>& lhs, const USQOneElectronOperator<RHSScalar, Components>& rhs) -> USQOneElectronOperator<sum_t<LHSScalar, RHSScalar>, Components> {

    using ResultScalar = sum_t<LHSScalar, RHSScalar>;

    auto F_min_a = lhs.allParameters(GQCP::Spin::alpha);
    auto F_min_b = lhs.allParameters(GQCP::Spin::beta);
    for (size_t i = 0; i < Components; i++) {
        F_min_a[i] -= rhs.parameters(GQCP::Spin::alpha, i);
        F_min_b[i] -= rhs.parameters(GQCP::Spin::beta, i);
    }

    return USQOneElectronOperator<ResultScalar, Components>(F_min_a, F_min_b);
}


}  // namespace GQCP
