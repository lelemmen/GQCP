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
#pragma once


#include "Mathematical/Representation/SquareMatrix.hpp"

#include <Eigen/Sparse>


namespace GQCP {


/**
 *  A templated private class for ONV basiss: its constructors are private, but ONV basiss are its friends
 *  It supports dense, sparse and matrix vector product storage for elements evaluated in the ONV basis.
 *
 *  @tparam _Matrix              the type of matrix in which the evaluations of the ONV basis will be stored
 */
template <typename _Matrix>
class EvaluationIterator {
public:
    using Matrix = _Matrix;

private:
    size_t index = 0;  // current position of the iterator in the dimension of the ONV basis
    size_t end;  // the dimension of the ONV basis

    Matrix matrix;  // matrix containing the evaluations

    /*
     *  CONSTRUCTORS
     */

    /**
     * @param dimension         the dimension of the ONV basis
     */
    EvaluationIterator(const size_t dimension) :  matrix(Matrix::Zero(dimension, dimension)), end(dimension) {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the row and the given index corresponds to the column
     * 
     *  @param column    column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addColumnwise(const size_t column, const double value) {
        this->matrix(this->index, column) += value;
    }

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the column and the given index corresponds to the row
     * 
     *  @param row       row index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addRowwise(const size_t row, const double value) {
        this->matrix(row, this->index) += value;
    }

    /**
     *  @return the evaluation that is stored
     */ 
    const Matrix& evaluation() const { return this->matrix; }

    /**
     *  Move to the next index in the iteration
     */
    void increment() {
        this->index++;
    }

    /**
     *  Tests if the iteration is finished, if true the index is reset to 0 
     * 
     *  @return true if the iteration is finished
     */
    bool is_finished() {
        if (index == end) {
            this->index = 0;
            return true;
        } else {
            return false;
        }
    }


    // Friend classes
    friend class ONVBasis;
    friend class SelectedONVBasis;
    friend class ProductONVBasis;
};



/**
 *  Sparse template specialization is required because insertions into an existing sparse matrix are expensive
 *  Elements should only be added to the matrix once all of them are evaluated in a vector of triplets
 */
template<>
class EvaluationIterator<Eigen::SparseMatrix<double>> {
    size_t index = 0;  // current position of the iterator in the dimension of the ONV basis
    size_t end;  // total dimension

    Eigen::SparseMatrix<double> matrix;  // matrix containing the evaluations
    std::vector<Eigen::Triplet<double>> triplet_vector;  // vector which temporarily contains the added values

    /*
     *  CONSTRUCTORS
     */
    
    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationIterator(const size_t dimension) : 
        matrix(Eigen::SparseMatrix<double>(dimension, dimension)), end(dimension) 
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Reserves an amount of memory for the triplet vector
     *
     *  @param n        the amount triplets that should be reserved
     */
    void reserve(const size_t n) {
        this->triplet_vector.reserve(n);
        this->matrix.reserve(n);
    }


    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the row and the given index corresponds to the column
     *  This function adds the values to a triplet vector
     *  to add the values to the sparse matrix, one should call "addToMatrix()"
     * 
     *  @param column    column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addColumnwise(const size_t column, const double value) {
        this->triplet_vector.emplace_back(this->index, column, value);
    }


    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the column and the given index corresponds to the row
     *  This function adds the values to a triplet vector
     *  to add the values to the sparse matrix, one should call "addToMatrix()"
     * 
     *  @param row       row index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addRowwise(const size_t row, const double value) {
        this->triplet_vector.emplace_back(row, this->index, value);
    }


    /**
     *  Fill the sparse matrix with the elements stored in the triplet vector
     *  the more elements that are already present in the matrix, the more expensive this operation becomes
     *  Therefore, it is ill-advised to call the method more than once,
     *  after filling, the triplet vector is cleared
     */
    void addToMatrix() {
        this->matrix.setFromTriplets(this->triplet_vector.begin(), this->triplet_vector.end());
        this->triplet_vector = {};
    }

    /**
     *  @return the evaluation that is stored
     * 
     *  @note the matrix will not be initialized if `addToMatrix()` has not been called
     */ 
    const Eigen::SparseMatrix<double>& evaluation() const { return this->matrix; }

    /**
     *  Move to the next index in the iteration
     */
    void increment() {
        this->index++;
    }

    /**
     *  Tests if the iteration is finished, if true the index is reset to 0 
     * 
     *  @return true if the iteration is finished
     */
    bool is_finished() {
        if (index == end) {
            this->index = 0;
            return true;
        } else {
            return false;
        }
    }

    /**
     *  @return the triplet vector
     */ 
    const std::vector<Eigen::Triplet<double>>& triplets() const { return triplet_vector; }

    // Friend classes
    friend class ONVBasis;
    friend class SelectedONVBasis;
    friend class ProductONVBasis;
};



/**
 *  Vector template specialization is required because of matvec evaluations are stored in a vector additions
 */
template<>
class EvaluationIterator<VectorX<double>> {
    size_t index = 0;  // current position of the iterator in the dimension of the ONV basis
    size_t end;  // total dimension

    VectorX<double> matvec;  // matvec containing the evaluations
    const VectorX<double>& coefficient_vector;  // vector with which is multiplied
    double sequential_double = 0;  // double which temporarily contains the sum of added values, it corresponds to the value in the matvec of the current index and allows a for a single write access each iteration instead of multiple
    double nonsequential_double = 0;  // double gathered from the coefficient for non-sequential matvec additions, this corresponds to the value of the coefficient vector of the current index, it allows for a single read operation each iteration


    /*
     *  CONSTRUCTORS
     */

    /**
     * @param dimension         the dimensions of the matrix (equal to that of the fock space)
     */
    EvaluationIterator(const VectorX<double>& coefficient_vector, const VectorX<double>& diagonal) : 
        end (coefficient_vector.rows()),
        coefficient_vector(coefficient_vector), 
        matvec(diagonal.cwiseProduct(coefficient_vector))
    {}

    EvaluationIterator(const VectorX<double>& coefficient_vector) : 
        coefficient_vector(coefficient_vector), 
        matvec(VectorX<double>::Zero(coefficient_vector.rows()))
    {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the row and the given index corresponds to the column
     * 
     *  @param column    column index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addColumnwise(const size_t column, const double value) {
        sequential_double += value * coefficient_vector(column);
    }

    /**
     *  Add a value to the matrix evaluation in which the current iterator index corresponds to the column and the given index corresponds to the row
     * 
     *  @param row       row index of the matrix
     *  @param value     the value which is added to a given position in the matrix
     */
    void addRowwise(const size_t row, const double value) {
        this->matvec(row) += value * this->nonsequential_double;
    }

    /**
     *  @return the evaluation that is stored
     */ 
    const VectorX<double>& evaluation() const { return this->matvec; }

    /**
     *  Move to the next index in the iteration, this is accompanied by an addition to the matvec and reset of the sequential double
     */
    void increment() {
        this->matvec(this->index) += sequential_double;
        this->sequential_double = 0;
        this->index++;
    }

    /**
     *  Tests if the iteration is finished, if true the index is reset to 0 
     *  If false the nonsequential_double is updated to the value of the current iteration
     * 
     *  @return true if the iteration is finished
     */
    bool is_finished() {
        if (index == end) {
            this->index = 0;
            return true;
        } else {
            this->nonsequential_double = coefficient_vector(this->index);
            return false;
        }
    }


    // Friend classes
    friend class ONVBasis;
    friend class SelectedONVBasis;
    friend class ProductONVBasis;
};

}  // namespace GQCP
