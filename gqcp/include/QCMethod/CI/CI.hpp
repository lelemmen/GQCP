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


#include "Mathematical/Algorithm/Algorithm.hpp"
#include "Mathematical/Algorithm/IterativeAlgorithm.hpp"
#include "Mathematical/Optimization/Eigenproblem/EigenproblemEnvironment.hpp"
#include "QCMethod/QCStructure.hpp"
#include "QCModel/CI/LinearExpansion.hpp"

#include <memory>


namespace GQCP {
namespace QCMethod {


/**
 *  The configuration interaction quantum chemical method.
 * 
 *  @tparam _Scalar             The scalar type of the expansion coefficients: real or complex.
 *  @tparam _ONVBasis           The type of ONV basis.
 */
template <typename _Scalar, typename _ONVBasis>
class CI {
public:
    // The scalar type of the expansion coefficients: real or complex.
    using Scalar = _Scalar;

    // The type of the ONV basis.
    using ONVBasis = _ONVBasis;


private:
    // The number of states that are searched for (including the ground state).
    size_t number_of_states;

    // The ONV basis with respect to which the configuration interaction is expressed.
    ONVBasis onv_basis;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param onv_basis                        The ONV basis with respect to which the configuration interaction is expressed.
     *  @param number_of_states                 The number of states that searched for (including the ground state).
     */
    CI(const ONVBasis& onv_basis, const size_t number_of_states = 1) :
        onv_basis {onv_basis},
        number_of_states {number_of_states} {}


    /*
     *  MARK: Optimization
     */

    /**
     *  Optimize the electronic structure model.
     * 
     *  @tparam Solver              the type of the solver
     * 
     *  @param solver               The solver that will try to optimize the parameters.
     */
    template <typename Solver>
    QCStructure<LinearExpansion<Scalar, ONVBasis>, Scalar> optimize(Solver& solver, EigenproblemEnvironment<Scalar>& environment) const {

        // The CI method's responsibility is to try to optimize the parameters of its method, given a solver and associated environment.
        solver.perform(environment);


        // Extract the requested number of eigenpairs from the environment and place them into the LinearExpansion wave function model.
        const auto eigenpairs = environment.eigenpairs(this->number_of_states);

        std::vector<LinearExpansion<Scalar, ONVBasis>> linear_expansions {};
        linear_expansions.reserve(number_of_states);

        std::vector<Scalar> energies {};
        energies.reserve(this->number_of_states);

        for (const auto& eigenpair : eigenpairs) {
            linear_expansions.emplace_back(onv_basis, eigenpair.eigenvector());
            energies.push_back(eigenpair.eigenvalue());
        }


        // Wrap all the requested number of states into a QCStructure.
        // Since we have already created a list of LinearExpansions, we only have to create a list of the corresponding energies.
        return QCStructure<LinearExpansion<Scalar, ONVBasis>, Scalar>(energies, linear_expansions);
    }
};


}  // namespace QCMethod
}  // namespace GQCP
