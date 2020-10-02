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


#include "Operator/SecondQuantized/ModelHamiltonian/HoppingMatrix.hpp"
#include "Operator/SecondQuantized/SQOneElectronOperator.hpp"
#include "Operator/SecondQuantized/SQTwoElectronOperator.hpp"


namespace GQCP {


/**
 *  The Hubbard model Hamiltonian.
 */
template <typename _Scalar>
class HubbardHamiltonian {
public:
    using Scalar = _Scalar;


private:
    HoppingMatrix<Scalar> H;  // the Hubbard hopping matrix


public:
    /*
     *  CONSTRUCTORS
     */

    /**
     *  @param H            the Hubbard hopping matrix
     */
    HubbardHamiltonian(const HoppingMatrix<Scalar>& H) :
        H {H} {}


    /*
     *  PUBLIC METHODS
     */

    /**
     *  @return the core Hamiltonian (i.e. resulting from the Hubbard hopping operator) as a one-electron operator
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQOneElectronOperator<Scalar>> core() const {

        const auto K = this->numberOfLatticeSites();
        SquareMatrix<double> h_par = SquareMatrix<double>::Zero(K);  // 'par' for 'parameters'

        // The one-electron hopping terms can be found on the off-diagonal elements of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            for (size_t q = p; q < K; q++) {
                if (p != q) {
                    h_par(p, q) = this->H(p, q);
                    h_par(q, p) = this->H(q, p);
                }
            }
        }

        return ScalarSQOneElectronOperator<double> {h_par};
    }


    /**
     *  @return the Hubbard hopping matrix for this Hubbard model Hamiltonian
     */
    const HoppingMatrix<Scalar>& hoppingMatrix() const { return this->H; }

    /**
     *  @return the number of lattice sites corresponding used in this Hubbard model Hamiltonian
     */
    size_t numberOfLatticeSites() const { return this->H.numberOfLatticeSites(); }

    /**
     *  @return the two-electron part of the Hamiltonian (resulting from the on-site repulsion) as a two-electron operator
     * 
     *  @note This method is only available for real scalars.
     */
    template <typename Z = Scalar>
    enable_if_t<std::is_same<Z, double>::value, ScalarSQTwoElectronOperator<Scalar>> twoElectron() const {

        const auto K = this->numberOfLatticeSites();

        QCRankFourTensor<double> g_par(K);
        g_par.setZero();

        // The two-electron on-site repulsion is found on the diagonal of the hopping matrix.
        for (size_t p = 0; p < K; p++) {
            g_par(p, p, p, p) = H(p, p);
        }

        return ScalarSQTwoElectronOperator<double> {g_par};
    }
};


}  // namespace GQCP
