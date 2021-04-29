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


#include "Operator/FirstQuantized/DiamagneticOperator.hpp"
#include "Operator/FirstQuantized/FQMolecularHamiltonian.hpp"
#include "Operator/FirstQuantized/OrbitalZeemanOperator.hpp"
#include "Operator/FirstQuantized/SpinZeemanOperator.hpp"


namespace GQCP {


/**
 *  The first-quantized, molecular electronic Pauli Hamiltonian for systems in a magnetic field.
 */
class FQMolecularPauliHamiltonian:
    public FQMolecularMagneticHamiltonian {
private:
    // The spin Zeeman operator.
    SpinZeemanOperator SZ;


public:
    /*
     *  MARK: Constructors
     */

    /**
     *  @param T            The kinetic energy operator.
     *  @param OZ           The orbital Zeeman operator.
     *  @param D            The diamagnetic operator.
     *  @param SZ           The spin Zeeman operator.
     *  @param V            The nuclear attraction operator.
     *  @param g            The two-electron repulsion operator.
     */
    FQMolecularPauliHamiltonian(const KineticOperator& T, const OrbitalZeemanOperator& OZ, const DiamagneticOperator& D, const SpinZeemanOperator& SZ, const NuclearAttractionOperator& V, const CoulombRepulsionOperator& g);


    /**
     *  Construct a `FQMolecularPauliHamiltonian` from a molecule and underlying homogeneous magnetic field.
     * 
     *  @param molecule         The molecule.
     *  @param B                The external, homogeneous magnetic field.
     */
    FQMolecularPauliHamiltonian(const Molecule& molecule, const HomogeneousMagneticField& B);


    /*
     *  MARK: Operators
     */

    /**
     *  @return The spin Zeeman operator.
     */
    const SpinZeemanOperator& spinZeeman() const { return this->SZ; }
};


}  // namespace GQCP
