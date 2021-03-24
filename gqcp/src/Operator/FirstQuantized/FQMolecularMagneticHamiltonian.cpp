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

#include "Operator/FirstQuantized/FQMolecularMagneticHamiltonian.hpp"


namespace GQCP {


/*
 *  MARK: Constructors
 */

/**
 *  @param T            The kinetic energy operator.
 *  @param P            The paramagnetic operator.
 *  @param D            The diamagnetic operator.
 *  @param V            The nuclear attraction operator.
 *  @param g            The two-electron repulsion operator.
 */
FQMolecularMagneticHamiltonian::FQMolecularMagneticHamiltonian(const KineticOperator& T, const ParamagneticOperator& P, const DiamagneticOperator& D, const NuclearAttractionOperator& V, const CoulombRepulsionOperator& g) :
    P {P},
    D {D},
    FQMolecularHamiltonian(T, V, g) {}


/**
 *  Construct a `FQMolecularMagneticHamiltonian` from a molecule and underlying homogeneous magnetic field.
 * 
 *  @param molecule         The molecule.
 *  @param B                The external, homogeneous magnetic field.
 */
FQMolecularMagneticHamiltonian::FQMolecularMagneticHamiltonian(const Molecule& molecule, const HomogeneousMagneticField& B) :
    FQMolecularMagneticHamiltonian(KineticOperator(), ParamagneticOperator(B), DiamagneticOperator(B), NuclearAttractionOperator(molecule.nuclearFramework()), CoulombRepulsionOperator()) {}


}  // namespace GQCP
