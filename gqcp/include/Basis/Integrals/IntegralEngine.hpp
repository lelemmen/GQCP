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


#include "Basis/Integrals/Interfaces/LibcintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibcintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintOneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Interfaces/LibintTwoElectronIntegralEngine.hpp"
#include "Basis/Integrals/OneElectronIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveAngularMomentumIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveCanonicalKineticEnergyIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveCoulombRepulsionIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveElectronicDipoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveElectronicQuadrupoleIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveLinearMomentumIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveNuclearAttractionIntegralEngine.hpp"
#include "Basis/Integrals/Primitive/PrimitiveOverlapIntegralEngine.hpp"
#include "Basis/Integrals/TwoElectronIntegralEngine.hpp"
#include "Utilities/aliases.hpp"


namespace GQCP {


/**
 *  A class that creates integral engines, like a factory class.
 */
class IntegralEngine {
public:
    /*
     *  MARK: GQCP ("In-house") integral engines
     * 
     *  These integral engines can only calculate integrals over Cartesian d-shells, not over spherical d-shells.
     */

    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the overlap operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The overlap operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the overlap operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const OverlapOperator& op) -> OneElectronIntegralEngine<PrimitiveOverlapIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveOverlapIntegralEngine<Shell>>(PrimitiveOverlapIntegralEngine<Shell>());
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the electronic dipole operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The electronic dipole operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the electronic dipole operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const ElectronicDipoleOperator& op) -> OneElectronIntegralEngine<PrimitiveElectronicDipoleIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveElectronicDipoleIntegralEngine<Shell>>(PrimitiveElectronicDipoleIntegralEngine<Shell>(op));
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the linear momentum operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The linear momentum operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the linear momentum operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const LinearMomentumOperator& op) -> OneElectronIntegralEngine<PrimitiveLinearMomentumIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveLinearMomentumIntegralEngine<Shell>>(PrimitiveLinearMomentumIntegralEngine<Shell>());
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the canonical kinetic energy operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The canonical kinetic energy operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the canonical kinetic energy operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const KineticOperator& op) -> OneElectronIntegralEngine<PrimitiveCanonicalKineticEnergyIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveCanonicalKineticEnergyIntegralEngine<Shell>>(PrimitiveCanonicalKineticEnergyIntegralEngine<Shell>());
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the angular momentum operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The angular momentum operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the angular momentum operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const AngularMomentumOperator& op) -> OneElectronIntegralEngine<PrimitiveAngularMomentumIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveAngularMomentumIntegralEngine<Shell>>(PrimitiveAngularMomentumIntegralEngine<Shell>(op));
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the electronic quadrupole operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The electronic quadrupole operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the electronic quadrupole operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const ElectronicQuadrupoleOperator& op) -> OneElectronIntegralEngine<PrimitiveElectronicQuadrupoleIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveElectronicQuadrupoleIntegralEngine<Shell>>(PrimitiveElectronicQuadrupoleIntegralEngine<Shell>(op));
    }


    /**
     *  Create an in-house one-electron integral engine that can calculate integrals over the nuclear attraction operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The nuclear attraction operator.
     * 
     *  @return A one-electron integral engine that can calculate integrals over the nuclear attraction operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const NuclearAttractionOperator& op) -> OneElectronIntegralEngine<PrimitiveNuclearAttractionIntegralEngine<Shell>> {

        return OneElectronIntegralEngine<PrimitiveNuclearAttractionIntegralEngine<Shell>>(PrimitiveNuclearAttractionIntegralEngine<Shell>(op));
    };


    /**
     *  Create an in-house two-electron integral engine that can calculate integrals over the Coulomb repulsion operator.
     * 
     *  @tparam Shell           The type of shell that the integrals should be calculated over.
     * 
     *  @param op               The Coulomb repulsion operator.
     * 
     *  @return A two-electron integral engine that can calculate integrals over the Coulomb repulsion operator.
     * 
     *  @note This integral engine can only calculate integrals over Cartesian d-shells, not spherical d-shells.
     */
    template <typename Shell>
    static auto InHouse(const CoulombRepulsionOperator& op) -> TwoElectronIntegralEngine<PrimitiveCoulombRepulsionIntegralEngine<Shell>> {

        return TwoElectronIntegralEngine<PrimitiveCoulombRepulsionIntegralEngine<Shell>>(PrimitiveCoulombRepulsionIntegralEngine<Shell>());
    }


    /*
     *  LIBINT
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libint integral library backend
     */
    static auto Libint(const CoulombRepulsionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintTwoElectronIntegralEngine<CoulombRepulsionOperator::NumberOfComponents>;

    /**
     *  @param op               the electronic electric dipole operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the electronic electric dipole operator using the Libint integral library backend
     */
    static auto Libint(const ElectronicDipoleOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<ElectronicDipoleOperator::NumberOfComponents>;

    /**
     *  @param op               the nuclear attraction operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libint integral library backend
     */
    static auto Libint(const NuclearAttractionOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<NuclearAttractionOperator::NumberOfComponents>;

    /**
     *  @param op               the kinetic operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libint integral library backend
     */
    static auto Libint(const KineticOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<KineticOperator::NumberOfComponents>;

    /**
     *  @param op               the overlap operator
     *  @param max_nprim        the maximum number of primitives per contracted Gaussian shell
     *  @param max_l            the maximum angular momentum of Gaussian shell
     * 
     *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libint integral library backend
     */
    static auto Libint(const OverlapOperator& op, const size_t max_nprim, const size_t max_l) -> LibintOneElectronIntegralEngine<OverlapOperator::NumberOfComponents>;


    /*
     *  LIBCINT
     */

    /**
     *  @param op               the Coulomb repulsion operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a two-electron integral engine that can calculate integrals over the Coulomb repulsion operator using the Libcint integral library backend
     */
    static auto Libcint(const CoulombRepulsionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintTwoElectronIntegralEngine<GTOShell, CoulombRepulsionOperator::NumberOfComponents, double>;

    /**
     *  @param op               the electron electronic dipole operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the electronic dipole operator using the Libcint integral library backend
     */
    static auto Libcint(const ElectronicDipoleOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, ElectronicDipoleOperator::NumberOfComponents, double>;

    /**
     *  @param op               the nuclear attraction operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the nuclear attraction operator using the Libcint integral library backend
     */
    static auto Libcint(const NuclearAttractionOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, NuclearAttractionOperator::NumberOfComponents, double>;

    /**
     *  @param op               the kinetic operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the kinetic operator using the Libcint integral library backend
     */
    static auto Libcint(const KineticOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, KineticOperator::NumberOfComponents, double>;

    /**
     *  @param op               the overlap operator
     *  @param shell_set        the ShellSet whose information should be converted to a RawContainer, which will serve as some kind of 'global' data for the libcint engine to use in all its calculate() calls
     * 
     *  @return a one-electron integral engine that can calculate integrals over the overlap operator using the Libcint integral library backend
     */
    static auto Libcint(const OverlapOperator& op, const ShellSet<GTOShell>& shell_set) -> LibcintOneElectronIntegralEngine<GTOShell, OverlapOperator::NumberOfComponents, double>;
};


}  // namespace GQCP
