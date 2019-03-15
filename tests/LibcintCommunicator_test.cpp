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
#define BOOST_TEST_MODULE "LibcintCommunicator"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>  // include this to get main(), otherwise the compiler will complain

#include "LibcintCommunicator.hpp"


#include "Molecule.hpp"
#include "Basis/AOBasis.hpp"
#include "LibintCommunicator.hpp"


BOOST_AUTO_TEST_CASE ( sandbox ) {

    GQCP::Atom h1 (1,  0.0, 0.0, 0.8);  // coordinates in bohr
    GQCP::Atom h2 (1,  0.0, 0.0, -0.8);
    GQCP::Molecule mol ({h1, h2});

    GQCP::AOBasis aobasis (mol, "STO-3G");
    auto S = GQCP::LibintCommunicator::get().calculateOverlapIntegrals(aobasis);

    GQCP::LibcintCommunicator libcint;
    BOOST_CHECK(S.isApprox(libcint.calculateOverlapIntegrals(), 1.0e-08));
}
