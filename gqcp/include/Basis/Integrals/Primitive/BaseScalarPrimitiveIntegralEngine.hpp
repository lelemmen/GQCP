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


#include <cstddef>


namespace GQCP {


/**
 *  A base class for primitive integral engines that can calculate integrals over scalar operators.
 */
class BaseScalarPrimitiveIntegralEngine {
public:
    /*
     *  MARK: Destructor
     */

    /**
     *  A pure virtual destructor in order to make this class abstract.
     */
    virtual ~BaseScalarPrimitiveIntegralEngine() = 0;


    /*
     *  MARK: Components
     */

    // The number of components the operator has. For a scalar operator, this is equal to 1.
    static constexpr size_t Components = 1;

    /**
     *  Prepare this engine's internal state such that it is able to calculate integrals over the given component of the operator.
     * 
     *  @param component                The index of the component of the operator.
     * 
     *  @note Since a scalar operator has only 1 component, this method has no effect.
     */
    void prepareStateForComponent(const size_t component);
};


}  // namespace GQCP
