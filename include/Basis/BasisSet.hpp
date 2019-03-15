#ifndef BasisSet_hpp
#define BasisSet_hpp


#include <vector>

#include "Basis/Shell.hpp"
#include "Molecule.hpp"


namespace GQCP {


class BasisSet : public std::vector<Shell> {
private:
    size_t number_of_basis_functions;


public:
    using std::vector<Shell>::vector;  // inherit base constructors


public:
    // CONSTRUCTORS
    BasisSet(const std::string& basisset_name, const Molecule& molecule);


    // GETTERS
    size_t numberOfBasisFunctions() const;


    // PUBLIC METHODS
    /**
     *  @return the number of shells in this basisset
     */
    size_t numberOfShells() const;

    /**
     *  @return an ordered vector of the unique atoms in this basisset
     */
    std::vector<Atom> atoms() const;

    /**
     *  @param shell_index      the index of the shell
     *
     *  @return the (total basis function) index that corresponds to the first basis function in the given shell
     */
    size_t basisFunctionIndex(size_t shell_index) const;
};


}  // namespace GQCP


#endif  /* BasisSet_hpp */
