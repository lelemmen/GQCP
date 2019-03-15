#ifndef BasisSet_hpp
#define BasisSet_hpp


#include <vector>

#include "Basis/Shell.hpp"
#include "Molecule.hpp"


namespace GQCP {


class BasisSet : public std::vector<Shell> {
public:
    using std::vector<Shell>::vector;  // inherit base constructors


public:
    // CONSTRUCTORS
    BasisSet(const std::string& basisset_name, const Molecule& molecule);


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
     *  @param atom     the atom to be found
     *
     *  @return the index 
     */
    size_t findAtom(const Atom& atom) const;
};


}  // namespace GQCP


#endif  /* BasisSet_hpp */
