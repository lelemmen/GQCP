#ifndef Shell_hpp
#define Shell_hpp


#include "Atom.hpp"
#include "Basis/Contraction.hpp"


namespace GQCP {


class Shell {
private:
    Atom atom;
    std::vector<double> exponents;  // shared for every contraction
    std::vector<Contraction> contractions;


public:
    // CONSTRUCTORS
    Shell(const Atom& atom, const std::vector<double>& exponents, const std::vector<Contraction>& contractions);


    // GETTERS
    const Atom& get_atom() const { return this->atom; }
    const std::vector<double>& get_exponents() const { return this->exponents; }
    const std::vector<Contraction>& get_contractions() const { return this->contractions; }


    // PUBLIC METHODS
    /**
     *  @return the number of contractions corresponding to this shell
     */
    size_t numberOfContractions() const;

    /**
     *  @return the number of basis functions that are in this shell
     */
    size_t numberOfBasisFunctions() const;
};


}  // namespace GQCP


#endif  /* Shell_hpp */
