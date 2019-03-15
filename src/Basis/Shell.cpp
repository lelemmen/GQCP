#include "Basis/Shell.hpp"


namespace GQCP {


/*
 *  CONSTRUCTORS
 */

Shell::Shell(const Atom& atom, const std::vector<double>& exponents, const std::vector<Contraction>& contractions) :
    atom (atom),
    exponents (exponents),
    contractions (contractions)
{
    for (const auto& contraction : this->contractions) {
        if (contraction.coefficients.size() != this->exponents.size()) {
            throw std::invalid_argument("Shell(Atom, std::vector<double>, std::vector<Contraction>): the exponents and contractions must match in size.");
        }
    }
}



/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of contractions corresponding to this shell
 */
size_t Shell::numberOfContractions() const {
    return this->contractions.size();
}


/**
 *  @return the number of basis functions that are in this shell
 */
size_t Shell::numberOfBasisFunctions() const {

    size_t nbf {};
    for (const auto& contraction : this->contractions) {
        nbf += contraction.numberOfBasisFunctions();
    }

    return nbf;
}


}  // namespace GQCP
