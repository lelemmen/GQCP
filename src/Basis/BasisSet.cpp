#include "Basis/BasisSet.hpp"


#include "LibintCommunicator.hpp"


#include <algorithm>



namespace GQCP {


BasisSet::BasisSet(const std::string& basisset_name, const Molecule& molecule) {
    const auto &atoms = molecule.get_atoms();

    // TO:DO no longer use libint2 to read this
    libint2::BasisSet libint_basis((basisset_name), LibintCommunicator::get().interface(molecule.get_atoms()));
    this->number_of_basis_functions = (static_cast<size_t>(libint_basis.nbf()));

    this->reserve(libint_basis.size());

    // copy libint2 shells to gqcp shells
    for (const auto &libint_shell : libint_basis) {

        std::vector<Contraction> contractions;
        contractions.reserve(libint_shell.contr.size());
        std::cout << libint_shell.contr.size();

        // copy libint2 contractions (contr) to gqcp contractions
        for (const auto &contraction : libint_shell.contr) {
            contractions.push_back({static_cast<size_t>(contraction.l), contraction.coeff});
        }

        // libint2 only stores the carthesian origin of the shell (not the atom)
        //  find the atom corresponding to the copied shell's origin
        Atom corresponding_atom;
        for (const Atom &atom : atoms) {
            Eigen::Map<const Eigen::Matrix<double, 3, 1>> libint_origin_map(libint_shell.O.data());
            if (atom.position.isApprox(libint_origin_map)) {
                corresponding_atom = atom;
                break;
            }
        }

        this->emplace_back(corresponding_atom, libint_shell.alpha, contractions);
    }
}


/*
 *  PUBLIC METHODS
 */

/**
 *  @return the number of shells in this basisset
 */
size_t BasisSet::numberOfShells() const {
    return this->size();
}


/**
 *  @return the number of basis functions in this basisset
 */
size_t BasisSet::numberOfBasisFunctions() const {
    size_t number_of_basis_functions = 0;
    for (const auto& shell : *this) {
        for (const auto& contraction : shell.get_contractions()) {
            number_of_basis_functions += contraction.numberOfBasisFunctions();
        }
    }
    return number_of_basis_functions;
}


/**
 *  @return an ordered vector of the unique atoms in this basisset
 */
std::vector<Atom> BasisSet::atoms() const {

    std::vector<Atom> atoms {};

    // Append every unique atom in this basisset's shells
    for (const auto& shell : *this) {
        auto atom = shell.get_atom();

        auto p = std::find(atoms.begin(), atoms.end(), atom);
        if (p == atoms.end()) {  // if unique
            atoms.push_back(atom);
        }
    }

    return atoms;
}


/**
 *  @param shell_index      the index of the shell
 *
 *  @return the (total basis function) index that corresponds to the first basis function in the given shell
 */
size_t BasisSet::basisFunctionIndex(size_t shell_index) const {

    size_t bf_index {};

    // Count the number of basis functions before 
    for (size_t i = 0; i < shell_index; i++) {
        bf_index += this->operator[](i).numberOfBasisFunctions();
    }

    return bf_index;
}


}  // namespace GQCP
