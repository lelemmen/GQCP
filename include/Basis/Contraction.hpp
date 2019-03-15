#ifndef Contraction_hpp
#define Contraction_hpp


#include <vector>


namespace GQCP {


/*
 *  name taken from libint, not a 'real' contraction
 */
struct Contraction {
    size_t l;  // angular momentum (x + y + z)
    std::vector<double> coefficients;  // contraction coefficients


    /**
     *  @return the length of the contraction, i.e. the number of contraction coefficients or the number of primitives it corresponds to
     */
    size_t length() const { return this->coefficients.size(); }
};


}  // namespace GQCP


#endif  /* Contraction_hpp */
