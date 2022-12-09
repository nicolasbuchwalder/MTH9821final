#include <iostream>
#include <iomanip>
#include "Option.hpp"
#include "MonteCarlo.hpp"

int main(int argc, const char * argv[]) {
    
    DivsTuple divs;    
    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::vanilla, 48., 50., .5, .25, .03, .0, divs, std::vector<double>());

    MonteCarlo mc;
    std::cout << "Done!" << std::endl;
    return 0;
}
