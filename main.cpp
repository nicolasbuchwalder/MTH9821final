//
//  main.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

/*
#include <iostream>
#include <iomanip>
#include "Option.hpp"
#include "FiniteDifference.hpp"

void PrintVector(const std::vector<double>& vec) {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}




int main(int argc, const char * argv[]) {

    DivsTuple divs;
    divs.push_back(std::make_tuple(true, 2./12., .01));

    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::vanilla, 48., 50., .5, .25, .03, .01, DivsTuple(), std::vector<double>());


    FiniteDifference fd;
    fd.set_params(o, 1000, 0.4);
    auto res = fd.price_option(Scheme::cn_lu, true);
    //fd.show_domain_params();
    //fd.show_grid(false);
    fd.print(res);
    return 0;
}
*/

#include <iostream>
#include <iomanip>
#include "Option.hpp"
#include "MonteCarlo.hpp"

int main(int argc, const char *argv[])
{

    DivsTuple divs;
    std::vector<double> add_params;
    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::vanilla, 48., 50., .5, .25, .03, .0, divs, add_params);

    MonteCarlo mc;
    std::cout << "Done!" << std::endl;
    return 0;
}
