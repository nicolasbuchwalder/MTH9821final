//
//  main.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//


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
    std::cout << std::fixed << std::setprecision(6);
    
    DivsTuple divs;
    divs.push_back(std::make_tuple(true, 5./12., .01));
    
    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::downout, 52., 50., 1., .2, .03, 0. , divs, std::vector<double>());
    
    
    std::vector<std::size_t> Ms{4};//, 16, 64, 256};
    for (auto M : Ms){
        FiniteDifference fd;
        fd.set_params(o, M, 0.4);
        auto res = fd.price_option(Scheme::eul_expl);
        //fd.show_domain_params();
        //fd.show_grid(false);
        fd.print(res);
    }
    
    return 0;
}


//#include <iostream>
//#include <iomanip>
//#include "Option.hpp"
//#include "MonteCarlo.hpp"
//
//int main(int argc, const char *argv[])
//{
//
//    DivsTuple divs;
//    std::vector<double> add_params;
//    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::vanilla, 48., 50., .5, .25, .03, .0, divs, add_params);
//
//    MonteCarlo mc;
//    std::cout << "Done!" << std::endl;
//    return 0;
//}
