//
//  main.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//

#include <iostream>
#include <iomanip>
#include "Option.hpp"
#include "FiniteDifferencePricer.hpp"

void PrintVector(const std::vector<double>& vec) {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}


int main(int argc, const char * argv[]) {
    
    DivsTuple divs;
    divs.push_back(std::make_tuple(true, 2./12., .01));
    
    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::vanilla, 48., 50., .5, .25, .03, .0, divs, std::vector<double>());
    
    
    return 0;
}
