//
//  main.cpp
//  MTH9815_HW12_FD
//
//  Created by Nicolas Buchwalder on 28.11.22.
//


#include <iostream>
#include <iomanip>
#include "Option.hpp"
#include "BinomialTree.hpp"
#include "TrinomialTree.hpp"
#include "MonteCarlo.hpp"
#include "FiniteDifference.hpp"


void PrintVector(const std::vector<double>& vec) {
    for (auto elem : vec) {
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
}


/**
 * Returns the exact price of a down-and-out Call option using Black-Scholes prices.
 * Change this function as necessary based on the option you are pricing.
 */
std::vector<double> price_option_closedform(OptionExercise ex, OptionPayoff payoff, OptionType type, double S, double K, double T, double sigma, double r, double q, DivsTuple divs, std::vector<double> add_params)
{
    Option o1(ex, payoff, type, S, K, T, sigma, r, q, divs, add_params);
    std::vector<double> price_BS1 = o1.price_european();

    Option o2(ex, payoff, type, pow(add_params[0], 2) / S, K, T, sigma, r, q, divs, add_params);
    std::vector<double> price_BS2 = o2.price_european();

    double a = (r - q) / pow(sigma, 2) - 0.5;
    double price = price_BS1[0] - pow(add_params[0] / S, 2 * a) * price_BS2[0];
    std::vector<double> prc{price};
    return prc;
}

int main(int argc, const char *argv[])
{
    std::cout << std::fixed << std::setprecision(6);
    
    Option o(OptionExercise::euro, OptionPayoff::put, OptionType::vanilla, 41., 39., 1., 0.25, 0.03, 0.005, DivsTuple(), std::vector<double>());
    TrinomialTree ttree;
    std::vector<int> Ns = {10, 20, 40, 80, 160, 320, 640, 1280};
    //std::vector<Method> methods = {average} //{standard, average, bbs, bbsr};
    for (auto N : Ns){
    //for (int N = 10; N <= 100; N++){
        ttree.set_params(o, N, 0.);
        std::cout << ttree.price_option(standard)[0] << std::endl;
    }

    
    
    
    
    
    
//    OptionExercise ex = OptionExercise::euro;
//    OptionPayoff payoff = OptionPayoff::call;
//    OptionType type = OptionType::vanilla;
//    double S = 39.0;
//    double K = 39.0;
//    double T = 0.75;
//    double sigma = 0.25;
//    double r = 0.02;
//    double q = 0.01;
//    DivsTuple divs;
//    std::vector<double> add_params{36.0};
//    RandomNumberGeneratorEnum rng = RandomNumberGeneratorEnum::BoxMullerRNG;
//
//    Option o(ex, payoff, type, S, K, T, sigma, r, q, divs, add_params);
//    std::vector<double> price_BS = price_option_closedform(ex, payoff, type, S, K, T, sigma, r, q, divs, add_params);
//
//    MonteCarlo mc{o, rng};
//
//    std::size_t len = 10000 * pow(2, 9);
//    std::size_t m1 = 200;
//    for (std::size_t k = 0; k <= 9; k++)
//    {
//        std::size_t n1 = 50 * pow(2, k);
//        std::vector<double> price_MC1 = mc.price_option(n1, m1, false);
//
//        std::size_t N1 = m1 * n1;
//        std::size_t m2 = ceil(pow(N1, 1.0 / 3) * pow(0.75, 2.0 / 3));
//        std::size_t n2 = floor(N1 / m2);
//        std::vector<double> price_MC2 = mc.price_option(n2, m2, false);
//
//        std::cout << N1 << "\t" << m1 << "\t" << n1 << "\t" << price_MC1[0] << "\t" << abs(price_BS[0] - price_MC1[0]) << "\t"
//                  << m2 << "\t" << n2 << "\t" << price_MC2[0] << "\t" << abs(price_BS[0] - price_MC2[0]) << std::endl;
//    }
    
    
    
    
    
    
//    DivsTuple divs;
//    divs.push_back(std::make_tuple(true, 5./12., .01));
//
//    Option o(OptionExercise::euro, OptionPayoff::call, OptionType::downout, 52., 50., 1., .2, .03, 0. , divs, std::vector<double>());
//
//
//    std::vector<std::size_t> Ms{4};//, 16, 64, 256};
//    for (auto M : Ms){
//        FiniteDifference fd;
//        fd.set_params(o, M, 0.4);
//        auto res = fd.price_option(Scheme::eul_expl);
//        //fd.show_domain_params();
//        //fd.show_grid(false);
//        fd.print(res);
//    }
    
    

    return 0;
}
