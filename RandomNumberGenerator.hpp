//
// Created by Xingyu Zhu on 2022/9/3.
//

#ifndef HW1_PRNG_HPP
#define HW1_PRNG_HPP

#include <cmath>
#include <numeric>
#include <utility>


class RandomNumberGenerator{
public:
    virtual double operator()() = 0;
};

class LinearCongruential: public RandomNumberGenerator {
protected:
    unsigned long long seed;
    unsigned long long a;
    unsigned long long c;
    unsigned long long k;

public:
    explicit LinearCongruential(unsigned long long seed_=1, unsigned long long a_=39373, long long c_=0, unsigned long k_= 2147483647):
            seed(seed_),a(a_),c(c_),k(k_){}

    double operator()() override {
        seed = (a * seed + c) % k;
        return (double)seed / (double)k;
    }
};

class InverseTransform :public RandomNumberGenerator {
protected:
    double a0 = 2.50662823884, a1 = -18.61500062529, a2 = 41.39119773534, a3 = -25.44106049637;
    double b0 = -8.47351093090, b1 = 23.08336743743, b2 = -21.06224101826, b3 = 3.13082909833;
    double c0 = 0.3374754822726147, c1 = 0.9761690190917186, c2 = 0.1607979714918209,
            c3 = 0.0276438810333863, c4 = 0.0038405729373609, c5 = 0.0003951896511919,
            c6 = 0.0000321767881768, c7 = 0.0000002888167364, c8 = 0.0000003960315187;
    LinearCongruential uniform_distribution;

public:

    explicit InverseTransform(LinearCongruential U_ = LinearCongruential()) :uniform_distribution(std::move(U_)) {}

    double inverse_norm(double u) const {
        double y = u - 0.5;
        double r, x;
        if (abs(y) < 0.42) {
            r = y * y;
            x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1);
        } else {
            r = u;
            if (y > 0)
                r = 1 - u;
            r = log(-log(r));
            x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
            if (y < 0)
                x = -x;
        }
        return x;
    }

    double operator()() override {
        double u = uniform_distribution();
        return inverse_norm(u);
    }
};

class AcceptanceRejection :public RandomNumberGenerator {
protected:
//    double c;
    LinearCongruential uniform_distribution;

public:
    explicit AcceptanceRejection(LinearCongruential U_ = LinearCongruential()):uniform_distribution(std::move(U_)) {}
//    explicit AcceptanceRejection(LinearCongruential U_ = LinearCongruential()):c(sqrt(2 *  M_E/ M_PI)), uniform_distribution(std::move(U_)) {}

    double operator()() override{
        double u1 = uniform_distribution();
        double u2 = uniform_distribution();
        double u3 = uniform_distribution();
        double x = -log(u1);
        while (u2 > exp(-0.5 * (x - 1)*(x - 1))) {
            u1 = uniform_distribution();
            u2 = uniform_distribution();
            u3 = uniform_distribution();
            x = -log(u1);
        }
        if (u3 <= 0.5)
            x = -x;
        return x;
    }
};

class BoxMuller :public RandomNumberGenerator {
protected:
    LinearCongruential uniform_distribution;
    bool flag;
    double prev;

public:
    explicit BoxMuller(LinearCongruential U_ = LinearCongruential()) :uniform_distribution(std::move(U_)), flag(false), prev(0.) {}

    double operator()() override {
        if (flag){
            flag = false;
            return prev;
        }
        double X = 2.;
        double u1, u2;
        while (X > 1) {
            u1 = uniform_distribution();
            u2 = uniform_distribution();
            u1 = 2 * u1 - 1;
            u2 = 2 * u2 - 1;
            X = u1 * u1 + u2 * u2;
        }
        double Y = sqrt(-2 * log(X) / X);
        prev = u2 * Y;
        return u1 * Y;
    }
};


#endif //HW1_PRNG_HPP
