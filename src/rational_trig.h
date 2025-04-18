// These constants include approximate discretizations of trig function pairs as rational numbers:
// int nc, ns, den;
//
// psi       = atan2(ns, nc)
// cos(psi) ~= dcos(psi) = nc / den
// sin(psi) ~= dsin(psi) = ns / den
//
// The error is:
// e         = (nc/den)^2 + (ns/den)^2 - 1
//
// For a tolerance of sqrt(|e|) <= 1E-2, the following values are
// sampled from the pareto frontier for values of the "angle gap"
// (i.e. maximum change in angle from one value of psi to the next)
//
// NOTE: only values in range (0, 90) deg are stored
// Other values for nc/ns are found by applying appropriate sign change
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <numbers>

using namespace std::numbers;

struct RationalTrig{
    unsigned short int den;  // Denominator of rational fraction

    // Numerators of rational fractions corresponding to
    // sin(theta) = ns/den, cos(theta) = nc/den
    //
    // ``fac'' stores greatest common divisor between nc and ns, resulting in integer after division
    std::vector<unsigned short int> ns, nc, fac;
    std::vector<double> tha;

    // Constructors
    RationalTrig() : den(1), ns({0, 1}), nc({1, 0}), fac({1, 1}), tha({0., 0.5*pi}) {};
    RationalTrig(unsigned short int den);
    RationalTrig(unsigned short int den, std::vector<unsigned short int> ns, std::vector<unsigned short int> nc, std::vector<unsigned short int> fac) : den(den), ns(ns), nc(nc), fac(fac) {
        // Calculate angle corresponding to each trig value
        tha.reserve(ns.size());
        for (unsigned short int idx = 0; idx < ns.size(); idx++) {
            tha.emplace_back(atan2(ns[idx], nc[idx]));
        }  // for
    };
    RationalTrig(unsigned short int den, std::vector<unsigned short int> ns, std::vector<unsigned short int> nc, std::vector<unsigned short int> fac, std::vector<double> tha) : den(den), ns(ns), nc(nc), fac(fac), tha(tha) {};

    std::string to_string() const {
        std::stringstream ss;
        ss << "RationalTrig[N/" << den << "]"; 
        return ss.str();
    };
};  // struct (RationalTrig)
