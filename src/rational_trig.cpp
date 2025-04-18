#include "rational_trig.h"
#include <iostream>

unsigned short int gcd(unsigned short int a, unsigned short int b) {
    int r;
    while (b != 0) {
        r = a % b;
        a = b;
        b = r;
    }
    return a;  // TODO
}  // gcd

RationalTrig::RationalTrig(unsigned short int den) : den(den) {
    std::vector<unsigned short int> n1, n2;
    n1.reserve(den);
    n2.reserve(den);

    unsigned short int den2 = den*den;
    for (unsigned short int idx = 0; idx <= den; idx++) {
        n1.emplace_back(idx);
        n2.emplace_back(static_cast<unsigned short int>(round(sqrt(den2 - idx*idx))));
    }  // for (populate n1, n2)

    // Create map to sort unique angle values in increasing order from 0 to 90 deg
    std::map<double, std::pair<unsigned short int, unsigned short int>> trig_map;

    std::vector<unsigned short int>::iterator n1_it = n1.begin();
    for (std::vector<unsigned short int>::iterator n2_it = n2.begin(); n2_it != n2.end(); n1_it++, n2_it++) {
        trig_map.insert({atan2(*n1_it, *n2_it), {*n1_it, *n2_it}});
        trig_map.insert({atan2(*n2_it, *n1_it), {*n2_it, *n1_it}});
    }  // for (create unique sorted map)

    // Initialize trig data from unique map values
    ns.reserve(trig_map.size());
    nc.reserve(trig_map.size());
    fac.reserve(trig_map.size());
    tha.reserve(trig_map.size());

    for (std::map<double, std::pair<unsigned short int, unsigned short int>>::iterator map_it = trig_map.begin(); map_it != trig_map.end(); map_it++) {
        tha.emplace_back(map_it->first);
        ns.emplace_back(map_it->second.first);
        nc.emplace_back(map_it->second.second);
        fac.emplace_back(gcd(map_it->second.first, map_it->second.second));
    }  // for (initialize trig values)
}  // constructor (RationalTrig)
