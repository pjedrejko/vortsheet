#pragma once
#include<vector>
#include <cstddef>
#include <string>

/*
stores nodes (positions or velocities) as SoA
*/

using NodeWeight = double;
using NodeIdx = int;

class Nodes {
public:
    std::vector<double> r; //radial coordinate
    std::vector<double> z; //vecrtical coordinate
    std::vector<double> g; //circulation density (per parameter)

    Nodes() = default;

    void reserve(size_t capacity) {
        r.reserve(capacity);
        z.reserve(capacity);
        g.reserve(capacity);
    }

    void clear() {
        r.clear();
        z.clear();
        g.clear();
    }

    void resize(size_t count) {
        r.resize(count);
        z.resize(count);
        g.resize(count);
    }
    
    void assign(size_t count, double val) {
        r.assign(count, val);
        z.assign(count, val);
        g.assign(count, val);
    }

    std::size_t size() const {
        return r.size();
    }

    NodeIdx add(double r_val, double z_val, double g_val) {
        r.push_back(r_val);
        z.push_back(z_val);
        g.push_back(g_val);
        return NodeIdx(size()) - 1;
    }

    void toFile  (int step, const std::string& prefix) const;
    void fromFile(int step, const std::string& prefix);

};
