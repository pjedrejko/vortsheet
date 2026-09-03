#include "Nodes.hpp"
#include "Params.hpp"
#include<stdexcept>
#include"Utilities.hpp"

void Nodes::toFile(int step, const std::string& prefix) const {
    std::string filename = "results/" + prefix + std::to_string(step);
    
    logg.formPrint(1, "writing to %s {", filename.c_str()).print("\n");

    FILE* f = std::fopen(filename.c_str(), "w");
    if (!f) 
        throw std::runtime_error("Could not open file: " + filename);

    std::fprintf(f, "%% %-21s %-23s %-23s\n", "r", "z", "g");
    for (size_t i = 0; i < r.size(); i++) {
        std::fprintf(f, "%.17e %.17e %.17e\n", r[i], z[i], g[i]);
    }

    std::fclose(f);

    logg.formPrint("%d nodes written", r.size()).print("\n");
    logg.formPrint(-1, "}\n");
}

void Nodes::fromFile(int step, const std::string& prefix) {
    std::string filename = "results/" + prefix + std::to_string(step);
   
    logg.formPrint(1, "reading from %s {", filename.c_str()).print("\n");
    
    FILE* f = std::fopen(filename.c_str(), "r");
    if (!f)
        throw std::runtime_error("Could not open file: " + filename);

    r.clear();
    z.clear();
    g.clear();

    int ch;
    while ((ch = std::fgetc(f)) != EOF) {
        if (ch == '%') {
            int dummy = std::fscanf(f, "%*[^\n]\n");
            (void)dummy;
        } else if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r') {
            std::ungetc(ch, f);
            break;
        }
    }

    double r_val, z_val, g_val;
    while (std::fscanf(f, "%lf %lf %lf", &r_val, &z_val, &g_val) == 3) {
        r.push_back(r_val);
        z.push_back(z_val);
        g.push_back(g_val);
    }

    std::fclose(f);

    logg.formPrint("%d nodes read", r.size()).print("\n");
    logg.formPrint(-1, "}\n");
}