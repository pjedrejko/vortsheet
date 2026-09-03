#include "Sheet.hpp"
#include<stdexcept>
#include<string>
#include"Utilities.hpp"

void Sheet::Integrals::compute(){
    time           += dt;
    step           += 1;
    parameterLength = computeParameterLength();
    physicalLength  = computePhysicalLength();
    circulation     = computeCirculation();
    impulse         = computeImpulse();
    volume          = computeVolume();
    axisThickness   = computeAxisThickness();
    dt              = computeTimeStep(circulation, axisThickness);
    //time           += dt;
}

double Sheet::Integrals::computeAxisThickness(){
    return sheet.x.z[sheet.G.topAxis] - sheet.x.z[sheet.G.bottomAxis];
}

double Sheet::Integrals::computeParameterLength(){
    double len = 0.0;
    const int n = sheet.x.size();

    for(NodeIdx i = 0; i < n; i++)
        len += sheet.G.W[i];

    return len;
}


double Sheet::Integrals::computePhysicalLength(){
    double len = 0.0;

    for(const Segment& s: sheet.G.seg){
        double dr = sheet.x.r[s.end] - sheet.x.r[s.start];
        double dz = sheet.x.z[s.end] - sheet.x.z[s.start];
        len += std::sqrt(sq(dr) + sq(dz)) * s.layers;
    }

    return len;
}

double Sheet::Integrals::computeCirculation(){
    double circ = 0.0;
    const int n = sheet.x.size();

    for(NodeIdx i = 0; i < n; i++)
        circ += sheet.x.g[i] * sheet.G.W[i];

    return circ;
}

double Sheet::Integrals::computeImpulse(){
    double imp = 0.0;
    const int n = sheet.x.size();

    for(NodeIdx i = 0; i < n; i++)
        imp += sheet.x.g[i] * sheet.G.W[i] * sq(sheet.x.r[i]);

    imp *= (2 * M_PI);
    return imp;
}

double Sheet::Integrals::computeVolume(){
    double vol = 0.0;
    for(const Segment& s: sheet.G.seg){
        if( !s.interfacial ) continue;
        double r  = (sheet.x.r[s.end] + sheet.x.r[s.start]) * 0.5;
        double dz = (sheet.x.z[s.end] - sheet.x.z[s.start]);
        vol += sq(r) * dz;
    }
    vol *= M_PI;
    return vol;
}

double Sheet::Integrals::computeTimeStep(double circulation, double axisThickness){
    double dt1 = Params::kt1 * std::abs(circulation) / axisThickness;
    double dt2 = Params::kt2 * Params::delta / std::abs(circulation);

    return std::min(dt1, dt2);
    //return 0.02;
}

void Sheet::Integrals::toFile(){

    const char* filepath = "results/integrals";
    logg.formPrint("writing to results/integrals... ");

    static bool firstCall = true;

    FILE* f = std::fopen(filepath, "a");
    if (!f)
        throw std::runtime_error("Could not open file: results/integrals");

    int colWidth = 25;
    if (firstCall) {
        std::fprintf(f, "%% %-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s %-*s\n",
            colWidth, "step", 
            colWidth, "time", 
            colWidth, "dt", 
            colWidth, "circulation", 
            colWidth, "axisThickness", 
            colWidth, "impulse", 
            colWidth, "volume",
            colWidth, "parameterLength",
            colWidth, "physicalLength");
       firstCall = false;
    }

    std::fprintf(f, "  %-*d %-*.17e %-*.17e %-*.17e %-*.17e %-*.17e %-*.17e %-*.17e %-*.17e\n",
            colWidth, step, 
            colWidth, time, 
            colWidth, dt, 
            colWidth, circulation, 
            colWidth, axisThickness, 
            colWidth, impulse, 
            colWidth, volume,
            colWidth, parameterLength,
            colWidth, physicalLength);
    std::fclose(f);    
    logg.print("done.\n");
}

void Sheet::Integrals::info() {
    logg.print("| step:          %8d |\n",     step);
    logg.print("| time:          %8.5f |\n",   time);
    logg.print("| nodes:         %8d |\n",     int(sheet.x.r.size()));
    logg.print("| segments:      %8d |\n",     int(sheet.G.seg.size()));
    logg.print("| dt:            %8.5f |\n",   dt);
    logg.print("| circulation:   %8.5f |\n",   circulation);
    logg.print("| axisThickness: %8.5f |\n",   axisThickness);
    logg.print("| impulse:       %8.5f |\n",   impulse);
    logg.print("| volume:        %8.5f |\n\n", volume);
}

void Sheet::Integrals::fromFile(int stp) {
    const char* filepath = "results/integrals";

    logg.formPrint(1, "reading step %d from results/integrals ", stp);

    FILE* f = std::fopen(filepath, "r");
    if (!f) 
        throw std::runtime_error("Could not open file: results/integrals");

    char line[1024];
    bool found = false;

    while (std::fgets(line, sizeof(line), f)) {
        //skip header lines or comments starting with '%'
        if (line[0] == '%' || line[0] == '#' || line[0] == '\n') {
            continue;
        }

        int readStep;
        double readTime, readDt, readCirc, readAxisThick, readImp, readVol, readParamLen, readPhysLen;

        int matched = std::sscanf(line, "%d %lf %lf %lf %lf %lf %lf %lf %lf",
                                  &readStep, 
                                  &readTime, 
                                  &readDt, 
                                  &readCirc, 
                                  &readAxisThick, 
                                  &readImp, 
                                  &readVol, 
                                  &readParamLen, 
                                  &readPhysLen);

        if (matched == 9 && readStep == stp) {
            step            = readStep;
            time            = readTime;
            dt              = readDt;
            circulation     = readCirc;
            axisThickness   = readAxisThick;
            impulse         = readImp;
            volume          = readVol;
            parameterLength = readParamLen;
            physicalLength  = readPhysLen;

            found = true;
            break; //stop after finding the targeted step
        }
    }

    std::fclose(f);

    if(!found)
        throw std::runtime_error("Could not find step " + std::to_string(stp) + " in results/integrals");

    logg.print("done.\n");
}
