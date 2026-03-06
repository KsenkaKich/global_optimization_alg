#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "discpp.h"
#include "structs.h"

using namespace std;

class OptimizationPlotter {
private:
    Dislin g;

public:
    void PlotAlgorithm(const std::vector<Trial>& trials, const Trial& best_point,
                      double a, double b, const std::string& algorithm_name) {
        
        if (trials.empty()) return;
        
        std::vector<Trial> sorted_trials = trials;
        std::sort(sorted_trials.begin(), sorted_trials.end());
        
        std::vector<double> x_all(sorted_trials.size());
        std::vector<double> y_all(sorted_trials.size());
        
        for (size_t i = 0; i < sorted_trials.size(); i++) {
            x_all[i] = sorted_trials[i].x;
            y_all[i] = sorted_trials[i].z;
        }
        
        double min_y = *std::min_element(y_all.begin(), y_all.end());
        double max_y = *std::max_element(y_all.begin(), y_all.end());
        double y_padding = (max_y - min_y) * 0.2;
        if (y_padding < 0.001) y_padding = 0.5;
        
        g.metafl("cons");
        g.scrmod("revers");
        g.disini();
        g.pagera();
        g.complx();
        
        g.axspos(450, 1800);
        g.axslen(2200, 1200);
        
        g.name("x", "x");
        g.name("f(x)", "y");
        
        g.labdig(-1, "x");
        g.ticks(5, "x");
        g.ticks(5, "y");
        
        g.titlin(const_cast<char*>(algorithm_name.c_str()), 1);
        
        std::string iter_str = "Iterations: " + std::to_string(trials.size());
        g.titlin(const_cast<char*>(iter_str.c_str()), 2);
        
        std::string best_info = "Best: x = " + std::to_string(best_point.x).substr(0, 8) +
                               ", f(x) = " + std::to_string(best_point.z).substr(0, 8) +
                               ", k = " + std::to_string(best_point.k);
        g.titlin(const_cast<char*>(best_info.c_str()), 3);
        
        int ic = g.intrgb(0.95, 0.95, 0.95);
        g.axsbgd(ic);
        
        g.graf(a, b, a, (b-a)/10, 
               min_y - y_padding, max_y + y_padding, 
               min_y - y_padding, (max_y - min_y + 2*y_padding)/10);
        
        g.setrgb(0.8, 0.8, 0.8);
        g.grid(3, 3);
        
        g.color("fore");
        g.height(50);
        g.title();
        
        g.setrgb(0.5, 0.5, 0.5);
        g.curve(x_all.data(), y_all.data(), sorted_trials.size());
        
        g.color("blue");
        g.incmrk(1);
        g.marker(1);
        g.hsymbl(15);
        g.curve(x_all.data(), y_all.data(), sorted_trials.size());
        
        g.disfin();
    }
};