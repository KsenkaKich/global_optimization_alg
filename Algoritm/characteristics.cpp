#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "structs.h"
#include "discpp.h"

#include "Hill/HillProblem.hpp"
#include "Hill/HillProblemFamily.hpp"
#include "Shekel/ShekelProblem.hpp"
#include "Shekel/ShekelProblemFamily.hpp"

using namespace std;

class CharacteristicsAnalyzer {
private:
    Dislin g;
    double eps;
    int Kmax;
    int num_tests;
    
public:
    CharacteristicsAnalyzer() : eps(1e-4), Kmax(300), num_tests(1000) {}
    
    vector<double> TestHillR(double r) {
        THillProblemFamily hillFam;
        int familySize = min(num_tests, hillFam.GetFamilySize());
        
        vector<int> iterations_needed(familySize);
        
        for (int i = 0; i < familySize; i++) {
            THillProblem hillProb(i);
            
            Task task;
            task.a = 0.0;
            task.b = 1.0;
            task.func = [&hillProb](double x) { 
                return hillProb.ComputeFunction({x}); 
            };
            
            GSASolver solver;
            solver.SetEps(eps);
            solver.SetKmax(Kmax);
            solver.SetR(r);
            solver.SetTask(task);
            
            solver.Solve();
            
            iterations_needed[i] = solver.GetTrials().size();
        }
        
        vector<double> success_rate(Kmax + 1, 0.0);
        for (int iter = 1; iter <= Kmax; iter++) {
            int solved = 0;
            for (int i = 0; i < familySize; i++) {
                if (iterations_needed[i] <= iter) {
                    solved++;
                }
            }
            success_rate[iter] = (double)solved / familySize;
        }
        
        return success_rate;
    }
    
    vector<double> TestShekelR(double r) {
        TShekelProblemFamily shekelFam;
        int familySize = min(num_tests, shekelFam.GetFamilySize());
        
        vector<int> iterations_needed(familySize);
        
        for (int i = 0; i < familySize; i++) {
            TShekelProblem shekelProb(i);
            
            Task task;
            task.a = 0.0;
            task.b = 10.0;
            task.func = [&shekelProb](double x) { 
                return shekelProb.ComputeFunction({x}); 
            };
            
            GSASolver solver;
            solver.SetEps(eps);
            solver.SetKmax(Kmax);
            solver.SetR(r);
            solver.SetTask(task);
            
            solver.Solve();
            
            iterations_needed[i] = solver.GetTrials().size();
        }
        
        vector<double> success_rate(Kmax + 1, 0.0);
        for (int iter = 1; iter <= Kmax; iter++) {
            int solved = 0;
            for (int i = 0; i < familySize; i++) {
                if (iterations_needed[i] <= iter) {
                    solved++;
                }
            }
            success_rate[iter] = (double)solved / familySize;
        }
        
        return success_rate;
    }
    
    void PlotCharacteristics(const vector<double>& success_rate2, const vector<double>& success_rate25,
                            const vector<double>& success_rate3, const vector<double>& success_rate35, const string& title) {
        
        int max_iterations = Kmax;
        
        g.metafl("cons");
        g.scrmod("revers");
        g.disini();
        g.pagera();
        g.complx();
        
        g.axspos(450, 1800);
        g.axslen(2200, 1200);
        
        g.name("number of trials (iterations)", "x");
        g.name("fraction of solved problems", "y");
        
        g.labdig(0, "x");
        g.ticks(10, "x");
        g.ticks(5, "y");
        
        g.titlin(const_cast<char*>(title.c_str()), 1);
        g.titlin("Operational Characteristics", 2);
        
        string subtitle = "blue: R = 2.0, green: R=2.5, orange: R = 3.0, purple: R = 3.5";
        g.titlin(const_cast<char*>(subtitle.c_str()), 3);
        
        int ic = g.intrgb(0.95, 0.95, 0.95);
        g.axsbgd(ic);
        
        g.graf(0.0, max_iterations, 0.0, max_iterations/10, 0.0, 1.0, 0.0, 0.1);
        
        g.setrgb(0.8, 0.8, 0.8);
        g.grid(3, 3);
        
        g.color("fore");
        g.height(50);
        g.title();
        
        g.thkcrv(3);

        vector<double> x_vals(max_iterations + 1);
        for (int i = 0; i <= max_iterations; i++) {
            x_vals[i] = i;
        }
        
        // R = 2.0 - синий
        g.setrgb(0.0, 0.0, 1.0);
        g.curve(x_vals.data(), success_rate2.data(), max_iterations + 1);
        
        // R = 2.5 - зеленый
        g.setrgb(0.0, 0.6, 0.0);
        g.curve(x_vals.data(), success_rate25.data(), max_iterations + 1);
        
        // R = 3.0 - оранжевый
        g.setrgb(1.0, 0.5, 0.0);
        g.curve(x_vals.data(), success_rate3.data(), max_iterations + 1);
        
        // R = 3.5 - фиолетовый
        g.setrgb(0.5, 0.0, 0.5);
        g.curve(x_vals.data(), success_rate35.data(), max_iterations + 1);
        
        g.disfin();
    }
};

int main() {
    CharacteristicsAnalyzer analyzer;
    
    cout << "select test type:" << endl;
    cout << "1. Hill" << endl;
    cout << "2. Shekel" << endl;
    cout << "choice: ";
    
    int choice;
    cin >> choice;
    
    vector<double> rate2 = (choice == 1) ? analyzer.TestHillR(2.0) : analyzer.TestShekelR(2.0);
    
    vector<double> rate25 = (choice == 1) ? analyzer.TestHillR(2.5) : analyzer.TestShekelR(2.5);
    
    vector<double> rate3 = (choice == 1) ? analyzer.TestHillR(3.0) : analyzer.TestShekelR(3.0);
    
    vector<double> rate35 = (choice == 1) ? analyzer.TestHillR(3.5) : analyzer.TestShekelR(3.5);
    
    string title = (choice == 1) ? "Hill Problem Family" : "Shekel Problem Family";
    analyzer.PlotCharacteristics(rate2, rate25, rate3, rate35, title);
    
    return 0;
}