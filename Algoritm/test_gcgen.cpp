#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include "structs.h"
#include "visualization.h"

#include "Hill/HillProblem.hpp"
#include "Hill/HillProblemFamily.hpp"
#include "Shekel/ShekelProblem.hpp"
#include "Shekel/ShekelProblemFamily.hpp"

using namespace std;

struct TestResult {
    string name;
    int iterations;
    double best_x;
    double best_f;
    double true_x;
    double true_f;
    double error_x;
    double error_f;
    vector<Trial> trials;
    function<double(double)> func;
    double a, b;
};

void saveResultsToCSV(const vector<TestResult>& results, const string& solver_name, const string& filename) {
    bool file_exists = false;
    ifstream check(filename);
    if (check.good()) file_exists = true;
    check.close();
    
    ofstream csv(filename, ios::app);
    
    if (!file_exists) {
        csv << "Solver,Function,Iterations,BestX,TrueX,ErrorX,RelativeError%,BestF,TrueF,ErrorF\n";
    }
    
    for (const auto& r : results) {
        double range = (r.name.find("Hill") != string::npos) ? 1.0 : 10.0;
        double err_percent = (r.error_x / range) * 100.0;
        
        csv << solver_name << ","
            << r.name << ","
            << r.iterations << ","
            << r.best_x << ","
            << r.true_x << ","
            << r.error_x << ","
            << err_percent << ","
            << r.best_f << ","
            << r.true_f << ","
            << r.error_f << "\n";
    }
    csv.close();
}

void runHillTests(Solver& solver, const string& solver_name, int family_count, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    // Hill Problem
    THillProblem hill(0);
    vector<double> hill_opt = hill.GetOptimumPoint();
    double true_x = hill_opt[0], true_f = hill.GetOptimumValue();
    
    Task task;
    task.a = a;
    task.b = b;
    task.func = [&hill](double x) { return hill.ComputeFunction({x}); };
    
    solver.SetEps(eps);
    solver.SetKmax(Kmax);
    if (solver_name == "GSA") dynamic_cast<GSASolver&>(solver).SetR(2.0);
    solver.SetTask(task);
    solver.Solve();
    
    Trial best = solver.GetBest();
    TestResult res;
    res.name = "Hill_0";
    res.iterations = solver.GetTrials().size();
    res.best_x = best.x;
    res.best_f = best.z;
    res.true_x = true_x;
    res.true_f = true_f;
    res.error_x = abs(best.x - true_x);
    res.error_f = abs(best.z - true_f);
    res.trials = solver.GetTrials();
    res.func = task.func;
    res.a = a;
    res.b = b;
    results.push_back(res);
    
    // Hill Family
    THillProblemFamily hillFam;
    int hillFamilySize = hillFam.GetFamilySize();
    int maxHillFam = min(family_count, hillFamilySize);
    
    for (int i = 0; i < maxHillFam; i++) {
        THillProblem hillProb(i);
        vector<double> opt_point = hillProb.GetOptimumPoint();
        double true_x = opt_point[0], true_f = hillProb.GetOptimumValue();
        
        Task fam_task;
        fam_task.a = a;
        fam_task.b = b;
        fam_task.func = [&hillProb](double x) { return hillProb.ComputeFunction({x}); };
        
        solver.SetTask(fam_task);
        solver.Solve();
        
        best = solver.GetBest();
        res.name = "HillFamily_" + to_string(i);
        res.iterations = solver.GetTrials().size();
        res.best_x = best.x;
        res.best_f = best.z;
        res.true_x = true_x;
        res.true_f = true_f;
        res.error_x = abs(best.x - true_x);
        res.error_f = abs(best.z - true_f);
        res.trials = solver.GetTrials();
        res.func = fam_task.func;
        res.a = a;
        res.b = b;
        results.push_back(res);
    }
}

void runShekelTests(Solver& solver, const string& solver_name, int family_count, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    // Shekel Problem
    TShekelProblem shekel(0);
    vector<double> shekel_opt = shekel.GetOptimumPoint();
    double true_x = shekel_opt[0], true_f = shekel.GetOptimumValue();
    
    Task task;
    task.a = a;
    task.b = b;
    task.func = [&shekel](double x) { return shekel.ComputeFunction({x}); };
    
    solver.SetEps(eps);
    solver.SetKmax(Kmax);
    if (solver_name == "GSA") dynamic_cast<GSASolver&>(solver).SetR(2.0);
    solver.SetTask(task);
    solver.Solve();
    
    Trial best = solver.GetBest();
    TestResult res;
    res.name = "Shekel_0";
    res.iterations = solver.GetTrials().size();
    res.best_x = best.x;
    res.best_f = best.z;
    res.true_x = true_x;
    res.true_f = true_f;
    res.error_x = abs(best.x - true_x);
    res.error_f = abs(best.z - true_f);
    res.trials = solver.GetTrials();
    res.func = task.func;
    res.a = a;
    res.b = b;
    results.push_back(res);
    
    // Shekel Family
    TShekelProblemFamily shekelFam;
    int shekelFamilySize = shekelFam.GetFamilySize();
    int maxShekelFam = min(family_count, shekelFamilySize);
    
    for (int i = 0; i < maxShekelFam; i++) {
        TShekelProblem shekelProb(i);
        vector<double> opt_point = shekelProb.GetOptimumPoint();
        double true_x = opt_point[0], true_f = shekelProb.GetOptimumValue();
        
        Task fam_task;
        fam_task.a = a;
        fam_task.b = b;
        fam_task.func = [&shekelProb](double x) { return shekelProb.ComputeFunction({x}); };
        
        solver.SetTask(fam_task);
        solver.Solve();
        
        best = solver.GetBest();
        res.name = "ShekelFamily_" + to_string(i);
        res.iterations = solver.GetTrials().size();
        res.best_x = best.x;
        res.best_f = best.z;
        res.true_x = true_x;
        res.true_f = true_f;
        res.error_x = abs(best.x - true_x);
        res.error_f = abs(best.z - true_f);
        res.trials = solver.GetTrials();
        res.func = fam_task.func;
        res.a = a;
        res.b = b;
        results.push_back(res);
    }
}

void printResults(const vector<TestResult>& results, const string& name) {
    cout << "\n" << name << " results\n";
    
    cout << left << setw(25) << "Function"
         << setw(12) << "Iters"
         << setw(14) << "Best x"
         << setw(14) << "True x"
         << setw(12) << "Error x"
         << setw(14) << "Best f"
         << setw(14) << "True f"
         << endl;

    for (const auto& r : results) {
        cout << left << setw(25) << r.name
             << setw(12) << r.iterations
             << setw(14) << fixed << setprecision(6) << r.best_x
             << setw(14) << r.true_x
             << setw(12) << scientific << setprecision(4) << r.error_x
             << setw(14) << fixed << setprecision(6) << r.best_f
             << setw(14) << r.true_f
             << endl;
    }
}

void printStats(const vector<TestResult>& results, const string& name, double threshold) {
    double total_iter = 0;
    int success = 0, hill_success = 0, shekel_success = 0;
    int hill_count = 0, shekel_count = 0;
    
    for (const auto& r : results) {
        total_iter += r.iterations;
        
        bool is_hill = (r.name.find("Hill") != string::npos);
        double range = is_hill ? 1.0 : 10.0;
        double relative_error = r.error_x / range;
        
        if (is_hill) {
            hill_count++;
            if (relative_error < threshold) { success++; hill_success++; }
        } else {
            shekel_count++;
            if (relative_error < threshold) { success++; shekel_success++; }
        }
    }
    
    vector<int> iters;
    for (const auto& r : results) iters.push_back(r.iterations);
    sort(iters.begin(), iters.end());
    
    cout << "\n" << name << " STATISTICS\n";
    cout << "Total tests: " << results.size() << endl;
    cout << "Hill tests: " << hill_count << ", Shekel tests: " << shekel_count << endl;
    cout << "Successful (Hill): " << hill_success << " (" << (100.0 * hill_success / hill_count) << "%)" << endl;
    cout << "Successful (Shekel): " << shekel_success << " (" << (100.0 * shekel_success / shekel_count) << "%)" << endl;
    cout << "Successful (all): " << success << " (" << (100.0 * success / results.size()) << "%)" << endl;
    cout << "Average iterations: " << (total_iter / results.size()) << endl;
    cout << "Min iterations: " << iters[0] << endl;
    cout << "Max iterations: " << iters[iters.size() - 1] << endl;
}

void visualSolver(const vector<TestResult>& results, const string& solver_name){
    OptimizationPlotter plotter;
    const TestResult* best_hill = nullptr;
    const TestResult* worst_hill = nullptr;
    const TestResult* best_shekel = nullptr;
    const TestResult* worst_shekel = nullptr;

    for (const auto& r : results) {
        if (r.name.find("Hill") != string::npos) {
            if (best_hill == nullptr || r.iterations < best_hill->iterations) best_hill = &r;
            if (worst_hill == nullptr || r.iterations > worst_hill->iterations) worst_hill = &r;
        } else if (r.name.find("Shekel") != string::npos) {
            if (best_shekel == nullptr || r.iterations < best_shekel->iterations) best_shekel = &r;
            if (worst_shekel == nullptr || r.iterations > worst_shekel->iterations) worst_shekel = &r;
        }
    }
    
    int choice;
    do {
        cout << "\nvisualization " << solver_name << "\n";
        cout << "1. Best Hill case\n";
        cout << "2. Worst Hill case\n";
        cout << "3. Best Shekel case\n";
        cout << "4. Worst Shekel case\n";
        cout << "0. Back\n";
        cout << "choice: ";
        cin >> choice;
        
        Trial trial;
        switch(choice) {
            case 1:
                trial.x = best_hill->best_x;
                trial.z = best_hill->best_f;
                trial.k = best_hill->iterations;
                plotter.PlotAlgorithm(best_hill->trials, trial, best_hill->a, best_hill->b, solver_name + " - Best Hill: " + best_hill->name);
                break;
            case 2:
                trial.x = worst_hill->best_x;
                trial.z = worst_hill->best_f;
                trial.k = worst_hill->iterations;
                plotter.PlotAlgorithm(worst_hill->trials, trial, worst_hill->a, worst_hill->b, solver_name + " - Worst Hill: " + worst_hill->name);
                break;
            case 3:
                trial.x = best_shekel->best_x;
                trial.z = best_shekel->best_f;
                trial.k = best_shekel->iterations;
                plotter.PlotAlgorithm(best_shekel->trials, trial, best_shekel->a, best_shekel->b, solver_name + " - Best Shekel: " + best_shekel->name);
                break;
            case 4:
                trial.x = worst_shekel->best_x;
                trial.z = worst_shekel->best_f;
                trial.k = worst_shekel->iterations;
                plotter.PlotAlgorithm(worst_shekel->trials, trial, worst_shekel->a, worst_shekel->b, solver_name + " - Worst Shekel: " + worst_shekel->name);
                break;
            case 0:
                break;
            default:
                cout << "Invalid choice\n";
        }
        
    } while (choice != 0);
}

int main() {
    const int Kmax = 1000;
    const double eps = 1e-6;
    const double err_count = 0.0001;
    const int fam_count = 20;
    
    vector<TestResult> gsa_results, scan_results;
    
    GSASolver gsa_solver;
    
    runHillTests(gsa_solver, "GSA", fam_count, gsa_results, Kmax, eps, 0.0, 1.0);
    runShekelTests(gsa_solver, "GSA", fam_count, gsa_results, Kmax, eps, 0.0, 10.0);
    
    ScanSolver scan_solver;
    
    runHillTests(scan_solver, "SCAN", fam_count, scan_results, Kmax, eps, 0.0, 1.0);
    runShekelTests(scan_solver, "SCAN", fam_count, scan_results, Kmax, eps, 0.0, 10.0);
    
    printResults(gsa_results, "GSA");
    printResults(scan_results, "SCAN");
    
    printStats(gsa_results, "GSA", err_count);
    printStats(scan_results, "SCAN", err_count);
    
    saveResultsToCSV(gsa_results, "GSA", "gcgen_results.csv");
    saveResultsToCSV(scan_results, "SCAN", "gcgen_results.csv");

    int main_choice;
    do {
        cout << "\n1. Visualize GSA Solver results\n";
        cout << "2. Visualize SCAN Solver results\n";
        cout << "0. Exit\n";
        cout << "choice: ";
        cin >> main_choice;
        
        switch(main_choice) {
            case 1:
                visualSolver(gsa_results, "GSA");
                break;
            case 2:
                visualSolver(scan_results, "SCAN");
                break;
            case 0:
                break;
            default:
                cout << "Invalid choice\n";
        }
    } while (main_choice != 0);
    
    return 0;
}