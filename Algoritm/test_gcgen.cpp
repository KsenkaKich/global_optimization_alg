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

void runSingleHillTest(Solver& solver, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    THillProblem hill(0);
    vector<double> hill_opt = hill.GetOptimumPoint();
    double true_x = hill_opt[0], true_f = hill.GetOptimumValue();
    
    Task task;
    task.a = a;
    task.b = b;
    task.func = [&hill](double x) { return hill.ComputeFunction({x}); };
    
    solver.SetEps(eps);
    solver.SetKmax(Kmax);
    
    GSASolver* gsa = dynamic_cast<GSASolver*>(&solver);
    if (gsa) gsa->SetR(2.0);
    
    solver.SetTask(task);
    solver.Solve();
    
    Trial best = solver.GetBest();
    TestResult res;
    res.name = "Hill";
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
}

void runHillFamilyTests(Solver& solver, int family_count, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    THillProblemFamily hillFam;
    int hillFamilySize = hillFam.GetFamilySize();
    int HillFam = min(family_count, hillFamilySize);
    solver.SetEps(eps);
    solver.SetKmax(Kmax); 
    
    for (int i = 0; i < HillFam; i++) {
        THillProblem hillProb(i);
        vector<double> opt_point = hillProb.GetOptimumPoint();
        double true_x = opt_point[0], true_f = hillProb.GetOptimumValue();
        
        Task fam_task;
        fam_task.a = a;
        fam_task.b = b;
        fam_task.func = [&hillProb](double x) { return hillProb.ComputeFunction({x}); };
        
        solver.SetTask(fam_task);
        solver.Solve();
        
        Trial best = solver.GetBest();
        TestResult res;
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

void runSingleShekelTest(Solver& solver, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    TShekelProblem shekel(0);
    vector<double> shekel_opt = shekel.GetOptimumPoint();
    double true_x = shekel_opt[0], true_f = shekel.GetOptimumValue();
    
    Task task;
    task.a = a;
    task.b = b;
    task.func = [&shekel](double x) { return shekel.ComputeFunction({x}); };
    
    solver.SetEps(eps);
    solver.SetKmax(Kmax);
    
    GSASolver* gsa = dynamic_cast<GSASolver*>(&solver);
    if (gsa) gsa->SetR(2.0);
    
    solver.SetTask(task);
    solver.Solve();
    
    Trial best = solver.GetBest();
    TestResult res;
    res.name = "Shekel";
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
}

void runShekelFamilyTests(Solver& solver, int family_count, vector<TestResult>& results, int Kmax, double eps, double a, double b) {
    TShekelProblemFamily shekelFam;
    int shekelFamilySize = shekelFam.GetFamilySize();
    int ShekelFam = min(family_count, shekelFamilySize);
    solver.SetEps(eps);
    solver.SetKmax(Kmax); 
    
    for (int i = 0; i < ShekelFam; i++) {
        TShekelProblem shekelProb(i);
        vector<double> opt_point = shekelProb.GetOptimumPoint();
        double true_x = opt_point[0], true_f = shekelProb.GetOptimumValue();
        
        Task fam_task;
        fam_task.a = a;
        fam_task.b = b;
        fam_task.func = [&shekelProb](double x) { return shekelProb.ComputeFunction({x}); };
        
        solver.SetTask(fam_task);
        solver.Solve();
        
        Trial best = solver.GetBest();
        TestResult res;
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
             << setw(14) << setprecision(6) << r.best_x
             << setw(14) << r.true_x
             << setw(12) << setprecision(4) << r.error_x
             << setw(14) << setprecision(6) << r.best_f
             << setw(14) << r.true_f
             << endl;
    }
}

void printStatistics(const vector<TestResult>& results, const string& name, double error_count) {
    double total_iter = 0;
    int success = 0;
    int hill_success = 0, shekel_success = 0;
    int hill_count = 0, shekel_count = 0;
    double min_iter = 1e9, max_iter = 0;
    double min_error = 1e9, max_error = 0;
    double total_error = 0;
    
    for (const auto& r : results) {
        total_iter += r.iterations;
        total_error += r.error_x;
        
        min_iter = min(min_iter, (double)r.iterations);
        max_iter = max(max_iter, (double)r.iterations);
        min_error = min(min_error, r.error_x);
        max_error = max(max_error, r.error_x);
        
        bool is_hill = (r.name.find("Hill") != string::npos);
        double range = is_hill ? 1.0 : 10.0;
        double relative_error = r.error_x / range;
        
        if (is_hill) {
            hill_count++;
            if (relative_error < error_count) { success++; hill_success++; }
        } else {
            shekel_count++;
            if (relative_error < error_count) { success++; shekel_success++; }
        }
    }
    
    double avg_iter = total_iter / results.size();
    double avg_error = total_error / results.size();
    
    cout << "\nSTATISTICS FOR: " << name << "\n";
    
    if (hill_count > 0) cout << "Hill tests: " << hill_count << endl;
    if (shekel_count > 0) cout << "Shekel tests: " << shekel_count << endl;
    
    if (hill_count > 0) 
        cout << "Hill success: " << hill_success << "/" << hill_count << " (" << (100.0 * hill_success / hill_count) << "%)\n";
    if (shekel_count > 0)
        cout << "Shekel success: " << shekel_success << "/" << shekel_count  << " (" << (100.0 * shekel_success / shekel_count) << "%)\n";
    
    cout << "\niterations\n";
    cout << "average: " << avg_iter << endl;
    cout << "min: " << min_iter << endl;
    cout << "max: " << max_iter << endl;
    
    cout << "\nerror in x\n";
    cout << "average: " << avg_error << endl;
    cout << "min: " << min_error << endl;
    cout << "max: " << max_error << endl;
}

void visualSolver(const vector<TestResult>& results, int test_choice) {
    OptimizationPlotter plotter;

    if (test_choice == 1 || test_choice == 2) {
        const TestResult& res = results[0];
        Trial trial;
        trial.x = res.best_x;
        trial.z = res.best_f;
        trial.k = res.iterations;
        plotter.PlotAlgorithm(res.trials, trial, res.a, res.b, res.name);
    }
    else if (test_choice == 3 || test_choice == 4) {
        const TestResult* best_result = &results[0];
        const TestResult* worst_result = &results[0];
        
        for (const auto& r : results) {
            if (r.iterations < best_result->iterations) best_result = &r;
            if (r.iterations > worst_result->iterations) worst_result = &r;
        }
        
        int choice;
        do {
            cout << "\nvisualization\n";
            cout << "1. best result (iterations: " << best_result->iterations << ")\n";
            cout << "2. worst result (iterations: " << worst_result->iterations << ")\n";
            cout << "0. back\n";
            cout << "choice: ";
            cin >> choice;
            
            Trial trial;
            switch(choice) {
                case 1:
                    trial.x = best_result->best_x;
                    trial.z = best_result->best_f;
                    trial.k = best_result->iterations;
                    plotter.PlotAlgorithm(best_result->trials, trial, best_result->a, best_result->b, "Best result: " + best_result->name);
                    break;
                case 2:
                    trial.x = worst_result->best_x;
                    trial.z = worst_result->best_f;
                    trial.k = worst_result->iterations;
                    plotter.PlotAlgorithm(worst_result->trials, trial, worst_result->a, worst_result->b, "Worst result: " + worst_result->name);
                    break;
                case 0:
                    break;
                default:
                    cout << "invalid choice\n";
            }
        } while (choice != 0);
    }
}

int main() {
    const int Kmax = 1000;
    const double eps = 1e-6;
    const double err_count = 0.0001;
    const int fam_count = 20;
    
    vector<TestResult> gsa_results, scan_results;
    
    int test_choice;
    cout << "1. Hill\n";
    cout << "2. Shekel\n";
    cout << "3. HillFamily\n";
    cout << "4. ShekelFamily \n";
    cout << "choice: ";
    cin >> test_choice;
    
    GSASolver gsa_solver;
    ScanSolver scan_solver;
    
    if (test_choice == 1) {
        runSingleHillTest(gsa_solver, gsa_results, Kmax, eps, 0.0, 1.0);
        runSingleHillTest(scan_solver, scan_results, Kmax, eps, 0.0, 1.0);
    }
    else if (test_choice == 2) {
        runSingleShekelTest(gsa_solver, gsa_results, Kmax, eps, 0.0, 10.0);
        runSingleShekelTest(scan_solver, scan_results, Kmax, eps, 0.0, 10.0);
    }
    else if (test_choice == 3) {
        runHillFamilyTests(gsa_solver, fam_count, gsa_results, Kmax, eps, 0.0, 1.0);
        runHillFamilyTests(scan_solver, fam_count, scan_results, Kmax, eps, 0.0, 1.0);
    }
    else if (test_choice == 4) {
        runShekelFamilyTests(gsa_solver, fam_count, gsa_results, Kmax, eps, 0.0, 10.0);
        runShekelFamilyTests(scan_solver, fam_count, scan_results, Kmax, eps, 0.0, 10.0);
    }
    
    printResults(gsa_results, "GSA");
    printResults(scan_results, "SCAN");
    
    if (test_choice == 1) {
        printStatistics(gsa_results, "GSA - Hill", err_count);
        printStatistics(scan_results, "SCAN - Hill", err_count);
    }
    else if (test_choice == 2) {
        printStatistics(gsa_results, "GSA - Shekel", err_count);
        printStatistics(scan_results, "SCAN - Shekel", err_count);
    }
    else if (test_choice == 3) {
        printStatistics(gsa_results, "GSA - HillFamily", err_count);
        printStatistics(scan_results, "SCAN - HillFamily", err_count);
    }
    else if (test_choice == 4) {
        printStatistics(gsa_results, "GSA - ShekelFamily", err_count);
        printStatistics(scan_results, "SCAN - ShekelFamily", err_count);
    }    
    
    int show_plot;
    cout << "\nshow plots? (1 - Yes, 0 - No): ";
    cin >> show_plot;
    
    if (show_plot == 1) {
        int main_choice;
        do {
            cout << "\n1. visualize GSA Solver results\n";
            cout << "2. visualize SCAN Solver results\n";
            cout << "0. exit\n";
            cout << "choice: ";
            cin >> main_choice;
            
            switch(main_choice) {
                case 1:
                    visualSolver(gsa_results, test_choice);
                    break;
                case 2:
                    visualSolver(scan_results, test_choice);
                    break;
                case 0:
                    break;
                default:
                    cout << "invalid choice\n";
            }
        } while (main_choice != 0);
    }
    
    return 0;
}