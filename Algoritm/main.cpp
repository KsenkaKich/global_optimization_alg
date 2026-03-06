#include <iostream>
#include <vector>
#include "structs.h"
#include "visualization.h"

double testFunction1(double x) {
    return (x - 2.0) * (x - 2.0) + 1.0; //x = 2, f(x) = 1
}

double testFunction2(double x) {
    return x * x + cos(18 * x); // x = 0.17, f(x) = -0.97
}

double testFunction3(double x) {
    return exp(x) + sin(17 * x); // x = -4.9, f(x) = -0.9
}

int main() {
    Task t;
    t.a = -5.0;
    t.b = 5.0;
    t.func = testFunction3;

    OptimizationPlotter plotter;

    Solver* s;
    std::cout << "GSA Solver" << std::endl;
    s = new GSASolver();
    s->SetEps(0.001);
    dynamic_cast<GSASolver*>(s)->SetR(2.0);
    s->SetKmax(200);
    s->SetTask(t);
    s->Solve();
    Trial result1 = s->GetBest();

    std::cout << "\nOptimization result:" << std::endl;
    std::cout << "Best point: x* = " << result1.x << std::endl;
    std::cout << "Minimum value: f(x*) = " << result1.z << std::endl;

    plotter.PlotAlgorithm(s->GetTrials(), result1, t.a, t.b, "GSA Algorithm");

    std::cout << "\nScan Solver" << std::endl;
    s = new ScanSolver();
    s->SetEps(0.001);
    s->SetKmax(200);
    s->SetTask(t);
    s->Solve();
    Trial result2 = s->GetBest();

    std::cout << "\nOptimization result:" << std::endl;
    std::cout << "Best point: x* = " << result2.x << std::endl;
    std::cout << "Minimum value: f(x*) = " << result2.z << std::endl;

    plotter.PlotAlgorithm(s->GetTrials(), result2, t.a, t.b, "Scan Algorithm");

    return 0;
}
