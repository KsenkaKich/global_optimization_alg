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

    std::cout << "GSA Solver" << std::endl;
    {
        GSASolver gsa;
        gsa.SetEps(0.001);
        gsa.SetKmax(200);
        gsa.SetR(2.0);
        gsa.SetP(2);
        gsa.SetTask(t);
        gsa.Solve();

        Trial best = gsa.GetBest();
        std::cout << "Best point: x* = " << best.x << std::endl;
        std::cout << "Minimum: f(x*) = " << best.z << std::endl;
        std::cout << "Iterations:      " << gsa.GetIterations() << std::endl;
        std::cout << "Trials:          " << gsa.GetTrialsCount() << std::endl;

        plotter.PlotAlgorithm(gsa.GetTrials(), best, t.a, t.b, "GSA Algorithm");
    }

    std::cout << "\nScan Solver" << std::endl;
    {
        ScanSolver scan;
        scan.SetEps(0.001);
        scan.SetKmax(200);
        scan.SetTask(t);
        scan.Solve();

        Trial best = scan.GetBest();
        std::cout << "Best point:  x* = " << best.x << std::endl;
        std::cout << "Minimum:   f(x*) = " << best.z << std::endl;
        std::cout << "Iterations:      " << scan.GetIterations() << std::endl;

        plotter.PlotAlgorithm(scan.GetTrials(), best, t.a, t.b, "Scan Algorithm");
    }

    return 0;
}
