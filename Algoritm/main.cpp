#include <iostream>
#include <vector>
#include "structs.h"

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

    Solver s;
    s.SetEps(0.001);
    s.SetR(2.0);
    s.SetKmax(200);
    s.SetTask(t);

    s.Solve();
    s.GetBest();

    return 0;
}