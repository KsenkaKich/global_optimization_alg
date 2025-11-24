#pragma once
#include <cmath>
#include <algorithm>
#include <functional>

struct Trial {
    double x;
    int k;
    double z;

    bool operator<(const Trial& other) const {
        return x < other.x;
    }
};

struct Task {
    double a, b;
    double (*func)(double);
};

class Solver {
private:
    double r;
    double eps;
    int Kmax;
    std::vector<Trial> Trials;
    Task task;
    Trial bestTrial;
public:
    Solver() : r(2.0), eps(0.001), Kmax(100) {}
    Solver(const Task& t, double r_val, double eps_val, int kmax_val)
        : task(t), r(r_val), eps(eps_val), Kmax(kmax_val) {
    }
    Solver(const Solver& other)
        : r(other.r), eps(other.eps), Kmax(other.Kmax),
        Trials(other.Trials), task(other.task), bestTrial(other.bestTrial) {
    }

    void SetTask(Task t) { task = t; }
    void SetR(double r_val) { r = r_val; }
    void SetEps(double eps_val) { eps = eps_val; }
    void SetKmax(int kmax_val) { Kmax = kmax_val; }

    // Functions to simplify Solve()
    void Initialize() { Trials.clear(); }

    void FirstTrial() {
        Trial trial0, trial1;
        trial0.k = 0;
        trial0.x = task.a;
        trial0.z = task.func(task.a);

        trial1.k = 1;
        trial1.x = task.b;
        trial1.z = task.func(task.b);

        Trials.push_back(trial0);
        Trials.push_back(trial1);
    }

    double EstimateM() {
        size_t n = Trials.size() - 1;
        double M = 0.0;
        for (size_t i = 1; i <= n; i++) {
            double MInter = std::abs((Trials[i].z - Trials[i - 1].z) /
                (Trials[i].x - Trials[i - 1].x));
            M = std::max(M, MInter);
        }
        if (M > 0) {
            return r * M;
        }
        else {
            return 1.0;
        }
    }

    std::vector<double> CalculateR(double m) {
        size_t n = Trials.size() - 1;
        std::vector<double> R(n);
        for (size_t i = 0; i < n; i++) {
            double delta_x = Trials[i + 1].x - Trials[i].x;
            double delta_z = Trials[i + 1].z - Trials[i].z;

            R[i] = m * delta_x + (delta_z * delta_z) / (m * delta_x)
                - 2.0 * (Trials[i + 1].z + Trials[i].z);
        }
        return R;
    }

    size_t FindMaxR(const std::vector<double>& R) {
        size_t n = Trials.size() - 1;
        size_t t = 0;
        double maxR = R[0];
        for (size_t i = 1; i < n; i++) {
            if (R[i] > maxR) {
                maxR = R[i];
                t = i;
            }
        }

        for (size_t i = 0; i < n; i++) {
            if (R[i] == maxR) {
                t = i;
                break;
            }
        }
        return t;
    }

    bool CheckStopCondition(size_t t) {
        return (Trials[t + 1].x - Trials[t].x) <= eps;
    }

    double CalculateNewX(size_t t, double m) {
        return (Trials[t + 1].x + Trials[t].x) / 2.0
            - (Trials[t + 1].z - Trials[t].z) / (2.0 * m);
    }

    Trial MakeNewTrial(double x, int k) {
        Trial newTrial;
        newTrial.k = k;
        newTrial.x = x;
        newTrial.z = task.func(x);
        return newTrial;
    }

    bool CheckPointExists(double x) {
        for (const auto& trial : Trials) {
            if (std::abs(trial.x - x) < 1e-15) {
                std::cout << "The point already exists" << std::endl;
                return true;
            }
        }
        return false;
    }

    void InsertNewTrial(const Trial& newTrial, size_t t) {
        Trials.insert(Trials.begin() + t + 1, newTrial);
    }

    void Solve() {
        Initialize();
        FirstTrial();

        int k = 1;
        bool Stop = false;

        while (k < Kmax && !Stop) {
            double m = EstimateM();
            std::vector<double> R = CalculateR(m);
            size_t t = FindMaxR(R);
            Stop = CheckStopCondition(t);
            double x = CalculateNewX(t, m);
            Trial newTrial = MakeNewTrial(x, k + 1);
            if (!CheckPointExists(newTrial.x)) {
                InsertNewTrial(newTrial, t);
                k++;
                std::cout << "k = " << k - 2 << " x = " << newTrial.x << std::endl;
            }
        }
    }

    void GetBest() {
        bestTrial = Trials[0];
        for (const auto& trial : Trials) {
            if (trial.z < bestTrial.z) {
                bestTrial = trial;
            }
        }

        std::cout << "\nOptimization result:" << std::endl;
        std::cout << "Best point: x* = " << bestTrial.x << std::endl;
        std::cout << "Minimum value: f(x*) = " << bestTrial.z << std::endl;
    }
};
