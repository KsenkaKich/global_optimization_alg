#pragma once
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <limits>
#include <omp.h>

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
    std::function<double(double)> func;
};

class Solver {
protected:
    double eps;
    int Kmax;
    std::vector<Trial> Trials;
    Task task;
    Trial bestTrial;
    int iter;

public:
    Solver() : eps(0.001), Kmax(100) {}
    Solver(const Task& t, double eps_val, int kmax_val)
        : task(t), eps(eps_val), Kmax(kmax_val) {
    }
    Solver(const Solver& other)
        : eps(other.eps), Kmax(other.Kmax),
        Trials(other.Trials), task(other.task), bestTrial(other.bestTrial) {
    }
    virtual ~Solver() {};

    virtual void SetTask(Task t) { task = t; }
    virtual void SetEps(double eps_val) { eps = eps_val; }
    virtual void SetKmax(int kmax_val) { Kmax = kmax_val; }
    const std::vector<Trial>& GetTrials() const { return Trials; }
    virtual int GetTrialsCount() const { return Trials.size(); }
    virtual int GetIterations() const { return iter; }

    void Initialize() { 
        Trials.clear(); 
        iter = 0;
    }

    void FirstTrial() {
        Trial trial0, trial1;
        trial0.k = 1;
        trial0.x = task.a;
        trial0.z = task.func(task.a);

        trial1.k = 2;
        trial1.x = task.b;
        trial1.z = task.func(task.b);

        Trials.push_back(trial0);
        Trials.push_back(trial1);

        iter = 2;
    }

    bool CheckStopCondition(size_t t) {
        double interval = Trials[t + 1].x - Trials[t].x;
        double range = task.b - task.a;
        if (interval <= range * eps) return true;
        if (interval <= eps) return true;
        
        return false;
    }

    void InsertNewTrial(const Trial& newTrial) {
        Trials.push_back(newTrial);
        if (newTrial.z < bestTrial.z) {
            bestTrial = newTrial;
        }
    }

    virtual void Solve() {
        Initialize();
        FirstTrial();
        Trial trial;
        trial.k = 0;
        trial.x = (task.a + task.b) / 2.0;
        trial.z = task.func(trial.x);
        
        Trials.push_back(trial);
    }

    virtual Trial GetBest() {
        bestTrial = Trials[0];
        for (const auto& trial : Trials) {
            if (trial.z < bestTrial.z) {
                bestTrial = trial;
            }
        }
        return bestTrial;
    }
};

class GSASolver : public Solver {
private:
    double r;
    int p;
    
    double EstimateM() {
        size_t n = Trials.size() - 1;
        double M = 0.0;
        for (size_t i = 1; i <= n; i++) {
            double MInter = std::abs((Trials[i].z - Trials[i - 1].z) / (Trials[i].x - Trials[i - 1].x));
            M = std::max(M, MInter);
        }
        if (M > 0) return r * M;
        return 1.0;
    }

    std::vector<double> CalculateR(double m) {
        size_t n = Trials.size() - 1;
        std::vector<double> R(n);
        for (size_t i = 0; i < n; i++) {
            double delta_x = Trials[i + 1].x - Trials[i].x;
            double delta_z = Trials[i + 1].z - Trials[i].z;
            R[i] = m * delta_x + (delta_z * delta_z) / (m * delta_x) - 2.0 * (Trials[i + 1].z + Trials[i].z);
        }
        return R;
    }
    
    std::vector<size_t> FindTopIntervals(const std::vector<double>& R) {
        size_t n = R.size();
        std::vector<size_t> indices(n);
        for (size_t i = 0; i < n; i++) indices[i] = i;
        
        std::sort(indices.begin(), indices.end(), [&R](size_t a, size_t b) { return R[a] > R[b]; });
        
        int take = std::min(p, (int)n);
        return std::vector<size_t>(indices.begin(), indices.begin() + take);
    }

    double CalculateNewX(size_t t, double m) {
        return (Trials[t + 1].x + Trials[t].x) / 2.0 - (Trials[t + 1].z - Trials[t].z) / (2.0 * m);
    }

    Trial MakeNewTrial(double x, int k) {
        Trial newTrial;
        newTrial.k = k;
        newTrial.x = x;
        newTrial.z = task.func(x);
        return newTrial;
    }

public:
    GSASolver() : Solver(), r(2.0), p(1) {}
    
    void SetR(double r_val) { r = r_val; }
    void SetP(int p_val) { p = p_val; }

    void Solve() override {
        Initialize();
        FirstTrial();

        iter = 2;
        bool Stop = false;

        while (iter < Kmax && !Stop) {
            std::sort(Trials.begin(), Trials.end());
            double m = EstimateM();
            std::vector<double> R = CalculateR(m);
            std::vector<size_t> top = FindTopIntervals(R);
            
            if (top.empty()) break;
            
            for (size_t t : top) {
                if (CheckStopCondition(t)) {
                    Stop = true;
                    break;
                }
            }
            
            if (!Stop) {
                std::vector<Trial> new_trials(top.size());
                #pragma omp parallel for
                for (int i = 0; i < (int)top.size(); ++i) {
                    size_t t = top[i];
                    double x = CalculateNewX(t, m);
                    new_trials[i] = MakeNewTrial(x, iter + 1);
                }
                
                for (const Trial& nt : new_trials) {
                    InsertNewTrial(nt);
                }
                iter++;
            }
        }
    }
};

class ScanSolver : public Solver {
public:
    ScanSolver() : Solver() {}
    
    void Solve() override {
        Initialize();
        FirstTrial();

        iter = 2;
        bool Stop = false;

        while (iter < Kmax && !Stop) {
            std::sort(Trials.begin(), Trials.end());

            size_t t = FindLongestInterval();
            Stop = CheckStopCondition(t);

            if (!Stop) {
                double x = CalculateMiddlePoint(t);
                Trial newTrial = MakeNewTrial(x, iter + 1);
                InsertNewTrial(newTrial);
                iter++;
            }
        }
    }
    
    size_t FindLongestInterval() {
        size_t longest_index = 0;
        double max_length = 0.0;
        
        for (size_t i = 0; i < Trials.size() - 1; i++) {
            double length = Trials[i + 1].x - Trials[i].x;
            if (length > max_length) {
                max_length = length;
                longest_index = i;
            }
        }
        return longest_index;
    }
    
    double CalculateMiddlePoint(size_t t) {
        return 0.5 * (Trials[t].x + Trials[t + 1].x);
    }
    
    Trial MakeNewTrial(double x, int k) {
        Trial newTrial;
        newTrial.k = k;
        newTrial.x = x;
        newTrial.z = task.func(x);
        return newTrial;
    }
};