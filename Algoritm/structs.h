#pragma once
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>

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
    std::ofstream outputFile; 

public:
    Solver() : eps(0.001), Kmax(100) {}
    Solver(const Task& t, double eps_val, int kmax_val)
        : task(t), eps(eps_val), Kmax(kmax_val) {
    }
    Solver(const Solver& other)
        : eps(other.eps), Kmax(other.Kmax),
        Trials(other.Trials), task(other.task), bestTrial(other.bestTrial) {
    }
    virtual ~Solver() {
        if (outputFile.is_open()) {
            outputFile.flush(); 
            outputFile.close();
        }
    };

    virtual void SetTask(Task t) { task = t; }
    virtual void SetEps(double eps_val) { eps = eps_val; }
    virtual void SetKmax(int kmax_val) { Kmax = kmax_val; }

    void Initialize() { Trials.clear(); }

    void FirstTrial() {
        Trial trial0, trial1;
        trial0.k = 1;
        trial0.x = task.a;
        trial0.z = task.func(task.a);
        std::cout << "k = " << trial0.k << " x = " << trial0.x << std::endl;

        trial1.k = 2;
        trial1.x = task.b;
        trial1.z = task.func(task.b);
        std::cout << "k = " << trial1.k << " x = " << trial1.x << std::endl;

        Trials.push_back(trial0);
        Trials.push_back(trial1);
    }

    bool CheckStopCondition(size_t t) {
        return (Trials[t + 1].x - Trials[t].x) <= eps;
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

    void OpenOutputFile(const std::string& filename) {
        outputFile.open(filename);
        if (!outputFile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }
        outputFile << "Iteration,x,f(x)\n";
        outputFile.flush();
    }

    void WriteTrialToFile(const Trial& trial) {
        if (outputFile.is_open()) {
            outputFile << std::fixed << std::setprecision(15);
            outputFile << trial.k << "," << trial.x << "," << trial.z << "\n";
            outputFile.flush();
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

class GSASolver : public Solver{
private:
    double r;
public:
    GSASolver(): Solver(), r(2.0) {}
    void SetR(double r_val) { r = r_val; }

    void Solve() override {
        Initialize();
        FirstTrial();

        OpenOutputFile("gsa_results.csv");
        WriteTrialToFile(Trials[0]);
        WriteTrialToFile(Trials[1]);

        int k = 2;
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
                WriteTrialToFile(newTrial);
                std::cout << "k = " << k << " x = " << newTrial.x << std::endl;
            }
        }
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
};

class ScanSolver : public Solver{
public:
    ScanSolver() : Solver() {}
    void Solve() override {
        
        Initialize();
        FirstTrial();

        OpenOutputFile("scan_results.csv");
        WriteTrialToFile(Trials[0]);
        WriteTrialToFile(Trials[1]);
        
        int k = 2;
        bool Stop = false;

        while (k < Kmax && !Stop) {
            size_t t = FindLongestInterval();
            double interval_length = Trials[t + 1].x - Trials[t].x;

            Stop = CheckStopCondition(t);

            if (!Stop) {
                double x = CalculateMiddlePoint(t);
                Trial newTrial = MakeNewTrial(x, k + 1);
                
                if (!CheckPointExists(newTrial.x)) {
                    InsertNewTrial(newTrial, t);
                    k++;
                    
                    WriteTrialToFile(newTrial);
                    std::cout << "k = " << k << " x = " << x << std::endl;
                }
            }
        }
    }
    size_t FindLongestInterval() {
        if (Trials.size() < 2) {
            return 0;
        }
        
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