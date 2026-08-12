 # Multidimensional Global Optimization Problems and Parallel Methods for Their Solution

Performed by:
Kichanova K.K., student of group 3823B1FI3,
Lobachevsky National Research Nizhny Novgorod State University.

Scientific advisor:
Barkalov K.A., Professor of the Department of MOSA (Mathematical Optimization and System Analysis),
Lobachevsky State University of Nizhny Novgorod, IITMM (Institute of Information Technologies, Mathematics and Mechanics).

## Objective of the work:  
To study and implement methods and algorithms for the efficient solution of multidimensional global optimization problems, based on the information-statistical approach and the principles of parallel computing.


## Mathematical Statement of the Global Optimization Problem

In the general formulation, the global optimization problem is stated as follows. Let an objective function be given:

$f(y) : Q \to \mathbb{R}$

where $Q \subset \mathbb{R}^N$ is the feasible region in the N-dimensional Euclidean space. It is required to find a point of global minimum $y^* \in Q$ such that: 

$f(y^*) = \min\_{y \in Q} f(y)$

Here, the vector $y = (y_1, y_2, \dots, y_N)$ represents a set of controlled parameters, and the value $f(y^*)$ is typically computed based on a mathematical model of the object or process under study.

## Development Environment and Tools

- Programming language: C++ (C++17 standard)
- Compiler: g++ (Linux)
- Integrated development environment: Visual Studio Code
- Libraries: DISLIN C++ for visualization, GCGen for generating families of test functions (Hill, Shekel)

## How to compile and run tests:
1. mkdir build
2. cd build
3. cmake ..
4. make
5. cd ..
6. ./test_gcgen or ./characteristics or ./optimization
