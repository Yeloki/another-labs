#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

void simpleIterationMethod(const vector<vector<double>>& A, const vector<double>& b, double tol, vector<double>& x, int& iterations) {
    int n = A.size();
    vector<double> x_old(n, 0.0);
    iterations = 0;
    
    while (true) {
        iterations++;
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x[i] -= A[i][j] * x_old[j];
                }
            }
            x[i] /= A[i][i];
        }
        
        double max_error = 0;
        for (int i = 0; i < n; ++i) {
            max_error = max(max_error, fabs(x[i] - x_old[i]));
        }
        
        if (max_error < tol) break;
        x_old = x;
    }
}

void gaussSeidelMethod(const vector<vector<double>>& A, const vector<double>& b, double tol, vector<double>& x, int& iterations) {
    int n = A.size();
    vector<double> x_old(n, 0.0);
    iterations = 0;
    
    while (true) {
        iterations++;
        for (int i = 0; i < n; ++i) {
            x[i] = b[i];
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    x[i] -= A[i][j] * x[j];
                }
            }
            x[i] /= A[i][i];
        }
        
        double max_error = 0;
        for (int i = 0; i < n; ++i) {
            max_error = max(max_error, fabs(x[i] - x_old[i]));
        }
        
        if (max_error < tol) break;
        x_old = x;
    }
}

int main() {
    vector<vector<double>> A = {
        {-14, 6, 1, -5},
        {-6, 27, 7, -6},
        {7, -5, -23, -8},
        {3, -8, -7, 26}
    };
    vector<double> b = {95, -41, 69, 27};

    double tol = 1e-6;
    vector<double> x(A.size(), 0.0);
    int iterations;

    cout << "Method of Simple Iterations: " << endl;
    simpleIterationMethod(A, b, tol, x, iterations);
    for (int i = 0; i < x.size(); ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "Number of iterations: " << iterations << endl;

    cout << endl;

    x.assign(A.size(), 0.0);
    cout << "Gauss-Seidel Method: " << endl;
    gaussSeidelMethod(A, b, tol, x, iterations);
    for (int i = 0; i < x.size(); ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    cout << "Number of iterations: " << iterations << endl;

    return 0;
}
