#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

void qrDecomposition(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size();
    Q = vector<vector<double>>(n, vector<double>(n));
    R = vector<vector<double>>(n, vector<double>(n));
    vector<vector<double>> A_copy = A;

    for (int k = 0; k < n; ++k) {
        double norm = 0.0;
        for (int i = 0; i < n; ++i) {
            norm += A_copy[i][k] * A_copy[i][k];
        }
        norm = sqrt(norm);
        
        R[k][k] = norm;
        for (int i = 0; i < n; ++i) {
            Q[i][k] = A_copy[i][k] / norm;
        }
        
        for (int j = k + 1; j < n; ++j) {
            double dot = 0.0;
            for (int i = 0; i < n; ++i) {
                dot += Q[i][k] * A_copy[i][j];
            }
            R[k][j] = dot;
            for (int i = 0; i < n; ++i) {
                A_copy[i][j] -= dot * Q[i][k];
            }
        }
    }
}

void qrAlgorithm(vector<vector<double>>& A, double tol, int& iterations) {
    int n = A.size();
    vector<vector<double>> Q, R;
    iterations = 0;

    while (true) {
        iterations++;
        qrDecomposition(A, Q, R);

        vector<vector<double>> A_new(n, vector<double>(n));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                A_new[i][j] = 0.0;
                for (int k = 0; k < n; ++k) {
                    A_new[i][j] += R[i][k] * Q[k][j];
                }
            }
        }

        double max_off_diagonal = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                max_off_diagonal = max(max_off_diagonal, fabs(A_new[i][j]));
            }
        }

        A = A_new;
        if (max_off_diagonal < tol) break;
    }
}

int main() {
    vector<vector<double>> A = {
        {1, 7, -1},
        {-2, 2, -2},
        {9, -7, 3}
    };
    
    double tol = 1e-6;
    int iterations;

    qrAlgorithm(A, tol, iterations);

    cout << "Eigenvalues:" << endl;
    for (int i = 0; i < A.size(); ++i) {
        cout << setw(10) << A[i][i] << endl;
    }

    cout << "Number of iterations: " << iterations << endl;

    return 0;
}
