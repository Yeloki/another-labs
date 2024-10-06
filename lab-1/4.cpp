#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void printMatrix(const vector<vector<double>> &A) {
  for (const auto &row : A) {
    for (const auto &val : row) {
      cout << setw(10) << val << " ";
    }
    cout << endl;
  }
}

void jacobiMethod(vector<vector<double>> &A,
                  vector<vector<double>> &eigenvectors, double tol,
                  int &iterations) {
  int n = A.size();
  eigenvectors = vector<vector<double>>(n, vector<double>(n, 0.0));
  for (int i = 0; i < n; ++i) {
    eigenvectors[i][i] = 1.0;
  }

  iterations = 0;
  double max_off_diagonal;

  do {
    iterations++;
    max_off_diagonal = 0;
    int p, q;

    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (fabs(A[i][j]) > max_off_diagonal) {
          max_off_diagonal = fabs(A[i][j]);
          p = i;
          q = j;
        }
      }
    }

    if (max_off_diagonal < tol)
      break;

    double theta = 0.5 * atan2(2.0 * A[p][q], A[q][q] - A[p][p]);
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    vector<vector<double>> R = eigenvectors;

    for (int i = 0; i < n; ++i) {
      R[i][p] = eigenvectors[i][p] * cos_theta - eigenvectors[i][q] * sin_theta;
      R[i][q] = eigenvectors[i][p] * sin_theta + eigenvectors[i][q] * cos_theta;
    }

    double a_pp = A[p][p];
    double a_qq = A[q][q];
    double a_pq = A[p][q];

    A[p][p] = cos_theta * cos_theta * a_pp -
              2.0 * cos_theta * sin_theta * a_pq + sin_theta * sin_theta * a_qq;
    A[q][q] = sin_theta * sin_theta * a_pp +
              2.0 * cos_theta * sin_theta * a_pq + cos_theta * cos_theta * a_qq;
    A[p][q] = 0.0;
    A[q][p] = 0.0;

    for (int i = 0; i < n; ++i) {
      if (i != p && i != q) {
        double a_ip = A[i][p];
        double a_iq = A[i][q];
        A[i][p] = A[p][i] = cos_theta * a_ip - sin_theta * a_iq;
        A[i][q] = A[q][i] = cos_theta * a_iq + sin_theta * a_ip;
      }
    }

    eigenvectors = R;
  } while (max_off_diagonal > tol);
}

int main() {
  vector<vector<double>> A = {{-3, -1, 3}, {-1, 8, 1}, {3, 1, 5}};

  vector<vector<double>> eigenvectors;
  double tol = 1e-6;
  int iterations;

  jacobiMethod(A, eigenvectors, tol, iterations);

  cout << "Eigenvalues:" << endl;
  for (int i = 0; i < A.size(); ++i) {
    cout << setw(10) << A[i][i] << endl;
  }

  cout << endl;
  cout << "Eigenvectors:" << endl;
  printMatrix(eigenvectors);

  cout << "Number of iterations: " << iterations << endl;

  return 0;
}
