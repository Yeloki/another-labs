#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void luDecomposition(vector<vector<double>> &A, int n, vector<int> &P) {
  for (int i = 0; i < n; ++i) {
    P[i] = i;
  }

  for (int i = 0; i < n; ++i) {
    double maxA = 0.0;
    int imax = i;

    for (int k = i; k < n; ++k) {
      if (fabs(A[k][i]) > maxA) {
        maxA = fabs(A[k][i]);
        imax = k;
      }
    }

    if (imax != i) {
      swap(P[i], P[imax]);
      for (int j = 0; j < n; ++j) {
        swap(A[i][j], A[imax][j]);
      }
    }

    for (int j = i + 1; j < n; ++j) {
      A[j][i] /= A[i][i];
      for (int k = i + 1; k < n; ++k) {
        A[j][k] -= A[j][i] * A[i][k];
      }
    }
  }
}

void luSolve(const vector<vector<double>> &LU, const vector<int> &P,
             const vector<double> &b, vector<double> &x, int n) {
  vector<double> y(n);

  for (int i = 0; i < n; ++i) {
    y[i] = b[P[i]];
    for (int j = 0; j < i; ++j) {
      y[i] -= LU[i][j] * y[j];
    }
  }

  for (int i = n - 1; i >= 0; --i) {
    x[i] = y[i];
    for (int j = i + 1; j < n; ++j) {
      x[i] -= LU[i][j] * x[j];
    }
    x[i] /= LU[i][n - 1];
  }
}

double luDeterminant(const vector<vector<double>> &LU, const vector<int> &P,
                     int n) {
  double det = 1.0;
  for (int i = 0; i < n; ++i) {
    det *= LU[i][i];
    if (P[i] != i) {
      det = -det;
    }
  }
  return det;
}

vector<vector<double>> luInverse(const vector<vector<double>> &LU,
                                 const vector<int> &P, int n) {
  vector<vector<double>> inverse(n, vector<double>(n));
  vector<double> b(n), x(n);

  for (int i = 0; i < n; ++i) {
    fill(b.begin(), b.end(), 0.0);
    b[i] = 1.0;
    luSolve(LU, P, b, x, n);
    for (int j = 0; j < n; ++j) {
      inverse[j][i] = x[j];
    }
  }

  return inverse;
}

int main() {
  const int n = 4;
  vector<vector<double>> A = {
      {-9, 8, 8, 6}, {-7, -9, 5, 4}, {-3, -1, 8, -1}, {3, -1, -4, -5}};
  vector<double> b = {-81, -50, -69, 48};

  vector<int> P(n);
  luDecomposition(A, n, P);

  vector<double> x(n);
  luSolve(A, P, b, x, n);

  cout << "Solution to the system:" << endl;
  for (int i = 0; i < n; ++i) {
    cout << "x" << i + 1 << " = " << x[i] << endl;
  }

  double det = luDeterminant(A, P, n);
  cout << "Determinant of the matrix: " << det << endl;

  vector<vector<double>> invA = luInverse(A, P, n);
  cout << "Inverse matrix:" << endl;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      cout << setw(12) << invA[i][j] << " ";
    }
    cout << endl;
  }

  return 0;
}
