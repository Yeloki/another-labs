#include <iostream>
#include <vector>

using namespace std;

void tridiagonalSolve(const vector<double> &a, const vector<double> &b,
                      const vector<double> &c, const vector<double> &d,
                      vector<double> &x, int n) {
  vector<double> P(n), Q(n);

  P[0] = -b[0] / a[0];
  Q[0] = d[0] / a[0];

  for (int i = 1; i < n; ++i) {
    double denominator = a[i] + c[i] * P[i - 1];
    P[i] = -b[i] / denominator;
    Q[i] = (d[i] - c[i] * Q[i - 1]) / denominator;
  }

  x[n - 1] = Q[n - 1];

  for (int i = n - 2; i >= 0; --i) {
    x[i] = P[i] * x[i + 1] + Q[i];
  }
}

int main() {
  int n = 5; // Number of equations

  // Non-zero elements of the tridiagonal matrix
  vector<double> a = {16, -16, 12, 12, 7};
  vector<double> b = {-8, 5, 3, -7, 0};
  vector<double> c = {0, -7, 4, -4, -1};

  // Right-hand side vector
  vector<double> d = {0, -123, -68, 104, 20};

  // Vector to store the solution
  vector<double> x(n);

  // Solve the system
  tridiagonalSolve(a, b, c, d, x, n);

  // Output the solution
  cout << "Solution to the system:" << endl;
  for (int i = 0; i < n; ++i) {
    cout << "x" << i + 1 << " = " << x[i] << endl;
  }

  return 0;
}
