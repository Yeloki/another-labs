#include <cmath>
#include <iostream>

using namespace std;

// Значение параметра a
const double a = 4.0;

// Система уравнений
// f1 = 0
double f1(double x1, double x2) {
  return (x1 * x1) / (a * a) + (x2 * x2) / ((a / 2) * (a / 2)) - 1;
}

// f2 = 0
double f2(double x1, double x2) { return a * x2 - exp(x1) - x1; }

// Частные производные для метода Ньютона
double df1dx1(double x1, double x2) { return 2 * x1 / (a * a); }

double df1dx2(double x1, double x2) { return 2 * x2 / ((a / 2) * (a / 2)); }

double df2dx1(double x1, double x2) { return -exp(x1) - 1; }

double df2dx2(double x1, double x2) { return a; }

// Метод простой итерации
void simpleIterationMethod(double x1_initial, double x2_initial, double tol,
                           int &iterations, double &x1, double &x2) {
  x1 = x1_initial;
  x2 = x2_initial;

  double x1_new, x2_new;

  iterations = 0;
  while (true) {
    iterations++;
    x1_new = sqrt(a * a * (1 - (x2 * x2) / ((a / 2) * (a / 2))));
    x2_new = (exp(x1) + x1) / a;

    if (fabs(x1_new - x1) < tol && fabs(x2_new - x2) < tol)
      break;

    x1 = x1_new;
    x2 = x2_new;
  }
}

// Метод Ньютона
void newtonMethod(double x1_initial, double x2_initial, double tol,
                  int &iterations, double &x1, double &x2) {
  x1 = x1_initial;
  x2 = x2_initial;

  iterations = 0;

  while (true) {
    double J[2][2] = {{df1dx1(x1, x2), df1dx2(x1, x2)},
                      {df2dx1(x1, x2), df2dx2(x1, x2)}};

    double F[2] = {-f1(x1, x2), -f2(x1, x2)};

    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    double dx1 = (F[0] * J[1][1] - F[1] * J[0][1]) / detJ;
    double dx2 = (J[0][0] * F[1] - J[1][0] * F[0]) / detJ;

    x1 += dx1;
    x2 += dx2;
    iterations++;

    if (fabs(dx1) < tol && fabs(dx2) < tol)
      break;
  }
}

int main() {
  double x1_initial = 2.0; // Начальное приближение для x1
  double x2_initial = 2.0; // Начальное приближение для x2
  double tol = 1e-6;       // Заданная точность
  int iterations;
  double x1, x2;

  // Метод простой итерации - решение
  simpleIterationMethod(x1_initial, x2_initial, tol, iterations, x1, x2);
  cout << "Метод простой итерации:" << endl;
  cout << "Решение: x1 = " << x1 << ", x2 = " << x2 << endl;
  cout << "Число итераций: " << iterations << endl;

  // Метод Ньютона - решение
  newtonMethod(x1_initial, x2_initial, tol, iterations, x1, x2);
  cout << "Метод Ньютона:" << endl;
  cout << "Решение: x1 = " << x1 << ", x2 = " << x2 << endl;
  cout << "Число итераций: " << iterations << endl;

  return 0;
}
