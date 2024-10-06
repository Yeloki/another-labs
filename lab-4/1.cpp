#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// Точное решение
double exactSolution(double x) { return 1.0 + log(fabs(x)); }

// Первая производная точного решения
double exactSolutionDerivative(double x) { return 1.0 / x; }

// Метод Эйлера
void eulerMethod(vector<double> &x, vector<double> &y, vector<double> &yp,
                 double h) {
  int n = x.size();
  for (int i = 1; i < n; ++i) {
    double x_prev = x[i - 1];
    double y_prev = y[i - 1];
    double yp_prev = yp[i - 1];
    // y'' = -(1/x) * y'
    double ypp = -(1.0 / x_prev) * yp_prev;
    y[i] = y_prev + h * yp_prev;
    yp[i] = yp_prev + h * ypp;
  }
}

// Метод Рунге-Кутты 4-го порядка
void rungeKuttaMethod(vector<double> &x, vector<double> &y, vector<double> &yp,
                      double h) {
  int n = x.size();
  for (int i = 1; i < n; ++i) {
    double xi = x[i - 1];
    double yi = y[i - 1];
    double ypi = yp[i - 1];

    double k1y = ypi;
    double k1yp = -(1.0 / xi) * ypi;

    double k2y = ypi + 0.5 * h * k1yp;
    double k2yp = -(1.0 / (xi + 0.5 * h)) * (ypi + 0.5 * h * k1yp);

    double k3y = ypi + 0.5 * h * k2yp;
    double k3yp = -(1.0 / (xi + 0.5 * h)) * (ypi + 0.5 * h * k2yp);

    double k4y = ypi + h * k3yp;
    double k4yp = -(1.0 / (xi + h)) * (ypi + h * k3yp);

    y[i] = yi + h * (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;
    yp[i] = ypi + h * (k1yp + 2 * k2yp + 2 * k3yp + k4yp) / 6.0;
  }
}

// Метод Адамса 4-го порядка
void adamsMethod(vector<double> &x, vector<double> &y, vector<double> &yp,
                 double h) {
  int n = x.size();
  vector<double> ypp(n, 0);
  // Заполняем начальные данных используя метод Рунге - Кутты 4-го порядка
  rungeKuttaMethod(x, y, yp, h);
  for (int i = 0; i < 4; ++i) {
    ypp[i] = -(1.0 / x[i]) * yp[i];
  }
  for (int i = 4; i < n; ++i) {
    y[i] =
        y[i - 1] +
        h * (55 * yp[i - 1] - 59 * yp[i - 2] + 37 * yp[i - 3] - 9 * yp[i - 4]) /
            24.0;
    yp[i] = yp[i - 1] + h *
                            (55 * ypp[i - 1] - 59 * ypp[i - 2] +
                             37 * ypp[i - 3] - 9 * ypp[i - 4]) /
                            24.0;
    ypp[i] = -(1.0 / x[i]) * yp[i];
  }
}

// Метод Рунге-Ромберга для оценки погрешности
double rungeRomberg(double I_h1, double I_h2, int p) {
  return (I_h2 - I_h1) / (pow(2, p) - 1);
}

int main() {
  double x0 = 1.0;
  double xn = 2.0;
  double h1 = 0.1;
  double h2 = 0.05;

  int n1 = (int)((xn - x0) / h1) + 1;
  int n2 = (int)((xn - x0) / h2) + 1;

  // Векторы для метода Эйлера
  vector<double> x1(n1), y1(n1), yp1(n1);
  for (int i = 0; i < n1; ++i) {
    x1[i] = x0 + i * h1;
  }
  y1[0] = 1.0;
  yp1[0] = 1.0;
  eulerMethod(x1, y1, yp1, h1);

  // Векторы для метода Рунге-Кутты
  vector<double> x2(n1), y2(n1), yp2(n1);
  for (int i = 0; i < n1; ++i) {
    x2[i] = x0 + i * h1;
  }
  y2[0] = 1.0;
  yp2[0] = 1.0;
  rungeKuttaMethod(x2, y2, yp2, h1);

  // Векторы для метода Адамса
  vector<double> x3(n1), y3(n1), yp3(n1);
  for (int i = 0; i < n1; ++i) {
    x3[i] = x0 + i * h1;
  }
  y3[0] = 1.0;
  yp3[0] = 1.0;
  adamsMethod(x3, y3, yp3, h1);

  // Векторы для вторых решений с шагом h2 (для метода Рунге - Ромберга)
  vector<double> x1_h2(n2), y1_h2(n2), yp1_h2(n2);
  for (int i = 0; i < n2; ++i) {
    x1_h2[i] = x0 + i * h2;
  }
  y1_h2[0] = 1.0;
  yp1_h2[0] = 1.0;
  eulerMethod(x1_h2, y1_h2, yp1_h2, h2);

  vector<double> x2_h2(n2), y2_h2(n2), yp2_h2(n2);
  for (int i = 0; i < n2; ++i) {
    x2_h2[i] = x0 + i * h2;
  }
  y2_h2[0] = 1.0;
  yp2_h2[0] = 1.0;
  rungeKuttaMethod(x2_h2, y2_h2, yp2_h2, h2);

  vector<double> x3_h2(n2), y3_h2(n2), yp3_h2(n2);
  for (int i = 0; i < n2; ++i) {
    x3_h2[i] = x0 + i * h2;
  }
  y3_h2[0] = 1.0;
  yp3_h2[0] = 1.0;
  adamsMethod(x3_h2, y3_h2, yp3_h2, h2);

  // Оценка погрешности по методу Рунге-Ромберга
  double R_euler = rungeRomberg(y1[n1 - 1], y1_h2[n2 - 1], 1);
  double R_rungeKutta = rungeRomberg(y2[n1 - 1], y2_h2[n2 - 1], 4);
  double R_adams = rungeRomberg(y3[n1 - 1], y3_h2[n2 - 1], 4);

  // Вывод результатов
  cout << "x       Euler Method    Runge-Kutta Method    Adams Method    Exact "
          "Solution"
       << endl;
  cout << fixed << setprecision(6);
  for (int i = 0; i < n1; ++i) {
    double exact = exactSolution(x1[i]);
    cout << x1[i] << "    " << y1[i] << "         " << y2[i] << "            "
         << y3[i] << "         " << exact << endl;
  }

  // Вывод оценок погрешности
  cout << "Погрешность метода Эйлера (h1 = 0.1): " << R_euler << endl;
  cout << "Погрешность метода Рунге-Кутты (h1 = 0.1): " << R_rungeKutta << endl;
  cout << "Погрешность метода Адамса (h1 = 0.1): " << R_adams << endl;

  return 0;
}