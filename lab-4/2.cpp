// CORRUPTED
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Точное решение
double exactSolution(double x) {
    return 1.0 + x + log(fabs(x));
}

// Функция для метода стрельбы
void rungeKutta(double x0, double y0, double yp0, double h, int n, vector<double>& x, vector<double>& y, vector<double>& yp) {
    x[0] = x0;
    y[0] = y0;
    yp[0] = yp0;
    for (int i = 1; i < n; ++i) {
        double k1y = yp[i - 1];
        double k1yp = (x[i - 1] * yp[i - 1] - y[i - 1]) / (x[i - 1] * x[i - 1] * log(fabs(x[i - 1])));

        double k2y = yp[i - 1] + 0.5 * h * k1yp;
        double k2yp = ((x[i - 1] + 0.5 * h) * (yp[i - 1] + 0.5 * h * k1yp) - (y[i - 1] + 0.5 * h * k1y)) / 
            ((x[i - 1] + 0.5 * h) * (x[i - 1] + 0.5 * h) * log(fabs(x[i - 1] + 0.5 * h)));

        double k3y = yp[i - 1] + 0.5 * h * k2yp;
        double k3yp = ((x[i - 1] + 0.5 * h) * (yp[i - 1] + 0.5 * h * k2yp) - (y[i - 1] + 0.5 * h * k2y)) / 
            ((x[i - 1] + 0.5 * h) * (x[i - 1] + 0.5 * h) * log(fabs(x[i - 1] + 0.5 * h)));

        double k4y = yp[i - 1] + h * k3yp;
        double k4yp = ((x[i - 1] + h) * (yp[i - 1] + h * k3yp) - (y[i - 1] + h * k3y)) / 
            ((x[i - 1] + h) * (x[i - 1] + h) * log(fabs(x[i - 1] + h)));

        y[i] = y[i - 1] + h * (k1y + 2 * k2y + 2 * k3y + k4y) / 6.0;
        yp[i] = yp[i - 1] + h * (k1yp + 2 * k2yp + 2 * k3yp + k4yp) / 6.0;
        x[i] = x[i - 1] + h;
    }
}

// Метод стрельбы
void shootingMethod(double x0, double xn, double y0, double ypn, double h, int n, vector<double>& x, vector<double>& y) {
    double yp0_1 = 0.5;  // Начальное приближение для y'(0)
    double yp0_2 = 1.0;  // Второе приближение для y'(0)

    vector<double> y1(n), y2(n), yp1(n), yp2(n);
    vector<double> x_vec(n);

    rungeKutta(x0, y0, yp0_1, h, n, x_vec, y1, yp1);
    rungeKutta(x0, y0, yp0_2, h, n, x_vec, y2, yp2);

    double shooting_error = 1e-6;
    double f1 = y1[n-1] - (ypn - y1[n-1]);
    double f2 = y2[n-1] - (ypn - y2[n-1]);
    
    while (fabs(f1 - f2) > shooting_error) {
        double yp0 = yp0_2 - f2 * (yp0_2 - yp0_1) / (f2 - f1);
        rungeKutta(x0, y0, yp0, h, n, x_vec, y, yp2);

        yp0_1 = yp0_2;
        yp0_2 = yp0;

        y1 = y2;
        y2 = y;

        f1 = f2;
        f2 = y[n-1] - (ypn - y[n-1]);
    }
}

// Конечно-разностный метод
void finiteDifferenceMethod(double x0, double xn, int n, vector<double>& x, vector<double>& y) {
    double h = (xn - x0) / (n - 1);

    vector<double> A(n), B(n), C(n), D(n);
    for (int i = 0; i < n; ++i) {
        x[i] = x0 + i * h;
    }

    for (int i = 1; i < n-1; ++i) {
        A[i] = (2 * x[i] * x[i] * log(fabs(x[i]))) - h * x[i];
        B[i] = -2 * (x[i] * x[i] * log(fabs(x[i])));
        C[i] = (2 * x[i] * x[i] * log(fabs(x[i]))) + h * x[i];
        D[i] = 0.0;
    }

    // Граничные условия
    A[0] = 0.0; B[0] = 1.0; C[0] = 0.0; D[0] = 0.0; // y'(0) = 0
    A[n-1] = 1.0; B[n-1] = -1.0; C[n-1] = 0.0; D[n-1] = 0.0; // y'(1) - y(1) = 0

    // Прогонка
    vector<double> alpha(n), beta(n);
    alpha[0] = -C[0] / B[0];
    beta[0] = D[0] / B[0];

    for (int i = 1; i < n; ++i) {
        alpha[i] = -C[i] / (A[i] * alpha[i-1] + B[i]);
        beta[i] = (D[i] - A[i] * beta[i-1]) / (A[i] * alpha[i-1] + B[i]);
    }

    y[n-1] = (D[n-1] - A[n-1] * beta[n-2]) / (A[n-1] * alpha[n-2] + B[n-1]);
    for (int i = n-2; i >= 0;--i) {
        y[i] = alpha[i] * y[i+1] + beta[i];
    }
}

// Метод Рунге-Ромберга для оценки погрешности
double rungeRomberg(double I_h1, double I_h2, int p) {
    return (I_h2 - I_h1) / (pow(2, p) - 1);
}

int main() {
    double x0 = 1.0;
    double xn = 2.0;
    int n1 = 11;  // Количество точек для h1 = 0.1
    int n2 = 21;  // Количество точек для h2 = 0.05

    // Векторы для метода стрельбы
    vector<double> x1(n1), y1(n1), x1_h2(n2), y1_h2(n2);

    // Векторы для конечно-разностного метода
    vector<double> x2(n1), y2(n1), x2_h2(n2), y2_h2(n2);

    // Решение методом стрельбы
    shootingMethod(x0, xn, 1.0, 0.0, 0.1, n1, x1, y1);
    shootingMethod(x0, xn, 1.0, 0.0, 0.05, n2, x1_h2, y1_h2);

    // Решение конечно-разностным методом
    finiteDifferenceMethod(x0, xn, n1, x2, y2);
    finiteDifferenceMethod(x0, xn, n2, x2_h2, y2_h2);

    // Оценка погрешности по методу Рунге-Ромберга
    double R_shooting = rungeRomberg(y1[n1-1], y1_h2[n2-1], 4);
    double R_fdm = rungeRomberg(y2[n1-1], y2_h2[n2-1], 2);

    // Вывод результатов
    cout << "x       Shooting Method    Finite Difference Method    Exact Solution" << endl;
    cout << fixed << setprecision(6);
    for (int i = 0; i < n1; ++i) {
        double exact = exactSolution(x1[i]);
        cout << x1[i] << "    " << y1[i] << "         " << y2[i] << "            " << exact << endl;
    }

    // Вывод оценок погрешности
    cout << "Погрешность метода стрельбы (h1 = 0.1): " << R_shooting << endl;
    cout << "Погрешность конечно-разностного метода (h1 = 0.1): " << R_fdm << endl;

    return 0;
}