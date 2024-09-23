#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Функция для построения и вычисления кубического сплайна
vector<vector<double>> cubicSpline(const vector<double>& x, const vector<double>& y) {
    int n = x.size() - 1;
    
    // Шаги между узлами
    vector<double> h(n);
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    
    // Альфа
    vector<double> alpha(n);
    for (int i = 1; i < n; ++i) {
        alpha[i] = (3 / h[i]) * (y[i + 1] - y[i]) - (3 / h[i - 1]) * (y[i] - y[i - 1]);
    }

    // Решение трехдиагональной системы методом прогонки
    vector<double> l(n + 1), mu(n + 1), z(n + 1);
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (int i = 1; i < n; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1.0;
    z[n] = 0.0;

    // Вычисление коэффициентов
    vector<double> a(n + 1), b(n), c(n + 1), d(n);
    for (int i = 0; i <= n; ++i) {
        a[i] = y[i];
    }
    c[n] = 0.0;
    for (int j = n - 1; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    vector<vector<double>> coefficients(n, vector<double>(4));
    for (int i = 0; i < n; ++i) {
        coefficients[i][0] = a[i];
        coefficients[i][1] = b[i];
        coefficients[i][2] = c[i];
        coefficients[i][3] = d[i];
    }

    return coefficients;
}

double evaluateSpline(const vector<vector<double>>& coefficients, const vector<double>& x, double x_star) {
    int n = x.size() - 1;
    int interval = 0;
    
    for (int i = 0; i < n; ++i) {
        if (x_star >= x[i] && x_star <= x[i + 1]) {
            interval = i;
            break;
        }
    }

    double dx = x_star - x[interval];
    return coefficients[interval][0] +
           coefficients[interval][1] * dx +
           coefficients[interval][2] * dx * dx +
           coefficients[interval][3] * dx * dx * dx;
}

int main() {
    vector<double> x = {1.0, 1.9, 2.8, 3.7, 4.6};
    vector<double> y = {2.8069, 1.8279, 1.6091, 1.5713, 1.5663};
    double x_star = 2.66666667;

    // Построение кубического сплайна
    vector<vector<double>> coefficients = cubicSpline(x, y);

    // Оценка значения сплайна в точке x_star
    double spline_result = evaluateSpline(coefficients, x, x_star);
    cout << "Интерполяция кубическим сплайном: " << fixed << setprecision(6) << spline_result << endl;

    return 0;
}
