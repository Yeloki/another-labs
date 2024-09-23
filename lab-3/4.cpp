
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Функции для численного дифференцирования
double firstDerivativeCentered(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        if (x[i] == x_star) {
            return (y[i + 1] - y[i - 1]) / (2 * h);
        }
    }
    return 0.0;
}

double secondDerivativeCentered(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        if (x[i] == x_star) {
            return (y[i + 1] - 2 * y[i] + y[i - 1]) / (h * h);
        }
    }
    return 0.0;
}

double firstDerivativeForward(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 0; i < n - 1; ++i) {
        if (x[i] == x_star && i + 1 < n) {
            return (y[i + 1] - y[i]) / h;
        }
    }
    return 0.0;
}

double secondDerivativeForward(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 0; i < n - 2; ++i) {
        if (x[i] == x_star && i + 2 < n) {
            return (y[i + 2] - 2 * y[i + 1] + y[i]) / (h * h);
        }
    }
    return 0.0;
}

double firstDerivativeBackward(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 1; i < n; ++i) {
        if (x[i] == x_star) {
            return (y[i] - y[i - 1]) / h;
        }
    }
    return 0.0;
}

double secondDerivativeBackward(const vector<double>& x, const vector<double>& y, double x_star, double h) {
    int n = x.size();
    for (int i = 2; i < n; ++i) {
        if (x[i] == x_star) {
            return (y[i] - 2 * y[i - 1] + y[i - 2]) / (h * h);
        }
    }
    return 0.0;
}

int main() {
    vector<double> x = {0.0, 0.2, 0.4, 0.6, 0.8};
    vector<double> y = {1.0, 1.4214, 1.8918, 2.4221, 3.0255};
    double x_star = 0.4;
    double h = 0.2;

    // Вычисление производных
    double first_derivative_centered = firstDerivativeCentered(x, y, x_star, h);
    double second_derivative_centered = secondDerivativeCentered(x, y, x_star, h);

    double first_derivative_forward = firstDerivativeForward(x, y, x_star, h);
    double second_derivative_forward = secondDerivativeForward(x, y, x_star, h);

    double first_derivative_backward = firstDerivativeBackward(x, y, x_star, h);
    double second_derivative_backward = secondDerivativeBackward(x, y, x_star, h);

    // Вывод результатов
    cout << "Центрированные разности:" << endl;
    cout << "Первая производная: " << fixed << setprecision(6) << first_derivative_centered << endl;
    cout << "Вторая производная: " << fixed << setprecision(6) << second_derivative_centered << endl;

    cout << "\nНаправленные вправо разности:" << endl;
    cout << "Первая производная: " << fixed << setprecision(6) << first_derivative_forward << endl;
    cout << "Вторая производная: " << fixed << setprecision(6) << second_derivative_forward << endl;

    cout << "\nНаправленные влево разности:" << endl;
    cout << "Первая производная: " << fixed << setprecision(6) << first_derivative_backward << endl;
    cout << "Вторая производная: " << fixed << setprecision(6) << second_derivative_backward << endl;

    return 0;
}