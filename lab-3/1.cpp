#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Рассчитываем значение функции y = cot(x) + x
double objectiveFunction(double x) {
    return 1.0 / tan(x) + x;
}

// Функция для вычисления интерполяционного многочлена Лагранжа
double lagrangeInterpolation(const vector<double>& x, const vector<double>& y, double x_star) {
    double result = 0.0;
    int n = x.size();

    for (int i = 0; i < n; ++i) {
        double term = y[i];
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                term *= (x_star - x[j]) / (x[i] - x[j]);
            }
        }
        result += term;
    }

    return result;
}

// Функция для вычисления интерполяционного многочлена Ньютона
double newtonInterpolation(const vector<double>& x, const vector<double>& y, double x_star) {
    int n = x.size();
    vector<vector<double>> divided_diff(n, vector<double>(n));
    
    // заполняем 0-й столбец
    for (int i = 0; i < n; ++i) {
        divided_diff[i][0] = y[i];
    }
    
    // заполняем разделённые разности
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
    
    // вычисляем интерполяцию в точке x_star
    double result = divided_diff[0][0];
    double product = 1.0;
    for (int i = 1; i < n; ++i) {
        product *= (x_star - x[i - 1]);
        result += divided_diff[0][i] * product;
    }

    return result;
}

int main() {
    // Данные таблицы
    vector<double> x_a = {M_PI / 8, 2 * M_PI / 8, 3 * M_PI / 8, 4 * M_PI / 8};
    vector<double> x_b = {M_PI / 8, M_PI / 3, 3 * M_PI / 8, 3 * M_PI / 2};
    
    vector<double> y_a, y_b;
    for (double xi : x_a) {
        y_a.push_back(objectiveFunction(xi));
    }
    for (double xi : x_b) {
        y_b.push_back(objectiveFunction(xi));
    }

    // Точка, в которой будем оценивать интерполяцию
    double x_star = 3 * M_PI / 16;

    // Интерполяция многочленом Лагранжа
    double lagrange_result_a = lagrangeInterpolation(x_a, y_a, x_star);
    double lagrange_result_b = lagrangeInterpolation(x_b, y_b, x_star);

    cout << "Интерполяция многочленом Лагранжа (a): " << fixed << setprecision(6) << lagrange_result_a << endl;
    cout << "Интерполяция многочленом Лагранжа (б): " << fixed << setprecision(6) << lagrange_result_b << endl;
    
    // Интерполяция многочленом Ньютона
    double newton_result_a = newtonInterpolation(x_a, y_a, x_star);
    double newton_result_b = newtonInterpolation(x_b, y_b, x_star);
    
    cout << "Интерполяция многочленом Ньютона (a): " << fixed << setprecision(6) << newton_result_a << endl;
    cout << "Интерполяция многочленом Ньютона (б): " << fixed << setprecision(6) << newton_result_b << endl;

    // Точное значение функции в точке x_star
    double exact_value = objectiveFunction(x_star);
    cout << "Точное значение функции: " << fixed << setprecision(6) << exact_value << endl;

    // Вычисление погрешности интерполяции
    double lagrange_error_a = fabs(lagrange_result_a - exact_value);
    double lagrange_error_b = fabs(lagrange_result_b - exact_value);
    double newton_error_a = fabs(newton_result_a - exact_value);
    double newton_error_b = fabs(newton_result_b - exact_value);
    
    cout << "Погрешность интерполяции Лагранжа (a): " << fixed << setprecision(6) << lagrange_error_a << endl;
    cout << "Погрешность интерполяции Лагранжа (б): " << fixed << setprecision(6) << lagrange_error_b << endl;
    cout << "Погрешность интерполяции Ньютона (a): " << fixed << setprecision(6) << newton_error_a << endl;
    cout << "Погрешность интерполяции Ньютона (б): " << fixed << setprecision(6) << newton_error_b << endl;

    return 0;
}