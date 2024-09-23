#include <iostream>
#include <cmath>

using namespace std;

// Функция для нелинейного уравнения
double f(double x) {
    return sin(x) - x * x + 1;
}

// Производная функции для метода Ньютона
double df(double x) {
    return cos(x) - 2 * x;
}

// Метод простой итерации: x = g(x)
double g(double x) {
    return sqrt(sin(x) + 1);
}

// Метод простой итерации
double simpleIterationMethod(double (*g)(double), double initial_guess, double tol, int& iterations) {
    double x_prev = initial_guess;
    double x_next = g(x_prev);
    iterations = 1;

    while (fabs(x_next - x_prev) > tol) {
        x_prev = x_next;
        x_next = g(x_prev);
        iterations++;
    }

    return x_next;
}

// Метод Ньютона
double newtonMethod(double (*f)(double), double (*df)(double), double initial_guess, double tol, int& iterations) {
    double x_prev = initial_guess;
    double x_next = x_prev - f(x_prev) / df(x_prev);
    iterations = 1;

    while (fabs(x_next - x_prev) > tol) {
        x_prev = x_next;
        x_next = x_prev - f(x_prev) / df(x_prev);
        iterations++;
    }

    return x_next;
}

int main() {
    double initial_guess = 1.0; // Запупчим с понятного начального приближения
    double tol = 1e-6; // Заданная точность
    int iterations;

    // Метод простой итерации - решение
    double root_simple_iteration = simpleIterationMethod(g, initial_guess, tol, iterations);
    cout << "Метод простой итерации:" << endl;
    cout << "Корень: " << root_simple_iteration << endl;
    cout << "Число итераций: " << iterations << endl;

    // Метод Ньютона - решение
    double root_newton = newtonMethod(f, df, initial_guess, tol, iterations);
    cout << "Метод Ньютона:" << endl;
    cout << "Корень: " << root_newton << endl;
    cout << "Число итераций: " << iterations << endl;

    return 0;
}
