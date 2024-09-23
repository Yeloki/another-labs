#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Функция для интегрирования
double f(double x) {
    return x / (pow(x, 4) + 81);
}

// Функция для вычисления значения интеграла методом левых прямоугольников
double leftRectangles(double (*func)(double), double a, double b, double h) {
    int n = (b - a) / h;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        sum += func(a + i * h);
    }
    return h * sum;
}

// Функция для вычисления значения интеграла методом трапеций
double trapezoidal(double (*func)(double), double a, double b, double h) {
    int n = (b - a) / h;
    double sum = 0.0;
    for (int i = 1; i < n; ++i) {
        sum += func(a + i * h);
    }
    return h * (func(a) + 2 * sum + func(b)) / 2.0;
}

// Функция для вычисления значения интеграла методом Симпсона
double simpson(double (*func)(double), double a, double b, double h) {
    int n = (b - a) / h;
    double sum_odd = 0.0;
    double sum_even = 0.0;
    for (int i = 1; i < n; i += 2) {
        sum_odd += func(a + i * h);
    }
    for (int i = 2; i < n; i += 2) {
        sum_even += func(a + i * h);
    }
    return h * (func(a) + 4 * sum_odd + 2 * sum_even + func(b)) / 3.0;
}

// Метод Рунге-Ромберга
double rungeRomberg(double I_h1, double I_h2, int p) {
    return (I_h2 - I_h1) / (pow(2, p) - 1);
}

int main() {
    double a = 0.0;
    double b = 2.0;
    double h1 = 0.5;
    double h2 = 0.25;

    // Вычисление интегралов по методам левых прямоугольников, трапеций и Симпсона при шаге h1
    double I_left_h1 = leftRectangles(f, a, b, h1);
    double I_trap_h1 = trapezoidal(f, a, b, h1);
    double I_simpson_h1 = simpson(f, a, b, h1);

    // Вычисление интегралов по методам левых прямоугольников, трапеций и Симпсона при шаге h2
    double I_left_h2 = leftRectangles(f, a, b, h2);
    double I_trap_h2 = trapezoidal(f, a, b, h2);
    double I_simpson_h2 = simpson(f, a, b, h2);

    // Оценка погрешности по методу Рунге-Ромберга
    double R_left = rungeRomberg(I_left_h1, I_left_h2, 1);
    double R_trap = rungeRomberg(I_trap_h1, I_trap_h2, 2);
    double R_simpson = rungeRomberg(I_simpson_h1, I_simpson_h2, 4);

    // Вывод результатов
    cout << "Интеграл методом левых прямоугольников: " << fixed << setprecision(6) << I_left_h2 << endl;
    cout << "Погрешность по методу Рунге-Ромберга: " << fixed << setprecision(6) << R_left << endl;

    cout << "Интеграл методом трапеций: " << fixed << setprecision(6) << I_trap_h2 << endl;
    cout << "Погрешность по методу Рунге-Ромберга: " << fixed << setprecision(6) << R_trap << endl;

    cout << "Интеграл методом Симпсона: " << fixed << setprecision(6) << I_simpson_h2 << endl;
    cout << "Погрешность по методу Рунге-Ромберга: " << fixed << setprecision(6) << R_simpson << endl;

    return 0;
}
