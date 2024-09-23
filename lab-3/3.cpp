
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Функция для решения системы линейных уравнений методом Гаусса
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        // Поиск главного элемента в колонке
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        // Перестановка строк
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        // Приведение к ступенчатому виду
        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Обратный ход метода Гаусса
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= A[k][i] * x[i];
        }
    }

    return x;
}

// Функция для решения нормальной системы уравнений для многочлена первой степени (линейная регрессия)
vector<double> fitLinear(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double sumX = 0, sumY = 0, sumXY = 0, sumXX = 0;

    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumXX += x[i] * x[i];
    }

    vector<vector<double>> A = {{static_cast<double>(n), sumX}, {sumX, sumXX}};
    vector<double> b = {sumY, sumXY};

    return gaussianElimination(A, b);
}

// Функция для решения нормальной системы уравнений для многочлена второй степени (квадратичная регрессия)
vector<double> fitQuadratic(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    double sumX = 0, sumX2 = 0, sumX3 = 0, sumX4 = 0;
    double sumY = 0, sumXY = 0, sumX2Y = 0;

    for (int i = 0; i < n; ++i) {
        sumX += x[i];
        double x2 = x[i] * x[i];
        double x3 = x2 * x[i];
        double x4 = x3 * x[i];
        sumX2 += x2;
        sumX3 += x3;
        sumX4 += x4;
        sumY += y[i];
        sumXY += x[i] * y[i];
        sumX2Y += x2 * y[i];
    }

    vector<vector<double>> A = {{static_cast<double>(n), sumX, sumX2}, {sumX, sumX2, sumX3}, {sumX2, sumX3, sumX4}};
    vector<double> b = {sumY, sumXY, sumX2Y};

    return gaussianElimination(A, b);
}

// Функция для вычисления суммы квадратов ошибок
double calculateSSE(const vector<double>& x, const vector<double>& y, const vector<double>& coeffs, int degree) {
    double sse = 0.0;

    for (int i = 0; i < x.size(); ++i) {
        double pred = 0.0;
        for (int j = 0; j <= degree; ++j) {
            pred += coeffs[j] * pow(x[i], j);
        }
        sse += pow(y[i] - pred, 2);
    }

    return sse;
}

int main() {
    vector<double> x = {1.0, 1.9, 2.8, 3.7, 4.6, 5.5};
    vector<double> y = {3.4142, 2.9818, 3.3095, 3.8184, 4.3599, 4.8318};

    // Нахождение коэффициентов для многочлена 1-ой и 2-ой степени
    vector<double> linearCoeff = fitLinear(x, y);
    vector<double> quadCoeff = fitQuadratic(x, y);

    cout << "Коэффициенты многочлена 1-ой степени: " << endl;
    cout << "a0 = " << linearCoeff[0] << ", a1 = " << linearCoeff[1] << endl;

    cout << "Коэффициенты многочлена 2-ой степени: " << endl;
    cout << "a0 = " << quadCoeff[0] << ", a1 = " << quadCoeff[1] << ", a2 = " << quadCoeff[2] << endl;

    // Расчет и вывод суммы квадратов ошибок
    double linearSSE = calculateSSE(x, y, linearCoeff, 1);
    double quadSSE = calculateSSE(x, y, quadCoeff, 2);
    cout << "Сумма квадратов ошибок для многочлена 1-ой степени: " << linearSSE << endl;
    cout << "Сумма квадратов ошибок для многочлена 2-ой степени: " << quadSSE << endl;

    return 0;
}