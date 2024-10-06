#include "common/matrix.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <iomanip>
#include <stdexcept>

Matrix::Matrix() : rows(0), cols(0), data(nullptr) {}

Matrix::Matrix(int rows, int cols)
    : rows(rows), cols(cols), data(new double[rows * cols]()) {}

Matrix::Matrix(const Matrix &other) { copyFrom(other); }

Matrix::Matrix(Matrix &&other) noexcept
    : rows(other.rows), cols(other.cols), data(other.data) {
  other.data = nullptr;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> init) {
  rows = init.size();
  if (rows > 0) {
    cols = init.begin()->size();
  } else {
    cols = 0;
  }
  data = new double[rows * cols];

  int i = 0;
  for (const auto &row : init) {
    if (row.size() != cols) {
      delete[] data;
      throw std::invalid_argument(
          "All rows must have the same number of columns");
    }
    std::copy(row.begin(), row.end(), data + i * cols);
    ++i;
  }
}

Matrix &Matrix::operator=(const Matrix &other) {
  if (this != &other) {
    delete[] data;
    copyFrom(other);
  }
  return *this;
}

Matrix &Matrix::operator=(Matrix &&other) noexcept {
  if (this != &other) {
    delete[] data;
    rows = other.rows;
    cols = other.cols;
    data = other.data;
    other.data = nullptr;
  }
  return *this;
}

Matrix::~Matrix() { delete[] data; }

double &Matrix::operator()(int row, int col) {
  if (row >= rows || col >= cols || row < 0 || col < 0) {
    throw std::out_of_range("Matrix indices out of range");
  }
  return data[row * cols + col];
}

const double &Matrix::operator()(int row, int col) const {
  if (row >= rows || col >= cols || row < 0 || col < 0) {
    throw std::out_of_range("Matrix indices out of range");
  }
  return data[row * cols + col];
}

Matrix Matrix::operator+(const Matrix &rhs) const {
  if (rows != rhs.rows || cols != rhs.cols) {
    throw std::invalid_argument("Matrix dimensions must match for addition");
  }
  Matrix result(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    result.data[i] = data[i] + rhs.data[i];
  }
  return result;
}

Matrix Matrix::operator-(const Matrix &rhs) const {
  if (rows != rhs.rows || cols != rhs.cols) {
    throw std::invalid_argument("Matrix dimensions must match for subtraction");
  }
  Matrix result(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    result.data[i] = data[i] - rhs.data[i];
  }
  return result;
}

Matrix Matrix::operator*(const Matrix &rhs) const {
  if (cols != rhs.rows) {
    throw std::invalid_argument(
        "Matrix dimensions must match for multiplication");
  }
  Matrix result(rows, rhs.cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < rhs.cols; ++j) {
      for (int k = 0; k < cols; ++k) {
        result(i, j) += (*this)(i, k) * rhs(k, j);
      }
    }
  }
  return result;
}

Matrix Matrix::operator*(double scalar) const {
  Matrix result(rows, cols);
  for (int i = 0; i < rows * cols; ++i) {
    result.data[i] = data[i] * scalar;
  }
  return result;
}

Matrix Matrix::transpose() const {
  Matrix result(cols, rows);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(j, i) = (*this)(i, j);
    }
  }
  return result;
}

double Matrix::determinant() const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square to compute determinant");
  }
  return determinantRecursive(*this);
}

double Matrix::determinantRecursive(const Matrix &mat) const {
  int n = mat.getRows();
  if (n == 1) {
    return mat(0, 0);
  }

  if (n == 2) {
    return (mat(0, 0) * mat(1, 1)) - (mat(0, 1) * mat(1, 0));
  }

  double det = 0;
  for (int p = 0; p < n; ++p) {
    Matrix subMatrix(n - 1, n - 1);
    for (int i = 1; i < n; ++i) {
      int col = 0;
      for (int j = 0; j < n; ++j) {
        if (j == p)
          continue;
        subMatrix(i - 1, col++) = mat(i, j);
      }
    }
    det += mat(0, p) * pow(-1, p) * determinantRecursive(subMatrix);
  }
  return det;
}

Matrix Matrix::inverse() const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square to compute inverse");
  }

  int n = rows;
  Matrix augmentedMatrix(n, 2 * n);

  // Создаем расширенную матрицу [A|I]
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      augmentedMatrix(i, j) = (*this)(i, j);
    }
    augmentedMatrix(i, n + i) = 1.0; // Заполняем единичную матрицу
  }

  // Прямой ход метода Гаусса для приведения к верхнетреугольному виду
  for (int i = 0; i < n; ++i) {
    if (augmentedMatrix(i, i) == 0) {
      // Если элемент на диагонали равен нулю, меняем строки местами
      bool swapped = false;
      for (int k = i + 1; k < n; ++k) {
        if (augmentedMatrix(k, i) != 0) {
          for (int j = 0; j < 2 * n; ++j) {
            std::swap(augmentedMatrix(i, j), augmentedMatrix(k, j));
          }
          swapped = true;
          break;
        }
      }
      if (!swapped) {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
      }
    }

    // Нормализуем текущую строку так, чтобы элемент на диагонали стал равен 1
    double diagValue = augmentedMatrix(i, i);
    for (int j = 0; j < 2 * n; ++j) {
      augmentedMatrix(i, j) /= diagValue;
    }

    // Вычитаем текущую строку из всех остальных, чтобы получить нули в текущем
    // столбце
    for (int k = 0; k < n; ++k) {
      if (k == i)
        continue;
      double factor = augmentedMatrix(k, i);
      for (int j = 0; j < 2 * n; ++j) {
        augmentedMatrix(k, j) -= factor * augmentedMatrix(i, j);
      }
    }
  }

  // Извлекаем обратную матрицу из расширенной матрицы
  Matrix inverseMatrix(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      inverseMatrix(i, j) = augmentedMatrix(i, n + j);
    }
  }

  return inverseMatrix;
}

void Matrix::print(std::ostream &os) const {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      os << (*this)(i, j) << " ";
    }
    os << std::endl;
  }
}

Matrix Matrix::identity(int size) {
  Matrix result(size, size);
  for (int i = 0; i < size; ++i) {
    result(i, i) = 1.0;
  }
  return result;
}

void Matrix::copyFrom(const Matrix &other) {
  rows = other.rows;
  cols = other.cols;
  data = new double[rows * cols];
  std::copy(other.data, other.data + (rows * cols), data);
}

bool Matrix::operator==(const Matrix &other) const {
  if (rows != other.rows || cols != other.cols) {
    return false;
  }
  for (int i = 0; i < rows * cols; ++i) {
    if (data[i] != other.data[i]) {
      return false;
    }
  }
  return true;
}

bool Matrix::operator!=(const Matrix &other) const { return !(*this == other); }

std::ostream &operator<<(std::ostream &os, const Matrix &mat) {
  int width = 7;     // Ширина столбца для вывода
  int precision = 4; // Точность вывода

  os << std::fixed
     << std::setprecision(precision); // Фиксированная точка и заданная точность
  for (int i = 0; i < mat.getRows(); ++i) {
    for (int j = 0; j < mat.getCols(); ++j) {
      os << std::setw(width) << mat(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}

std::tuple<Matrix, Matrix, Matrix> Matrix::luDecomposition() const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square for LU decomposition");
  }

  int n = rows;
  Matrix L = Matrix::identity(n); // Нижняя треугольная матрица
  Matrix U = *this; // Складываем в U текущую матрицу
  Matrix P = Matrix::identity(
      n); // Матрица перестановок (матрица выбора главного элемента)

  for (int i = 0; i < n; ++i) {
    // Выбираем главный элемент и выполним перестановки строк
    double maxElem = std::abs(U(i, i));
    int pivotRow = i;
    for (int k = i + 1; k < n; ++k) {
      if (std::abs(U(k, i)) > maxElem) {
        maxElem = std::abs(U(k, i));
        pivotRow = k;
      }
    }

    if (pivotRow != i) {
      // Перестановка строк в U
      for (int j = 0; j < n; ++j) {
        std::swap(U(i, j), U(pivotRow, j));
      }
      // Перестановка строк в P
      for (int j = 0; j < n; ++j) {
        std::swap(P(i, j), P(pivotRow, j));
      }
      // Перестановка строк в L
      if (i > 0) {
        for (int j = 0; j < i; ++j) {
          std::swap(L(i, j), L(pivotRow, j));
        }
      }
    }

    // Разложение на L и U
    for (int k = i + 1; k < n; ++k) {
      L(k, i) = U(k, i) / U(i, i);
      for (int j = i; j < n; ++j) {
        U(k, j) -= L(k, i) * U(i, j);
      }
    }
  }

  return std::make_tuple(P, L, U);
}
Matrix Matrix::solve(const Matrix &B) const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix A must be square for solving AX = B");
  }

  if (B.getRows() != rows) {
    throw std::invalid_argument(
        "Matrix B must have the same number of rows as matrix A");
  }

  // Получаем LU-разложение
  Matrix P, L, U;
  std::tie(P, L, U) = luDecomposition();

  // Применяем перестановку к B: PB
  Matrix PB = P * B;

  // Решаем LY = PB прямым ходом
  Matrix Y(rows, 1);
  for (int i = 0; i < rows; ++i) {
    Y(i, 0) = PB(i, 0);
    for (int j = 0; j < i; ++j) {
      Y(i, 0) -= L(i, j) * Y(j, 0);
    }
    Y(i, 0) /= L(i, i);
  }

  // Решаем UX = Y обратным ходом
  Matrix X(rows, 1);
  for (int i = rows - 1; i >= 0; --i) {
    X(i, 0) = Y(i, 0);
    for (int j = i + 1; j < cols; ++j) {
      X(i, 0) -= U(i, j) * X(j, 0);
    }
    X(i, 0) /= U(i, i);
  }

  return X;
}

Matrix Matrix::tridiagonalSolve(const Matrix &a, const Matrix &b,
                                const Matrix &c, const Matrix &d) {
  if (b.getCols() != 1 || c.getCols() != 1 || d.getCols() != 1 ||
      a.getCols() != 1) {
    throw std::invalid_argument("All input matrices must be column vectors.");
  }

  int n = b.getRows();
  Matrix x(n, 1); // Решение

  // Прогоночные коэффициенты
  Matrix p(n, 1); // Решение
  Matrix q(n, 1); // Решение

  // Инициализация первых прогонов
  p(0, 0) = -c(0, 0) / b(0, 0);
  q(0, 0) = d(0, 0) / b(0, 0);

  // Прямой ход прогонки
  for (int i = 1; i < n; ++i) {
    double denom = b(i, 0) + a(i, 0) * p(i - 1, 0);
    p(i, 0) = -c(i, 0) / denom;
    q(i, 0) = (d(i, 0) - a(i, 0) * q(i - 1, 0)) / denom;
  }

  // Обратный ход прогонки
  x(n - 1, 0) = q(n - 1, 0);
  for (int i = n - 2; i >= 0; --i) {
    x(i, 0) = p(i, 0) * x(i + 1, 0) + q(i, 0);
  }

  return x;
}

Matrix Matrix::jacobiSolve(const Matrix &b, double tolerance,
                           int maxIterations) const {
  if (b.getRows() != rows || b.getCols() != 1) {
    throw std::invalid_argument(
        "Right-hand side vector has incorrect dimensions.");
  }

  Matrix x(rows, 1); // Начальное приближение
  Matrix x_new(rows, 1);

  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    for (int i = 0; i < rows; ++i) {
      double sum = 0;
      for (int j = 0; j < cols; ++j) {
        if (i != j) {
          sum += (*this)(i, j) * x(j, 0);
        }
      }
      x_new(i, 0) = (b(i, 0) - sum) / (*this)(i, i);
    }

    // Проверка условия остановки
    double maxError = 0;
    for (int i = 0; i < rows; ++i) {
      maxError = std::max(maxError, std::abs(x_new(i, 0) - x(i, 0)));
    }

    if (maxError < tolerance) {
      return x_new;
    }

    x = x_new;
  }

  throw std::runtime_error("Jacobi method did not converge within the maximum "
                           "number of iterations.");
}

Matrix Matrix::gaussSeidelSolve(const Matrix &b, double tolerance,
                                int maxIterations) const {
  if (b.getRows() != rows || b.getCols() != 1) {
    throw std::invalid_argument(
        "Right-hand side vector has incorrect dimensions.");
  }

  Matrix x(rows, 1); // Начальное приближение

  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    Matrix x_old = x;

    for (int i = 0; i < rows; ++i) {
      double sum = 0;
      for (int j = 0; j < cols; ++j) {
        if (i != j) {
          sum += (*this)(i, j) * x(j, 0);
        }
      }
      x(i, 0) = (b(i, 0) - sum) / (*this)(i, i);
    }

    // Проверка условия остановки
    double maxError = 0;
    for (int i = 0; i < rows; ++i) {
      maxError = std::max(maxError, std::abs(x(i, 0) - x_old(i, 0)));
    }

    if (maxError < tolerance) {
      return x;
    }
  }

  throw std::runtime_error("Gauss-Seidel method did not converge within the "
                           "maximum number of iterations.");
}

std::tuple<Matrix, Matrix> Matrix::jacobiRotation(double tolerance,
                                                  int maxIterations) const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square.");
  }

  int n = rows;
  Matrix A(*this); // Копируем текущую матрицу
  Matrix V = Matrix::identity(
      n); // Единичная матрица для накопления собственных векторов

  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    // Находим максимальный внедиагональный элемент
    double maxOffDiagonal = 0;
    int p = 0, q = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i != j && std::abs(A(i, j)) > maxOffDiagonal) {
          maxOffDiagonal = std::abs(A(i, j));
          p = i;
          q = j;
        }
      }
    }

    // Проверяем условие остановки
    if (maxOffDiagonal < tolerance) {
      break;
    }

    // Вычисляем параметры вращения
    double phi = 0.5 * std::atan2(2 * A(p, q), A(p, p) - A(q, q));
    double c = std::cos(phi);
    double s = std::sin(phi);

    Matrix rot = Matrix::identity(n); // Получаем матрицу вращения

    rot(p, p) = c;
    rot(p, q) = -s;
    rot(q, p) = s;
    rot(q, q) = c;

    V = V * rot;

    A = rot.transpose() * A * rot;
  }

  // Собственные значения находились на диагонали A, собственные векторы в V
  Matrix eigenvalues(n, 1);
  for (int i = 0; i < n; ++i) {
    eigenvalues(i, 0) = A(i, i);
  }

  return std::make_tuple(eigenvalues, V);
}

Matrix Matrix::getRow(size_t row_idx) const {
  Matrix result(1, cols);
  for (int i(0); i < cols; ++i) {
    result(0, i) = this->operator()(row_idx, i);
  }
  return result;
}

Matrix Matrix::getCol(size_t col_idx) const {
  Matrix result(rows, 1);
  for (int i(0); i < rows; ++i) {
    result(i, 0) = this->operator()(i, col_idx);
  }
  return result;
}

std::tuple<Matrix, Matrix> Matrix::QRDecomposition() const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square.");
  }

  const auto sign = [](auto v) -> int {
    if (v > 0)
      return 1;
    if (v < 0)
      return -1;
    return 0;
  };

  int n = rows;
  Matrix Q = Matrix::identity(n);
  Matrix R(*this);

  for (int k = 0; k < n - 1; ++k) {
    Matrix H = Matrix::identity(n);

    double a = R(k, k);
    double b = R(k + 1, k);
    double r = std::sqrt(a * a + b * b);
    double c = a / r;
    double s = -b / r;

    H(k, k) = c;
    H(k, k + 1) = -s;
    H(k + 1, k) = s;
    H(k + 1, k + 1) = c;

    R = H * R;
    Q = Q * H.transpose();
  }

  return std::make_tuple(Q, R);
}

Matrix Matrix::QRAlgorithm(double tolerance, int maxIterations) const {
  if (rows != cols) {
    throw std::invalid_argument("Matrix must be square.");
  }

  Matrix A(*this);

  for (int iteration = 0; iteration < maxIterations; ++iteration) {
    Matrix Q, R;
    std::tie(Q, R) = A.QRDecomposition();

    A = R * Q;

    // Проверка на сходимость
    bool converged = true;
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < i; ++j) {
        if (std::abs(A(i, j)) > tolerance) {
          converged = false;
          break;
        }
      }
    }

    if (converged) {
      break;
    }
  }

  Matrix eigenvalues(rows, 1);
  for (int i = 0; i < rows; ++i) {
    eigenvalues(i, 0) = A(i, i);
  }

  return eigenvalues;
}