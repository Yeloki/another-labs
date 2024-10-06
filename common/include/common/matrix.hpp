#pragma once
#include <initializer_list>
#include <ostream>
#include <vector>

class Matrix {
public:
  Matrix();
  Matrix(int rows, int cols);
  Matrix(const Matrix &other); // Копирующий конструктор
  Matrix(Matrix &&other) noexcept; // Перемещающий конструктор

  Matrix(std::initializer_list<std::initializer_list<double>> init);

  Matrix &operator=(const Matrix &other); // Оператор копирующего присваивания
  Matrix &
  operator=(Matrix &&other) noexcept; // Оператор перемещающего присваивания
  ~Matrix();

  int getRows() const { return rows; }
  int getCols() const { return cols; }
  Matrix getRow(size_t row_idx) const;
  Matrix getCol(size_t col_idx) const;

  double &operator()(int row, int col); // Оператор доступа к элементу
  const double &operator()(int row, int col) const;

  Matrix operator+(const Matrix &rhs) const; // Сложение матриц
  Matrix operator-(const Matrix &rhs) const; // Вычитание матриц
  Matrix operator*(const Matrix &rhs) const; // Умножение матриц
  Matrix operator*(double scalar) const; // Скалярное умножение

  Matrix transpose() const; // Транспонирование

  double determinant() const; // Определитель

  Matrix inverse() const; // Обратная матрица
  std::tuple<Matrix, Matrix, Matrix> luDecomposition() const;
  void print(std::ostream &os) const; // Вывод

  static Matrix identity(int size); // Единичная матрица

  Matrix solve(const Matrix &B) const;
  static Matrix tridiagonalSolve(const Matrix &a, const Matrix &b,
                                 const Matrix &c, const Matrix &d);

  bool operator==(const Matrix &other) const;
  bool operator!=(const Matrix &other) const;

  Matrix jacobiSolve(const Matrix &b, double tolerance = 1e-10,
                     int maxIterations = 1000) const;

  Matrix gaussSeidelSolve(const Matrix &b, double tolerance = 1e-10,
                          int maxIterations = 1000) const;

  std::tuple<Matrix, Matrix> jacobiRotation(double tolerance = 1e-10,
                                            int maxIterations = 1000) const;

  std::tuple<Matrix, Matrix> QRDecomposition() const;
  Matrix QRAlgorithm(double tolerance = 1e-10, int maxIterations = 1000) const;

private:
  int rows, cols;
  double *data;

  void copyFrom(const Matrix &other);
  double determinantRecursive(const Matrix &mat) const;
};

std::ostream &operator<<(std::ostream &os, const Matrix &mat);
