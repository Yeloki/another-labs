#include "common/matrix.hpp" // включаем заголовочный файл для работы с матрицами
#include <cassert>
#include <iostream>

int main() {
  // Демонстрация работы с классом Matrix

  // Создаем матрицу 2x2
  Matrix mat1(2, 2);
  mat1(0, 0) = 1;
  mat1(0, 1) = 2;
  mat1(1, 0) = 3;
  mat1(1, 1) = 4;

  std::cout << "Matrix 1:\n" << mat1 << std::endl;

  // Создаем еще одну матрицу 2x2
  Matrix mat2(2, 2);
  mat2(0, 0) = 5;
  mat2(0, 1) = 6;
  mat2(1, 0) = 7;
  mat2(1, 1) = 8;

  std::cout << "Matrix 2:\n" << mat2 << std::endl;

  // Сложение матриц
  Matrix sum = mat1 + mat2;
  std::cout << "Sum of Matrix 1 and Matrix 2:\n" << sum << std::endl;

  // Вычитание матриц
  Matrix diff = mat1 - mat2;
  std::cout << "Difference of Matrix 1 and Matrix 2:\n" << diff << std::endl;

  // Умножение матриц
  Matrix prod = mat1 * mat2;
  std::cout << "Product of Matrix 1 and Matrix 2:\n" << prod << std::endl;

  // Скалярное умножение
  Matrix scalarProd = mat1 * 2;
  std::cout << "Matrix 1 multiplied by 2:\n" << scalarProd << std::endl;

  // Транспонирование
  Matrix transposed = mat1.transpose();
  std::cout << "Transposed Matrix 1:\n" << transposed << std::endl;

  // Создание единичной матрицы
  Matrix identityMat = Matrix::identity(3);
  std::cout << "3x3 Identity Matrix:\n" << identityMat << std::endl;

  // Для демонстрации вычисления определителя и обратной матрицы,
  // убедитесь, что ваш класс Matrix реализует эти методы.
  // Здесь приведена заглушка, настоящую реализацию добавьте сами.
  try {
    double det = mat1.determinant(); // Вызов функции для вычисления
                                     // определителя (нужно реализовать)
    std::cout << "Determinant of Matrix 1: " << det << std::endl;

    Matrix inverseMat = mat1.inverse(); // Вызов функции для вычисления обратной
                                        // матрицы (нужно реализовать)
    std::cout << "Inverse of Matrix 1:\n" << inverseMat << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  // Создаем матрицу 3x3
  Matrix mat3(3, 3);
  mat3(0, 0) = 6;
  mat3(0, 1) = 1;
  mat3(0, 2) = 1;
  mat3(1, 0) = 4;
  mat3(1, 1) = -2;
  mat3(1, 2) = 5;
  mat3(2, 0) = 2;
  mat3(2, 1) = 8;
  mat3(2, 2) = 7;

  std::cout << "Matrix 3:\n" << mat3 << std::endl;

  // Вычисление определителя
  try {
    double det =
        mat3.determinant(); // Вызов функции для вычисления определителя
    std::cout << "Determinant of Matrix 3: " << det << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
  }

  Matrix mat4 = {
      {1, 2, 3}, //
      {4, 5, 6}, //
      {7, 8, 9}  //
  };

  std::cout << "Matrix 4:\n" << mat4 << std::endl;

  {
    Matrix mat = {
        {2, 3, 1, 5}, {6, 13, 5, 19}, {2, 19, 10, 23}, {4, 10, 11, 31}};

    std::cout << "Matrix:\n" << mat << std::endl;

    Matrix P, L, U;
    try {
      std::tie(P, L, U) =
          mat.luDecomposition(); // Вызов метода для LU-разложения
      std::cout << "P (Permutation Matrix):\n" << P << std::endl;
      std::cout << "L (Lower Triangular Matrix):\n" << L << std::endl;
      std::cout << "U (Upper Triangular Matrix):\n" << U << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }
  {

    // Пример матрицы А
    Matrix A = {{2, 3, 1, 5}, {6, 13, 5, 19}, {2, 19, 10, 23}, {4, 10, 11, 31}};

    // Вектор (матрица-столбец) B
    Matrix B = {{1}, {0}, {0}, {1}};

    // Решение системы AX = B
    Matrix X;
    try {
      X = A.solve(B);
      std::cout << "Solution X:\n" << X << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }
  {
    // Матрицы для трехдиагональной системы:
    // a - поддиагональная часть
    // b - диагональ
    // c - наддиагональная часть
    // d - правая часть уравнения
    // Матрица на входе:
    // 3 1 0
    // -1 4 2
    // 0 2 3
    //
    // Правая часть: 1 -3 5
    //
    // Ответ: 1 -2 3
    Matrix a({{0}, // a[0] не используется
              {-1},
              {2}});

    Matrix b({{3}, //
              {4},
              {3}});

    Matrix c({
        {1}, //
        {2},
        {0} // c[n-1] не используется
    });

    Matrix d({{1}, //
              {-3},
              {5}});

    Matrix sol({{1}, //
                {-2},
                {3}});

    // Решение системы AX = D методом прогонки
    Matrix x;

    try {
      x = Matrix::tridiagonalSolve(a, b, c, d); // Применяем метод прогонки

      std::cout << "Solution X:\n" << x << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }
  {
    // Пример матрицы A:
    Matrix A({{4, 1, 2}, //
              {3, 5, 1},
              {1, 1, 3}});

    // Вектор правой части b:
    Matrix b({{4}, //
              {7},
              {3}});

    // Решение системы AX = B методом Якоби
    Matrix x;
    try {
      x = A.jacobiSolve(b); // Применяем метод Якоби
      std::cout << "Solution Z:\n" << x << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }

  {
    // Пример матрицы A:
    Matrix A({{4, 1, 2}, {3, 5, 1}, {1, 1, 3}});

    // Вектор правой части b:
    Matrix b({{4}, {7}, {3}});

    // Решение системы AX = B методом Гаусса-Зейделя
    Matrix x;
    try {
      x = A.gaussSeidelSolve(b); // Применяем метод Гаусса-Зейделя
      std::cout << "Solution X:\n" << x << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }

  {
    // Пример симметричной матрицы A:
    Matrix A({{4, 1, 2}, {1, 3, 0}, {2, 0, 2}});

    // Поиск собственных значений и собственных векторов методом вращений
    Matrix eigenvalues, eigenvectors;
    try {
      std::tie(eigenvalues, eigenvectors) =
          A.jacobiRotation(); // Применяем метод вращений
      std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;
      std::cout << "Eigenvectors:\n" << eigenvectors << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }

  {
    // Пример симметричной матрицы A:
    Matrix A({
        {4, 1, 2}, //
        {1, 3, 0}, //
        {2, 0, 2}  //
    });

    // Поиск собственных значений методом QR-разложения
    Matrix eigenvalues;
    try {
      eigenvalues = A.QRAlgorithm(); // Применяем метод QR-разложения
      std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;
    } catch (const std::exception &e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
  }

  return 0;
}
