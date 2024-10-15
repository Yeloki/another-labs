

import numpy as np
import math


def L2_norm(X):
    """
    Count ||X||_2
    """
    n = X.shape[0]
    l2_norm = 0
    for i in range(n):
        l2_norm += X[i] * X[i]
    return math.sqrt(l2_norm)


def solve_iterative(A, b, eps):
    """
    Uses iterative method to solve Ax=b
    Returns x and number of iterations
    """
    n = A.shape[0]

    # Step 1. Ax=b -> x = alpha * x + beta
    alpha = np.zeros_like(A, dtype="float")
    beta = np.zeros_like(b, dtype="float")
    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -A[i][j] / A[i][i]

        beta[i] = b[i] / A[i][i]

    # Step 2. Iterating
    iterations = 0
    cur_x = np.copy(beta)
    converge = False
    while not converge:
        prev_x = np.copy(cur_x)
        cur_x = alpha @ prev_x + beta
        iterations += 1
        converge = L2_norm(prev_x - cur_x) <= eps
    return cur_x, iterations


def seidel_multiplication(alpha, x, beta):
    """
    Count alhpa * x + beta for seidel method
    """
    res = np.copy(x)
    for i in range(alpha.shape[0]):
        res[i] = beta[i]
        for j in range(alpha.shape[1]):
            res[i] += alpha[i][j] * res[j]
    return res


def solve_seidel(A, b, eps):
    """
    Uses Seidel method to solve Ax=b
    Returns x and number of iterations
    """
    n = A.shape[0]

    # Step 1. Ax=b -> x = alpha * x + beta
    alpha = np.zeros_like(A, dtype="float")
    beta = np.zeros_like(b, dtype="float")
    for i in range(n):
        for j in range(n):
            if i == j:
                alpha[i][j] = 0
            else:
                alpha[i][j] = -A[i][j] / A[i][i]

        beta[i] = b[i] / A[i][i]

    # Step 2. Iterating
    iterations = 0
    cur_x = np.copy(beta)
    converge = False
    while not converge:
        prev_x = np.copy(cur_x)
        cur_x = seidel_multiplication(alpha, prev_x, beta)
        iterations += 1
        converge = L2_norm(prev_x - cur_x) <= eps
    return cur_x, iterations


def lab13():
    n = int(input("Enter the size of matrix: "))

    print("Enter matrix A")
    A = [[] for _ in range(n)]
    for i in range(n):
        row = list(map(int, input().split()))
        A[i] = row
    A = np.array(A, dtype="float")
    print("Enter vector b")
    b = np.array(list(map(int, input().split())), dtype="float")
    eps = float(input("Enter epsilon: "))

    print("Iteration method")
    x_iter, i_iter = solve_iterative(A, b, eps)
    print(x_iter)
    print("Iterations:", i_iter)
    print()

    print("Seidel method")
    x_seidel, i_seidel = solve_seidel(A, b, eps)
    print(x_seidel)
    print("Iterations:", i_seidel)


def read_and_run13(infile, outfile):
    f = open(infile, encoding='utf-8', mode='r')
    if outfile is not None:
        out = open(outfile, encoding='utf-8', mode='w')
    else:
        out = None
    data = list(map(lambda x: x.rstrip('\n').lstrip('\n'), (f.readlines())))
    
    n = int(data[0])

    A = [[] for _ in range(n)]
    for i in range(n):
        row = list(map(int, data[i + 1].split()))
        A[i] = row
    A = np.array(A, dtype="float")
    b = np.array(list(map(int, data[n + 1].split())), dtype="float")
    eps = float(data[n + 2])

    print("Iteration method", file=out)
    x_iter, i_iter = solve_iterative(A, b, eps)
    print(x_iter, file=out)
    print("Iterations:", i_iter, file=out)
    print(file=out)

    print("Seidel method", file=out)
    x_seidel, i_seidel = solve_seidel(A, b, eps)
    print(x_seidel, file=out)
    print("Iterations:", i_seidel, file=out)
    

if __name__ == "__main__":
    lab13()
