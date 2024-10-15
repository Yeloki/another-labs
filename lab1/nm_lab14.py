

import numpy as np


def find_max_upper_element(X):
    """
    Find coords of max element by absolute value above the main diagonal
    Returns i, j of max element
    """
    n = X.shape[0]
    i_max, j_max = 0, 1
    max_elem = abs(X[0][1])
    for i in range(n):
        for j in range(i + 1, n):
            if abs(X[i][j]) > max_elem:
                max_elem = abs(X[i][j])
                i_max = i
                j_max = j
    return i_max, j_max


def matrix_norm(X):
    """
    Calculates L2 norm for elements above the main diagonal
    """
    norm = 0
    n = len(X)
    for i in range(n):
        for j in range(i + 1, n):
            norm += X[i][j] * X[i][j]
    return np.sqrt(norm)


def rotation_method(A, eps):
    """
    Find eigen values and eigen vectors using rotation method
    Returns eigen values, eigen vectors, number of iterations
    """
    n = A.shape[0]
    A_i = np.copy(A)
    eigen_vectors = np.eye(n)
    iterations = 0

    while matrix_norm(A_i) > eps:
        i_max, j_max = find_max_upper_element(A_i)
        if A_i[i_max][i_max] - A_i[j_max][j_max] == 0:
            phi = np.pi / 4
        else:
            phi = 0.5 * np.arctan(
                2 * A_i[i_max][j_max] / (A_i[i_max][i_max] - A_i[j_max][j_max])
            )

        # create rotation matrix
        U = np.eye(n)
        U[i_max][j_max] = -np.sin(phi)
        U[j_max][i_max] = np.sin(phi)
        U[i_max][i_max] = np.cos(phi)
        U[j_max][j_max] = np.cos(phi)

        A_i = U.T @ A_i @ U
        eigen_vectors = eigen_vectors @ U
        iterations += 1

    eigen_values = np.array([A_i[i][i] for i in range(n)])
    return eigen_values, eigen_vectors, iterations


def lab14():
    n = int(input("Enter the size of matrix: "))

    print("Enter matrix A")
    A = [[] for _ in range(n)]
    for i in range(n):
        row = list(map(int, input().split()))
        A[i] = row
    A = np.array(A, dtype="float")
    eps = float(input("Enter epsilon: "))

    eig_values, eig_vectors, iters = rotation_method(A, eps)
    print("Eigen values:", eig_values)
    print("Eigen vectors")
    print(eig_vectors)
    print("Iterations:", iters)


def read_and_run14(infile, outfile):
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
    eps = float(data[n + 1])

    eig_values, eig_vectors, iters = rotation_method(A, eps)
    print("Eigen values:", eig_values, file=out)
    print("Eigen vectors",file=out)
    print(eig_vectors,file=out)
    print("Iterations:", iters,file=out)

if __name__ == "__main__":
    lab14()
