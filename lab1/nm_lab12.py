
def read_tridiagonal_matrix(n, data=None):
    """
    Get tri-diagonal matrix with n rows by reading only not-null elements
    """
    if data is None:
        A = [[0 for _ in range(n)] for _ in range(n)]
        A[0][0], A[0][1] = map(int, input().split())
        for i in range(1, n - 1):
            A[i][i - 1], A[i][i], A[i][i + 1] = map(int, input().split())
        A[n - 1][n - 2], A[n - 1][n - 1] = map(int, input().split())
        return A
    else:
        A = [[0 for _ in range(n)] for _ in range(n)]
        A[0][0], A[0][1] = map(int, data[1].split())
        for i in range(1, n - 1):
            A[i][i - 1], A[i][i], A[i][i + 1] = map(int, data[i + 1].split())
        A[n - 1][n - 2], A[n - 1][n - 1] = map(int, data[-2].split())
        return A


def tridiagonal_solve(A, b):
    """
    Solves Ax=b, where A - tri-diagonal matrix
    Returns x
    """
    n = len(A)
    # Step 1. Forward
    v = [0 for _ in range(n)]
    u = [0 for _ in range(n)]
    v[0] = A[0][1] / -A[0][0]
    u[0] = b[0] / A[0][0]
    for i in range(1, n - 1):
        v[i] = A[i][i + 1] / (-A[i][i] - A[i][i - 1] * v[i - 1])
        u[i] = (A[i][i - 1] * u[i - 1] - b[i]) / (-A[i][i] - A[i][i - 1] * v[i - 1])
    v[n - 1] = 0
    u[n - 1] = (A[n - 1][n - 2] * u[n - 2] - b[n - 1]) / (
        -A[n - 1][n - 1] - A[n - 1][n - 2] * v[n - 2]
    )

    # Step 2. Backward
    x = [0 for _ in range(n)]
    x[n - 1] = u[n - 1]
    for i in range(n - 1, 0, -1):
        x[i - 1] = v[i - 1] * x[i] + u[i - 1]
    return x


def read_and_run12(infile, outfile):
    f = open(infile, encoding='utf-8', mode='r')
    if outfile is not None:
        out = open(outfile, encoding='utf-8', mode='w')
    else:
        out = None
    data = list(map(lambda x: x.rstrip('\n').lstrip('\n'), (f.readlines())))
    n = int(data[0])
    A = read_tridiagonal_matrix(n, data)
    b = list(map(int, data[-1].split()))

    print("Solution", file=out)
    x = tridiagonal_solve(A, b)
    print(x, file=out)

def lab12():
    n = int(input("Enter the number of equations: "))
    print("Enter not null elements of tri-diagonal matrix")
    A = read_tridiagonal_matrix(n)
    print("Enter vector b")
    b = list(map(int, input().split()))

    print("Solution")
    x = tridiagonal_solve(A, b)
    print(x)


if __name__ == "__main__":
    lab12()
