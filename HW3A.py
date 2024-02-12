import math
import copy
import random

# Code below written and checked with help of Chat GPT
# Following code written to check matrixes for symmetric or positive definite, and use stated methods to solve.
def print_matrix(matrix):
    """
    Prints the elements of a matrix
    :param matrix: The matrix that it prints
    :return: Nothing
    """
    for row in matrix:
        print(row)

def transpose(matrix):
    """
    Transposes a matrix
    :param matrix: The matrix that is transposed
    :return: The transposed matrix
    """
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def matrix_multiply(matrix1, matrix2):
    """
    Multiplies two matrices
    :param matrix1: Matrix 1
    :param matrix2: Matrix 2
    :return: The Resulting Matrix
    """
    result = [[0 for _ in range(len(matrix2[0]))] for _ in range(len(matrix1))]
    for i in range(len(matrix1)):
        for j in range(len(matrix2[0])):
            for k in range(len(matrix2)):
                result[i][j] += matrix1[i][k] * matrix2[k][j]
    return result

def is_symmetric(matrix):
    """
    Checks if the matrix is symmetric
    :param matrix: The Matrix that is checked
    :return: True or False statement depending on if matrix is symmetric
    """
    return all(matrix[i][j] == matrix[j][i] for i in range(len(matrix)) for j in range(len(matrix[0])))

def is_positive_definite(matrix):
    """
    Checks if the matrix is positive definite
    :param matrix: Matrix that is checked
    :return: True or False depending on if matrix is positive definite
    """
    try:
        for i in range(len(matrix)):
            submatrix = [row[:i+1] for row in matrix[:i+1]]
            if determinant(submatrix) <= 0:
                return False
    except ValueError:
        return False
    return True

def cholesky_decomposition(matrix):
    """
    Preforms Cholesky decomposition of matrix
    :param matrix: Matrix to go through Cholesky Decomposition
    :return: Lower traingular matrix of the decompisition
    """
    n = len(matrix)
    L = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i+1):
            if i == j:
                sum_val = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = math.sqrt(matrix[i][i] - sum_val)
            else:
                sum_val = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (matrix[i][j] - sum_val) / L[j][j]
    return L

def forward_substitution(L, b):
    """
    Preforms forward subsitution to solve a lower triangle matrix
    :param L: Lower traingular matrix
    :param b: Column vectors of constants
    :return: The Solution vector
    """
    n = len(L)
    y = [0] * n
    for i in range(n):
        sum_val = sum(L[i][j] * y[j] for j in range(i))
        y[i] = (b[i][0] - sum_val) / L[i][i]
    return y


def backward_substitution(U, y):
    """
    Preforms backward substitution to solve an upper triangle matrix
    :param U: Upper triangular matrix
    :param y: The solution vector obtained from forward subsitutions
    :return: The Solution vector
    """
    n = len(U)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(U[i][j] * x[j] for j in range(i + 1, n))) / U[i][i]
    return x

def determinant(matrix):
    """
    Calculates the determinate of a matrix
    :param matrix: The Matrix to be caculated
    :return: The Determinate
    """
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    elif n == 2:
        return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
    else:
        det = 0
        for j in range(n):
            minor = [row[:j] + row[j+1:] for row in matrix[1:]]
            det += (-1)**j * matrix[0][j] * determinant(minor)
        return det

def doolittle_decomposition(matrix):
    """
    Preforms Doolittle decomposition of a matrix
    :param matrix: The Matrix to go through Doolittle decomposition
    :return: A tuple containing the lower triangle matrix and the upper triangle matrix
    """
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i, n):
            sum_val = sum((L[i][k] * U[k][j]) for k in range(i))
            U[i][j] = matrix[i][j] - sum_val

        for j in range(i, n):
            if i == j:
                L[i][j] = 1.0
            else:
                sum_val = sum((L[j][k] * U[k][i]) for k in range(i))
                L[j][i] = (matrix[j][i] - sum_val) / U[i][i]

    return L, U

def solve_cholesky(matrix, b):
    """
    Solves a system of linear equations using Cholesky decomp.
    :param matrix: The Coefficient matrix
    :param b: The Column Vector of constants
    :return: The solution vector
    """
    L = cholesky_decomposition(matrix)
    y = forward_substitution(L, b)
    x = backward_substitution(transpose(L), y)
    return x

def solve_doolittle(matrix, b):
    """
    Solves a system of linear equations using Doolittles decomp.
    :param matrix: The Coefficient matrix
    :param b: The columns vector of constants
    :return: The solution vector
    """
    L, U = doolittle_decomposition(matrix)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return x

def main():
    """
    The main function to demostrate solving linear equations using both methods
    :return:Solution
    """
    # Problem 1 from HW3A
    A1 = [[1, -1, 3, 2], [-1, 5, -5, -2], [3, -5, 19, 3], [2, -2, 3, 21]]
    b1 = [15], [-35], [94], [1]
    print("Problem 1:")
    print("Matrix A1:")
    print_matrix(A1)
    print("Vector b1:", b1)

    if is_symmetric(A1) and is_positive_definite(A1):
        print("Matrix A1 is symmetric and positive definite. Using Cholesky method.")
        x1 = solve_cholesky(A1, b1)
    else:
        print("Matrix A1 is not symmetric and positive definite. Using Doolittle method.")
        x1 = solve_doolittle(A1, b1)

    print("Solution vector x1:", x1)

    # Problem 2 from HW3A
    A2 = [[4, 2, 4, 0], [2, 3, 3, 2], [4, 3, 6, 3], [0, 2, 3, 9]]
    b2 = [20], [36], [60], [122]

    print("\nProblem 2:")
    print("Matrix A2:")
    print_matrix(A2)
    print("Vector b2:", b2)

    if is_symmetric(A2) and is_positive_definite(A2):
        print("Matrix A2 is symmetric and positive definite. Using Cholesky method.")
        x2 = solve_cholesky(A2, b2)
    else:
        print("Matrix A2 is not symmetric and positive definite. Using Doolittle method.")
        x2 = solve_doolittle(A2, b2)

    print("Solution vector x2:", x2)

if __name__ == "__main__":
    main()
