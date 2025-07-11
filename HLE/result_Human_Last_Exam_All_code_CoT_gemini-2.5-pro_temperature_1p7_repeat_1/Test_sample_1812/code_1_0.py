import numpy as np

def solve_determinant():
    """
    Computes and explains the determinant of the given matrix.
    """
    # The matrix from the problem description
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Calculate the determinant using numpy
    determinant = np.linalg.det(A)

    print("The given matrix is A:")
    print(A)
    print("\nTo compute the determinant of A, we use the formula:")
    print("det(A) = a(ei - fh) - b(di - fg) + c(dh - eg)\n")
    print("Substituting the values from matrix A:")
    print("det(A) = 2 * ((2 * 2) - (-1 * -1)) - (-1) * ((-1 * 2) - (-1 * -1)) + (-1) * ((-1 * -1) - (2 * -1))")
    print("\nStep-by-step simplification:")
    print("det(A) = 2 * (4 - 1) + 1 * (-2 - 1) - 1 * (1 - (-2))")
    print("det(A) = 2 * (3) + 1 * (-3) - 1 * (3)")
    print("det(A) = 6 - 3 - 3")

    # Print the final result, rounded to the nearest integer
    print(f"\nThe final computed determinant is: {int(round(determinant))}")

solve_determinant()