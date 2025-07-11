import numpy as np

def solve_determinant():
    """
    Computes the determinant of the given 3x3 matrix and prints the calculation.
    """
    # The adjacency matrix is given as A.
    # The value m_1/2 = 13 is not needed for this calculation.
    A = np.array([
        [2, -1, -1],
        [-1, 2, -1],
        [-1, -1, 2]
    ])

    # Extract elements for the formula
    a, b, c = A[0, 0], A[0, 1], A[0, 2]
    d, e, f = A[1, 0], A[1, 1], A[1, 2]
    g, h, i = A[2, 0], A[2, 1], A[2, 2]

    # Calculate the determinant using the formula.
    # We can also use np.linalg.det(A) for a direct answer.
    determinant = np.linalg.det(A)
    
    # Per the user's request, we will output each number in the final equation.
    # The formula is: a(ei - fh) - b(di - fg) + c(dh - eg)
    # The exact value of the determinant is an integer, so we cast it.
    final_result = int(round(determinant))
    
    # We construct the string representing the equation.
    equation = (
        f"{a} * ({e}*{i} - ({f})*{h}) - "
        f"({b}) * ({d}*{i} - ({f})*{g}) + "
        f"({c}) * ({d}*{h} - {e}*{g}) = {final_result}"
    )

    print("The determinant of the adjacency matrix is computed as follows:")
    print(equation)

solve_determinant()
<<<0>>>