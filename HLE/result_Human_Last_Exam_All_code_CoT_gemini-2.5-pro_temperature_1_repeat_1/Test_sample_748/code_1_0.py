import sympy

def analyze_continuity(M, M_k, k, point_name):
    """
    Analyzes the continuity of the minimal polynomial map at a given point M.

    Args:
        M (sympy.Matrix): The point (matrix) at which to check continuity.
        M_k (sympy.Matrix): A sequence of matrices converging to M as k -> oo.
        k (sympy.Symbol): The symbolic variable for the sequence.
        point_name (str): A name for the matrix M for printing.
    """
    # Use sympy's variable X for the polynomial
    X = sympy.Symbol('X')
    
    # Calculate the minimal polynomial of the limit matrix M
    try:
        pi_M = M.minpoly(X)
    except Exception as e:
        pi_M = f"Could not compute: {e}"

    # Calculate the minimal polynomial of the sequence matrix M_k
    try:
        pi_M_k = M_k.minpoly(X)
    except Exception as e:
        pi_M_k = f"Could not compute: {e}"

    # Calculate the limit of the minimal polynomial of the sequence
    try:
        limit_pi_M_k = sympy.limit(pi_M_k, k, sympy.oo)
    except Exception as e:
        limit_pi_M_k = f"Could not compute: {e}"

    # Check for continuity
    is_continuous = (pi_M == limit_pi_M_k)

    # Print the results in a clear format
    print(f"--- Analysis at Point {point_name} ---")
    print(f"M = \n{sympy.pretty(M)}\n")
    print(f"Minimal polynomial of M: π_M(X) = {sympy.expand(pi_M)}")
    print("-" * 20)
    print(f"Sequence M_k = \n{sympy.pretty(M_k)}\n")
    print(f"Minimal polynomial of M_k: π_M_k(X) = {sympy.expand(pi_M_k)}")
    print(f"Limit of π_M_k(X) as k -> ∞: lim(π_M_k) = {sympy.expand(limit_pi_M_k)}")
    print("-" * 20)
    print(f"Is lim(π_M_k) == π_M? {is_continuous}")
    if is_continuous:
        print("Conclusion: The map appears to be CONTINUOUS at this point.")
    else:
        print("Conclusion: The map appears to be DISCONTINUOUS at this point.")
    print("\n" + "="*40 + "\n")


def main():
    """
    Main function to run the demonstrations.
    """
    # Define a symbolic variable k for our sequences
    k = sympy.Symbol('k', positive=True, real=True)
    
    # --- Case 1: Point of Discontinuity (Derogatory Matrix) ---
    # M_d is a scalar matrix. It is derogatory.
    # Its minimal polynomial has degree 1.
    M_d = sympy.Matrix([[2, 0], [0, 2]])
    
    # M_d_k is a sequence of diagonalizable matrices with distinct eigenvalues
    # that converges to M_d as k -> infinity.
    M_d_k = sympy.Matrix([[2 + 1/k, 0], [0, 2 - 1/k]])
    
    analyze_continuity(M_d, M_d_k, k, "M_d (Derogatory)")

    # --- Case 2: Point of Continuity (Non-derogatory Matrix) ---
    # M_c is a non-diagonalizable Jordan block. It is non-derogatory.
    # Its minimal polynomial is equal to its characteristic polynomial, degree 2.
    M_c = sympy.Matrix([[2, 1], [0, 2]])
    
    # M_c_k is a sequence of matrices converging to M_c.
    # For any finite k, the eigenvalues are distinct (2 +/- 1/sqrt(k)),
    # so M_c_k is non-derogatory.
    M_c_k = sympy.Matrix([[2, 1], [1/k, 2]])
    
    analyze_continuity(M_c, M_c_k, k, "M_c (Non-derogatory)")

if __name__ == "__main__":
    main()
