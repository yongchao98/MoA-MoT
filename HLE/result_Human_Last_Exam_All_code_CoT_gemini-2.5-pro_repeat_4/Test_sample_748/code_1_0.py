import sympy

def demonstrate_discontinuity():
    """
    This function demonstrates the discontinuity of the map M -> pi_M
    at a derogatory matrix.
    """
    # Define the symbol for the polynomial variable
    x = sympy.Symbol('X')
    n = 3

    # 1. Define a derogatory matrix M0.
    # The 3x3 identity matrix is derogatory because its minimal polynomial,
    # X - 1, has degree 1, which is less than n=3.
    M0 = sympy.eye(n)
    
    # 2. Compute and display the minimal polynomial of M0.
    pi_M0 = M0.minpoly(x)

    print("Step 1: Analyze a derogatory matrix M0.")
    print("Let M0 be the 3x3 identity matrix:")
    sympy.pprint(M0)
    print(f"\nThe minimal polynomial of M0 is pi_M0(X) = {sympy.pretty(pi_M0)}.")
    print("-" * 50)

    # 3. Define a sequence of non-derogatory matrices Mk that converges to M0.
    # We choose Mk to be a Jordan block with eigenvalue 1 and off-diagonal 
    # elements 1/k. As k -> infinity, 1/k -> 0, so Mk -> M0.
    print("Step 2: Construct a sequence of matrices Mk converging to M0.")
    print("Let's define Mk as follows:")
    k_sym = sympy.Symbol('k')
    Mk_symbolic = sympy.Matrix([[1, 1/k_sym, 0], [0, 1, 1/k_sym], [0, 0, 1]])
    sympy.pprint(Mk_symbolic)
    print("\nAs k -> infinity, Mk converges to M0.")
    print("-" * 50)
    
    # 4. Analyze the minimal polynomials of Mk for a few values of k.
    print("Step 3: Analyze the minimal polynomials of Mk.")
    # For any k > 0, Mk is a single Jordan block, so it is non-derogatory.
    # Its minimal polynomial is equal to its characteristic polynomial, (X-1)^3.
    pi_Mk_general = Mk_symbolic.minpoly(x)
    
    print("For any k, the minimal polynomial of Mk is:")
    print(f"pi_Mk(X) = {sympy.pretty(pi_Mk_general)}")
    print(f"Equation: {sympy.pretty(pi_Mk_general)} = 0")
    print("\nThis polynomial does not depend on k and is not equal to pi_M0(X).")
    print("-" * 50)

    # 5. Conclude by comparing the limits.
    print("Step 4: Compare the limits to show discontinuity.")
    print(f"The map value at the limit point M0 is theta(M0) = pi_M0(X) = {sympy.pretty(pi_M0)}")
    print(f"The limit of the map values along the sequence is lim_{{k->inf}} theta(Mk) = {sympy.pretty(pi_Mk_general)}")
    print("\nSince lim_{k->inf} theta(Mk) != theta(M0), the map is discontinuous at M0.")


if __name__ == '__main__':
    demonstrate_discontinuity()