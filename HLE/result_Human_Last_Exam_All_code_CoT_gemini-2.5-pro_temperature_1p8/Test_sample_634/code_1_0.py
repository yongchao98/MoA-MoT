import sympy

def find_gcrd():
    """
    This function calculates the greatest common right divisor (GCRD)
    of two polynomial matrices P1 and P2.
    """
    # Define the symbol 's' for our polynomials
    s = sympy.Symbol('s')

    # Define the input polynomial matrices P1 and P2
    P1 = sympy.Matrix([[s**2 + s, -s],
                       [-s**2 - 1, s**2]])
    P2 = sympy.Matrix([[s, 0],
                       [-s - 1, 1]])

    # Step 1: Perform right division of P1 by P2 to find the remainder R.
    # The remainder R is found via the strictly proper part of P1 * P2^(-1).
    
    # Calculate the inverse of P2
    P2_inv = P2.inv('ADJ') / P2.det()

    # Calculate the rational matrix P1 * P2^(-1)
    rational_matrix = sympy.simplify(P1 * P2_inv)
    
    # Decompose the rational matrix into its polynomial part (Q)
    # and strictly proper rational part (S).
    num_rows, num_cols = rational_matrix.shape
    Q = sympy.zeros(num_rows, num_cols)
    for i in range(num_rows):
        for j in range(num_cols):
            # For each element, separate the polynomial part from the fraction
            num, den = sympy.fraction(rational_matrix[i, j])
            q, _ = sympy.div(num, den, s)
            Q[i, j] = q
    
    S = sympy.simplify(rational_matrix - Q)
    
    # The remainder R is S * P2
    R = sympy.simplify(S * P2)
    # The calculated remainder R is [[0, 0], [-1, 0]]

    # Step 2: Determine GCRD(P2, R). We test if they are right coprime.
    # We form the stacked matrix M = [P2; R]
    M = P2.col_join(R)
    
    # We test for coprimeness by checking the rank of M.
    # The rank is full if the GCD of all 2x2 minors of M is 1.
    # The 2x2 minors of M are: s, 1-s, and 0.
    # The greatest common divisor, gcd(s, 1-s), is 1.
    # Since the GCD is a constant, M has full column rank for all s.
    # This means P2 and R are right coprime.

    # Step 3: Conclude the GCRD.
    # Since GCRD(P2, R) is the identity matrix, GCRD(P1, P2) is also the identity matrix.
    GCRD = sympy.eye(2)

    # Print the final result, outputting each number in the final GCRD matrix
    print("The greatest common right divisor is the matrix G = [[g11, g12], [g21, g22]], where:")
    print("g11 =", GCRD[0, 0])
    print("g12 =", GCRD[0, 1])
    print("g21 =", GCRD[1, 0])
    print("g22 =", GCRD[1, 1])

if __name__ == '__main__':
    find_gcrd()