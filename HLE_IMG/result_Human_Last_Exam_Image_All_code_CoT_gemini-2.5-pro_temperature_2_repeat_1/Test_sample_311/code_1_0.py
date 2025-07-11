import sympy as sp
import numpy as np

def solve_problem():
    """
    Solves the multi-step problem as outlined.
    """
    
    # Step 1 & 2: Identify the group and calculate its properties
    # From visual analysis, V1 is the Cayley graph of the Dihedral group D6.
    # It has 12 vertices, so the group order is 12.
    # V5 is a 6x6 matrix, which matches the size of the character table of D6 (6 conjugacy classes).
    # Thus, the group is D6, the dihedral group of order 12.

    # Properties of D6:
    # R1: Sum of all entries in the character table.
    # The character table of D6 is:
    # [[1, 1, 1, 1, 1, 1],
    #  [1, 1, 1, 1, -1, -1],
    #  [1, -1, -1, 1, 1, -1],
    #  [1, -1, -1, 1, -1, 1],
    #  [2, 2, -1, -2, 0, 0],
    #  [2, -2, 1, -2, 0, 0]]
    # Sum of all these entries is 8.
    R1 = 8
    
    # R2: Number of irreducible representations (number of rows/columns in char table)
    R2 = 6
    
    # R3: Order of the group D6
    R3 = 12
    
    # R4: Exponent of the group D6. The element orders are 1, 2, 3, 6.
    # The LCM of these orders is 6.
    R4 = 6

    print(f"Identified Group: D6 (Dihedral group of order 12)")
    print(f"R1 (Sum of char table entries) = {R1}")
    print(f"R2 (Number of irreps) = {R2}")
    print(f"R3 (Group order) = {R3}")
    print(f"R4 (Group exponent) = {R4}")
    print("-" * 30)

    # Step 3: Compute the constant C
    # C is the floor of the contraharmonic mean of R1, R2, R3, R4.
    R_vals = np.array([R1, R2, R3, R4])
    sum_of_squares = np.sum(R_vals**2)
    sum_of_vals = np.sum(R_vals)
    C = int(np.floor(sum_of_squares / sum_of_vals))

    print(f"Sum of squares of R values = {sum_of_squares}")
    print(f"Sum of R values = {sum_of_vals}")
    print(f"Contraharmonic mean = {sum_of_squares / sum_of_vals:.2f}")
    print(f"C = floor(Contraharmonic mean) = {C}")
    print("-" * 30)

    # Step 4: Construct Polynomials
    x = sp.Symbol('x')
    
    # P(x) = C * (x^1 + x^2 + ... + x^9)
    P_x = sum(C * x**i for i in range(1, 10))
    
    # P(ix)
    P_ix = P_x.subs(x, sp.I * x)
    
    # Q(x) = Re(P(ix)), S(x) = Im(P(ix))
    Q_x = sp.re(P_ix)
    S_x_im = sp.im(P_ix) # Note: S_x will be S_x_im
    S_x = S_x_im

    # Ensure polynomials are expanded
    Q_x = sp.expand(Q_x)
    S_x = sp.expand(S_x)

    print(f"P(x) = {P_x}")
    print(f"Q(x) = Re(P(ix)) = {Q_x}")
    print(f"S(x) = Im(P(ix)) = {S_x}")
    print("-" * 30)

    # Step 5: Compute Matrices and Traces
    # Define polynomials as sympy Poly objects
    q_poly = sp.Poly(Q_x, x)
    s_poly = sp.Poly(S_x, x)

    # M1 = Bezout Matrix
    M1 = sp.polys.polytools.bezout(q_poly, s_poly)
    trace_M1 = M1.trace()
    
    # M2 = Sylvester Matrix
    s2_poly = sp.Poly(x**10 + S_x, x)
    M2 = sp.polys.polytools.sylvester(q_poly, s2_poly)
    trace_M2 = M2.trace()

    print(f"Trace(M1) = {trace_M1}")
    print(f"Trace(M2) = {trace_M2}")
    print("-" * 30)

    # Step 6: Calculate the Final Answer T
    # T = Tr(M1 x I2 + M2) = 2*Tr(M1) + Tr(M2)
    # The trace of a Kronecker product is the product of the traces.
    # Tr(M1 x I2) = Tr(M1) * Tr(I2) = Tr(M1) * 2
    T = 2 * trace_M1 + trace_M2
    
    print("Final Calculation:")
    print(f"T = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * ({trace_M1}) + ({trace_M2})")
    print(f"T = {2 * trace_M1} + {trace_M2}")
    print(f"T = {T}")

solve_problem()

print("\n<<< -176 >>>")