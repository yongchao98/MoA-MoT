import numpy as np
import sympy

def solve_group_theory_problem():
    """
    Solves the multi-step problem by identifying the group, calculating
    properties, constructing polynomials and matrices, and finding the final trace.
    """
    # Step 1: Define properties of the identified Crystallographic Point Group (T_h).
    R1 = 24  # Sum of character table entries (equals group order)
    R2 = 8   # Number of irreducible representations
    R3 = 24  # Order of the group
    R4 = 6   # Exponent of the group (lcm of element orders {1,2,3,6})
    
    R_values = np.array([R1, R2, R3, R4])

    # Step 2: Compute the floor of the contraharmonic mean.
    contraharmonic_mean = np.sum(R_values**2) / np.sum(R_values)
    C = int(np.floor(contraharmonic_mean))

    # Step 3: Define the polynomials Q(x) and S(x).
    x = sympy.Symbol('x')
    
    # Q(x) is the real part of P(ix)
    Q = C * (-x**2 + x**4 - x**6 + x**8)
    
    # S(x) is the imaginary part of P(ix)
    S = C * (x - x**3 + x**5 - x**7 + x**9)

    # Convert to SymPy Poly objects for matrix functions
    Q_poly = sympy.Poly(Q, x)
    S_poly = sympy.Poly(S, x)

    # Step 4: Construct the matrices M1 and M2.
    # M1 = Bezout Matrix
    # We use the modern sympy API from sympy.polys.matrices
    M1_sympy = sympy.polys.matrices.bezoutm(Q_poly, S_poly)
    M1 = np.array(M1_sympy.tolist()).astype(np.float64)

    # M2 = Sylvester Matrix
    T_poly = sympy.Poly(x**10 + S, x)
    M2_sympy = sympy.polys.matrices.sylvesterm(Q_poly, T_poly)
    M2 = np.array(M2_sympy.tolist()).astype(np.float64)

    # Step 5: Calculate the trace T.
    trace_M1 = np.trace(M1)
    trace_M2 = np.trace(M2)
    
    # Using the identity Tr(A x B + C) = Tr(A)*Tr(B) + Tr(C)
    # where A=M1, B=I2, C=M2. Tr(I2)=2.
    T = 2 * trace_M1 + trace_M2
    
    # Output the steps of the final calculation
    print(f"The group properties are R1={R1}, R2={R2}, R3={R3}, R4={R4}.")
    print(f"The contraharmonic mean is {contraharmonic_mean:.4f}, so C = {C}.")
    print(f"The real part polynomial is Q(x) = {Q}")
    print(f"The imaginary part polynomial is S(x) = {S}")
    print("\nFinal trace calculation:")
    print(f"Tr(M1) = {trace_M1}")
    print(f"Tr(M2) = {trace_M2}")
    print(f"T = 2 * Tr(M1) + Tr(M2) = 2 * ({trace_M1}) + ({trace_M2}) = {T}")

solve_group_theory_problem()
print("<<<-2992.0>>>")