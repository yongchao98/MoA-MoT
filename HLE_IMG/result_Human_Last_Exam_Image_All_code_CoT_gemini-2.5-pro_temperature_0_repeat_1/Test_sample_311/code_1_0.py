import sympy
import numpy as np

def solve():
    """
    This function solves the problem by following the outlined plan.
    """
    # Step 1 & 2: Define the R values for the group D_4h.
    # R3: Order of the group D4h = D4 x C2 is 8 * 2 = 16.
    R3 = 16
    # R2: Number of conjugacy classes of D4h is num_classes(D4) * num_classes(C2) = 5 * 2 = 10.
    R2 = 10
    # R4: Exponent of D4h is lcm(exponent(D4), exponent(C2)) = lcm(4, 2) = 4.
    R4 = 4
    # R1: Sum of character table entries is the number of elements g with g^2=e.
    # In D4h, these are the identity (1) and elements of order 2 (11). So, 1 + 11 = 12.
    R1 = 12

    # Step 3: Calculate the constant C.
    R_values = [R1, R2, R3, R4]
    r_sq_sum = sum(r**2 for r in R_values)
    r_sum = sum(R_values)
    C = int(np.floor(r_sq_sum / r_sum))

    # Step 4: Construct the polynomials P(x), Q(x), and S(x).
    x = sympy.Symbol('x')
    P = sum(C * x**i for i in range(1, 10))
    
    # Substitute x with i*x and separate real and imaginary parts.
    P_ix = P.subs(x, sympy.I * x)
    Q = sympy.expand(sympy.re(P_ix))
    S = sympy.expand(sympy.im(P_ix))

    # Step 5: Construct matrices M1 and M2.
    # Convert symbolic expressions to Poly objects for matrix functions.
    Q_poly = sympy.Poly(Q, x)
    S_poly = sympy.Poly(S, x)
    S2_poly = sympy.Poly(x**10 + S, x)

    # M1 is the Bezout matrix of Q and S.
    M1_sympy = sympy.polys.matrices.bezout(Q_poly, S_poly)
    M1 = M1_sympy.to_Matrix()

    # M2 is the Sylvester matrix of Q and (x^10 + S).
    M2_sympy = sympy.polys.matrices.sylvester(Q_poly, S2_poly)
    M2 = M2_sympy.to_Matrix()

    # Step 6: Calculate the traces.
    trace_M1 = sympy.trace(M1)
    trace_M2 = sympy.trace(M2)

    # Step 7: Compute the final result T.
    T = 2 * trace_M1 + trace_M2

    # Print the required intermediate and final values.
    print(f"The four characteristic numbers for the group D4h are:")
    print(f"R1 (Sum of character table entries) = {R1}")
    print(f"R2 (Number of irreducible representations) = {R2}")
    print(f"R3 (Order of the group) = {R3}")
    print(f"R4 (Exponent of the group) = {R4}")
    print(f"The constant C is floor(({R1}^2 + {R2}^2 + {R3}^2 + {R4}^2) / ({R1} + {R2} + {R3} + {R4})) = {C}")
    print(f"The trace of M1 is: {trace_M1}")
    print(f"The trace of M2 is: {trace_M2}")
    print(f"The final equation is T = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * {trace_M1} + {trace_M2} = {T}")
    
    return T

result = solve()
print(f"\nFinal Answer: {result}")
<<<128>>>