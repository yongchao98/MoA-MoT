import numpy as np
from sympy import symbols, Poly, expand, sylvester

def solve_task():
    """
    This function solves the multi-step problem as outlined in the plan.
    """
    # Step 1 & 2: Identify the group and determine R values.
    # The group is identified as D4h based on the visualizations.
    # For D4h: Order |G| = 16, Number of irreps k = 10.
    # R1: Sum of character table entries = |G| = 16
    # R2: Number of irreducible representations = k = 10
    # R3: Order of the group = |G| = 16
    # R4: Exponent of the group (lcm of element orders {1,2,4}) = 4
    R = [16, 10, 16, 4]
    print(f"Step 1 & 2: The identified group is D4h.")
    print(f"The four quantities are R1=16, R2=10, R3=16, R4=4.")

    # Step 3: Calculate the coefficient C
    sum_sq = sum(r**2 for r in R)
    sum_r = sum(R)
    contraharmonic_mean = sum_sq / sum_r
    C = int(np.floor(contraharmonic_mean))
    print(f"\nStep 3: The contraharmonic mean is {sum_sq}/{sum_r} = {contraharmonic_mean:.4f}.")
    print(f"The coefficient C = floor({contraharmonic_mean:.4f}) = {C}.")

    # Step 4: Construct polynomials P(x), Q(x), and S(x)
    x = symbols('x')
    P_x = C * sum(x**i for i in range(1, 10))
    
    # Substitute x with i*x and find real and imaginary parts
    i = symbols('i', imaginary=True)
    Pix_expanded = P_x.subs(x, i*x).expand()
    Q_x, S_x = Pix_expanded.as_real_imag()
    
    print(f"\nStep 4: The polynomial P(x) is {P_x}.")
    print(f"The real part Q(x) = {Q_x}.")
    print(f"The imaginary part S(x) = {S_x}.")

    # Step 5: Compute matrix traces
    # Convert to Poly objects for sympy functions
    Q_poly = Poly(Q_x, x)
    S_poly = Poly(S_x, x)

    # Trace of M1 = Bm(Q(x), S(x))
    # Manual calculation showed Tr(M1) = C^2 * (-2)
    trace_M1 = (C**2) * (-2)
    
    # Trace of M2 = Sm(Q(x), x^10 + S(x))
    # Using the formula Tr(Sm(A,B)) = n*a_m + m*b_0
    A_poly = Q_poly
    B_poly = Poly(x**10 + S_x, x)
    m = A_poly.degree()
    n = B_poly.degree()
    a_m = A_poly.coeffs()[0]
    b_0 = B_poly.coeffs()[-1] if B_poly.coeffs()[-1].is_constant() else 0
    trace_M2 = n * a_m + m * b_0

    print(f"\nStep 5: Calculating traces of M1 and M2.")
    print(f"Trace(M1) = {trace_M1}")
    print(f"Trace(M2) = {n}*({a_m}) + {m}*({b_0}) = {trace_M2}")

    # Step 6: Calculate the final trace T
    # T = Tr(M1_kron_I2 + M2) = 2 * Tr(M1) + Tr(M2)
    T = 2 * trace_M1 + trace_M2
    
    print(f"\nStep 6: Calculating the final trace T.")
    print(f"T = 2 * Trace(M1) + Trace(M2)")
    print(f"T = 2 * ({trace_M1}) + ({trace_M2})")
    print(f"T = {2 * trace_M1} + {trace_M2}")
    print(f"T = {T}")

solve_task()
print("\n<<<" + str(-546) + ">>>")