import sympy
import numpy as np
from sympy.polys.matrices import bezout, sylvester
from sympy import Symbol, Poly, re, im, I

def solve_task():
    """
    Solves the multi-step problem by calculating coefficients, constructing polynomials,
    and finding the trace of the resulting matrix expression.
    """

    # Step 1 & 2: Define the Rji quantities for each group type.
    # Group type 1: Deduced from 10-dimensional visualizations (V1, V6, V7, V9) -> D5
    # R1: Sum of Char Table = #elements of order <= 2 = 1 (e) + 5 (refl) = 6
    # R2: #Irreps = 4
    # R3: Order = 10
    # R4: Exponent = lcm(1, 2, 5) = 10
    R_d5 = [6, 4, 10, 10]

    # Group type 2: Deduced from 16-dimensional visualizations (V2, V4, V5, V8) -> D4h
    # R1: Sum of Char Table = #elements of order <= 2 = 12
    # R2: #Irreps = 10
    # R3: Order = 16
    # R4: Exponent = 4
    R_d4h = [12, 10, 16, 4]

    # Group type 3: Deduced from 4-dimensional visualization (V3) -> D2 (Klein-4)
    # R1: Sum of Char Table = #elements of order <= 2 = 4
    # R2: #Irreps = 4 (abelian group)
    # R3: Order = 4
    # R4: Exponent = 2
    R_d2 = [4, 4, 4, 2]

    R_map = {
        1: R_d5, 2: R_d4h, 3: R_d2, 4: R_d4h, 5: R_d4h,
        6: R_d5, 7: R_d5, 8: R_d4h, 9: R_d5
    }

    # Step 3: Calculate Ci for each visualization
    C = [0] * 10  # Use 1-based indexing
    print("Calculating coefficients C_i:")
    for i in range(1, 10):
        R = R_map[i]
        sum_R = sum(R)
        sum_R_sq = sum(x**2 for x in R)
        contraharmonic_mean = sum_R_sq / sum_R
        C[i] = int(np.floor(contraharmonic_mean))
        print(f"C_{i} = floor( ({R[0]}^2+{R[1]}^2+{R[2]}^2+{R[3]}^2) / ({R[0]}+{R[1]}+{R[2]}+{R[3]}) ) = floor( {sum_R_sq} / {sum_R} ) = {C[i]}")

    # Step 4: Construct Polynomials P(x), Q(x), S(x)
    x = Symbol('x')
    P_poly_terms = [C[i] * x**i for i in range(1, 10)]
    P_poly = sum(P_poly_terms)
    
    # Substitute x with ix to find P(ix)
    P_ix = P_poly.subs(x, I * x)

    # Separate into real (Q) and imaginary (S) parts
    Q_poly = re(P_ix).expand()
    S_poly = im(P_ix).expand()

    print("\nConstructed Polynomials:")
    print(f"P(x) = {P_poly}")
    print(f"Q(x) = RealPart[P(ix)] = {Q_poly}")
    print(f"S(x) = ImaginaryPart[P(ix)] = {S_poly}")
    
    # Step 5: Matrix Computations
    Q_sympy_poly = Poly(Q_poly, x)
    S_sympy_poly = Poly(S_poly, x)

    # Bezout Matrix M1
    M1 = bezout(Q_sympy_poly, S_sympy_poly)
    trace_M1 = M1.trace()

    # Sylvester Matrix M2
    P2_poly = x**10 + S_poly
    P2_sympy_poly = Poly(P2_poly, x)
    
    # We can compute the trace of the Sylvester matrix using a known formula:
    # Tr(Sm(P,Q)) = deg(P)*coeff(Q,deg(Q)) + deg(Q)*coeff(P,deg(P))
    # In our case, P = Q_sympy_poly, Q = P2_sympy_poly
    deg_Q_poly = Q_sympy_poly.degree()
    lead_coeff_Q_poly = Q_sympy_poly.LC()
    deg_P2_poly = P2_sympy_poly.degree()
    lead_coeff_P2_poly = P2_sympy_poly.LC()
    trace_M2 = deg_Q_poly * lead_coeff_P2_poly + deg_P2_poly * lead_coeff_Q_poly

    # Alternatively, construct the matrix and find its trace
    # M2 = sylvester(Q_sympy_poly, P2_sympy_poly)
    # trace_M2 = M2.trace()
    
    print("\nMatrix Trace Calculations:")
    print(f"Tr(M1) = Tr(Bezout[Q(x), S(x)]) = {trace_M1}")
    print(f"Tr(M2) = Tr(Sylvester[Q(x), x^10 + S(x)]) = {trace_M2}")

    # Step 6: Calculate Final Trace T
    T = 2 * trace_M1 + trace_M2

    print("\nFinal Calculation:")
    print(f"T = Tr(M1 x I_2 + M2) = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * {trace_M1} + {trace_M2}")
    print(f"T = {2 * trace_M1} + {trace_M2}")
    print(f"T = {T}")
    
    # Final answer format
    print(f"\n<<<{T}>>>")

solve_task()