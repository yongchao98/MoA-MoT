import math
import sympy

def solve():
    """
    This function performs the calculations based on the analysis of the problem.
    """
    # Step 1 & 2: Define the properties of the identified group, D5.
    R = [10, 4, 10, 10]
    print(f"The four characteristic quantities R_j for the group D5 are: {R}")

    # Step 3: Calculate the coefficient C
    sum_R = sum(R)
    sum_R_sq = sum(x**2 for x in R)
    contraharmonic_mean = sum_R_sq / sum_R
    C = math.floor(contraharmonic_mean)
    print(f"The contraharmonic mean is {sum_R_sq}/{sum_R} = {contraharmonic_mean:.4f}")
    print(f"The coefficient C = floor({contraharmonic_mean:.4f}) = {C}")
    print("-" * 20)

    # Step 4: Define the polynomials Q(x) and S(x)
    x = sympy.Symbol('x')
    P = C * sum(x**i for i in range(1, 10))
    Pix = P.subs(x, sympy.I * x).expand()
    Q = sympy.re(Pix)
    S = sympy.im(Pix)
    print(f"P(x) = {P}")
    print(f"Q(x) = Re(P(ix)) = {Q}")
    print(f"S(x) = Im(P(ix)) = {S}")
    print("-" * 20)

    # Step 5: Calculate the traces of M1 and M2
    
    # Trace of M2 = Sylvester(Q, x^10 + S)
    # deg(Q) = 8, deg(x^10+S) = 10
    # lead_coeff(Q) = 9, const_coeff(x^10+S) = 0
    tr_M2 = 10 * C + 8 * 0
    print(f"Trace(M2) = 10 * (lead_coeff(Q)) + 8 * (const_coeff(x^10+S))")
    print(f"Trace(M2) = 10 * {C} + 8 * 0 = {tr_M2}")

    # Trace of M1 = Bezout(Q, S)
    # We use the formula for the trace of a Bezout matrix.
    # Tr(Bez(Q0, S0)) = -Tr(Bez(S0, Q0))
    # Tr(Bez(S0, Q0)) = sum_{k=0 to n-1} (s_{k+1}q_k - s_k q_{k+1})
    # where S0 = S/C and Q0 = Q/C
    Q0_poly = sympy.Poly(Q / C, x)
    S0_poly = sympy.Poly(S / C, x)
    
    q_coeffs = Q0_poly.all_coeffs()
    s_coeffs = S0_poly.all_coeffs()
    
    # Pad with zeros to make them same length for easier indexing
    deg_s = S0_poly.degree()
    deg_q = Q0_poly.degree()
    
    q_coeffs_map = {deg_q - i: c for i, c in enumerate(q_coeffs)}
    s_coeffs_map = {deg_s - i: c for i, c in enumerate(s_coeffs)}

    trace_bez_S0_Q0 = 0
    for k in range(deg_s): # deg_s = 9, so k from 0 to 8
        s_k_plus_1 = s_coeffs_map.get(k + 1, 0)
        q_k = q_coeffs_map.get(k, 0)
        s_k = s_coeffs_map.get(k, 0)
        q_k_plus_1 = q_coeffs_map.get(k + 1, 0)
        trace_bez_S0_Q0 += (s_k_plus_1 * q_k - s_k * q_k_plus_1)

    # We need Tr(Bez(Q0, S0)) which is the negative of the above
    trace_bez_Q0_S0 = -trace_bez_S0_Q0
    
    tr_M1 = (C**2) * trace_bez_Q0_S0
    print(f"Trace(M1) = C^2 * Trace(Bezout(Q/C, S/C)) = {C**2} * {trace_bez_Q0_S0} = {tr_M1}")
    print("-" * 20)

    # Step 6: Final Calculation
    # T = Tr(M1 kron I2 + M2) = 2 * Tr(M1) + Tr(M2)
    T = 2 * tr_M1 + tr_M2
    
    print("Final Calculation:")
    print(f"T = Tr(M1 x I_2 + M2) = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * ({tr_M1}) + {tr_M2}")
    print(f"T = {2 * tr_M1} + {tr_M2}")
    print(f"T = {T}")
    
    return T

if __name__ == '__main__':
    final_answer = solve()
    print(f"\nFinal Answer: {final_answer}")
    print(f'<<<{final_answer}>>>')
