import numpy as np

def solve_final_trace():
    """
    This function calculates the final trace T based on the derived properties and formulas.
    """
    # Step 1: Define group properties determined from the visualizations.
    # R_1i: Sum of character table entries = 8
    # R_2i: Number of irreducible representations = 5
    # R_3i: Order of the group = 8
    # R_4i: Exponent of the group = 4
    R_values = [8, 5, 8, 4]

    # Step 2: Calculate the constant C.
    sum_R = sum(R_values)
    sum_sq_R = sum(x**2 for x in R_values)
    contraharmonic_mean = sum_sq_R / sum_R
    C = int(np.floor(contraharmonic_mean))

    # Step 3: Define coefficients for polynomials Q(x) and S(x).
    # Q(x) = C * (-x^2 + x^4 - x^6 + x^8)
    # S(x) = C * (x - x^3 + x^5 - x^7 + x^9)
    # Coefficients q_k for x^k and s_k for x^k are stored in arrays.
    q_coeffs = np.zeros(11)
    s_coeffs = np.zeros(11)
    
    q_coeffs[2] = -C; q_coeffs[4] = C; q_coeffs[6] = -C; q_coeffs[8] = C
    s_coeffs[1] = C; s_coeffs[3] = -C; s_coeffs[5] = C; s_coeffs[7] = -C; s_coeffs[9] = C
    
    # Step 4: Calculate the trace of the Bezout matrix M1.
    # The Bezout matrix M1 is 9x9 (indexed 0 to 8).
    # Diagonal entry b_ii = sum_k (s_{i-k}*q_{i+k+1} - q_{i-k}*s_{i+k+1})
    trace_M1 = 0
    for i in range(9):
        b_ii = 0
        for k in range(i + 1):
            if (i + k + 1) < len(q_coeffs):
                 term = s_coeffs[i-k] * q_coeffs[i+k+1] - q_coeffs[i-k] * s_coeffs[i+k+1]
                 b_ii += term
        trace_M1 += b_ii

    # Step 5: Calculate the trace of the Sylvester matrix M2.
    # M2 = Sm(Q(x), x^10 + S(x)). deg(Q)=8, deg(x^10+S)=10.
    # Trace(M2) = deg(U)*lead_coeff(Q) + deg(Q)*const_term(U)
    deg_Q = 8
    deg_U = 10
    lead_coeff_Q = q_coeffs[8]
    const_term_U = s_coeffs[0]  # x^10 has 0 const term, S(x) has 0 const term
    trace_M2 = deg_U * lead_coeff_Q + deg_Q * const_term_U

    # Step 6: Compute the final trace T.
    # T = Tr(M1 x I2 + M2) = 2*Tr(M1) + Tr(M2)
    T = 2 * trace_M1 + trace_M2
    
    print(f"The constant C is: {C}")
    print(f"Trace(M1) = {trace_M1}")
    print(f"Trace(M2) = {trace_M2}")
    print(f"The final equation is: T = 2 * {int(trace_M1)} + {int(trace_M2)}")
    print(f"The final trace T is: {int(T)}")
    
solve_final_trace()

<<<-372>>>