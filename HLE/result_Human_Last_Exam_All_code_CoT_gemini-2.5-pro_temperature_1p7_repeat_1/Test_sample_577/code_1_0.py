import sympy

def solve_chern_class():
    """
    Computes the total Chern class of the bundle S^{(2,1)}T(CP^4).
    """
    # h is the generator of the cohomology ring of CP^4
    h = sympy.Symbol('h')

    # Step 1 & 2: Define ch(T) using K-theory relation T = 5*O(1) - O
    # ch(O(k)) = exp(k*h)
    ch_T = 5 * sympy.exp(h) - 1
    
    # We work in the ring Z[h]/<h^5>, so we expand up to degree 4.
    series_limit = 5
    ch_T_series = ch_T.series(h, 0, series_limit).removeO()

    # Step 3 & 4: Define ch(psi^3(T))
    # psi^3(T) = 5*O(3) - O, so ch(psi^3(T)) = 5*exp(3h) - 1
    ch_psi3_T = 5 * sympy.exp(3*h) - 1
    ch_psi3_T_series = ch_psi3_T.series(h, 0, series_limit).removeO()

    # Step 5: Compute ch(S^{(2,1)}T)
    # The formula is ch(S^(2,1)T) = 1/3 * (ch(T)^3 - ch(psi^3 T))
    ch_T_cubed_series = sympy.poly((ch_T_series**3).expand(), h).trunc(series_limit)
    
    ch_S_series = sympy.Rational(1,3) * (ch_T_cubed_series - ch_psi3_T_series)
    ch_S_poly = sympy.Poly(ch_S_series, h).trunc(series_limit)

    # Step 6: Convert Chern character coefficients to Chern classes
    # ch(E) = sum(Ch_k * h^k)
    # The power sums of roots are p_k = k! * Ch_k
    
    Ch = [ch_S_poly.coeff_monomial(h**k) for k in range(series_limit)]
    
    p = [sympy.factorial(k) * Ch[k] for k in range(series_limit)]
    
    # Use Newton's sums to find the Chern classes c_k
    # c_k are coefficients of h^k in the total Chern class polynomial
    c = {}
    c[0] = 1 # c_0 is always 1
    
    # c_1 = p_1
    c[1] = p[1]
    
    # 2*c_2 = c_1*p_1 - p_2
    c[2] = (c[1] * p[1] - p[2]) / 2
    
    # 3*c_3 = c_2*p_1 - c_1*p_2 + p_3
    c[3] = (c[2] * p[1] - c[1] * p[2] + p[3]) / 3
    
    # 4*c_4 = c_3*p_1 - c_2*p_2 + c_1*p_3 - p_4
    c[4] = (c[3] * p[1] - c[2] * p[2] + c[1] * p[3] - p[4]) / 4

    # Step 7: Print the final polynomial
    print("The total Chern class of S^{(2,1)}T(CP^4) is:")
    
    # Constructing the output string ensuring integers are displayed correctly.
    # Each coefficient c[k] is guaranteed to be an integer.
    c_int = {k: int(v) for k, v in c.items()}
    
    # The problem asks to output each number in the final equation.
    # e.g., print(f"c_1 = {c_int[1]}") etc. to build the string
    # For a clean final answer we build the polynomial string directly.
    final_poly = (f"{c_int[0]} + {c_int[1]}*h + {c_int[2]}*h^2 + "
                  f"{c_int[3]}*h^3 + {c_int[4]}*h^4")
    print(final_poly)

    # Output for the final answer block
    return final_poly

if __name__ == '__main__':
    result = solve_chern_class()
    # The required format is <<<answer content>>> at the very end.
    # We will print the computed polynomial in this format.
    print(f"\n<<<{result}>>>")