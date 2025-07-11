import sympy

def solve_energy_shift():
    """
    Calculates the second-order energy shift due to the relativistic kinetic
    energy correction using the Unsöld approximation.
    """
    # Define state quantum numbers
    n = 3
    l = 2

    # Use sympy's rational numbers for precision
    S = sympy.S

    # --- Step 1: Calculate expectation values in atomic units (a_0 = 1) ---
    # Unperturbed energy E_n = -1/(2n^2)
    E_n = -S(1) / (2 * n**2)

    # <1/r> = 1/n^2
    exp_r_neg1 = S(1) / (n**2)

    # <1/r^2> = 1 / (n^3 * (l + 1/2))
    exp_r_neg2 = S(1) / (n**3 * (l + S(1)/2))

    # <1/r^3> = 1 / (n^3 * l * (l + 1) * (l + 1/2))
    exp_r_neg3 = S(1) / (n**3 * l * (l + 1) * (l + S(1)/2))

    # <1/r^4> = [3n^2 - l(l+1)] / [2n^5 * (l-1/2)l(l+1/2)(l+1)(l+3/2)]
    num_r_neg4 = 3 * n**2 - l * (l + 1)
    den_r_neg4 = 2 * n**5 * (l - S(1)/2) * l * (l + S(1)/2) * (l + 1) * (l + S(3)/2)
    exp_r_neg4 = num_r_neg4 / den_r_neg4

    # Expectation values of powers of the potential V = -1/r
    exp_V = -exp_r_neg1
    exp_V2 = exp_r_neg2
    exp_V3 = -exp_r_neg3
    exp_V4 = exp_r_neg4

    # --- Step 2: Calculate <H'> and <(H')^2> in atomic units ---
    # H' = -p^4/(8c^2) = -(H_0 - V)^2 / (2c^2) in a.u. (m=1)
    # The speed of light c in atomic units is 1/alpha
    # We will compute the coefficient C such that <H'>_au = C * alpha^2
    # <H'> = < -(H_0 - V)^2 / (2c^2) > = -alpha^2/2 * <(E_n - V)^2>
    
    # <(E_n - V)^2> = E_n^2 - 2*E_n*<V> + <V^2>
    exp_H_prime_val = E_n**2 - 2 * E_n * exp_V + exp_V2

    # <(H')^2> = < (-(H_0 - V)^2 / (2c^2))^2 > = alpha^4/4 * <(E_n - V)^4>
    # <(E_n - V)^4> = E_n^4 - 4*E_n^3*<V> + 6*E_n^2*<V^2> - 4*E_n*<V^3> + <V^4>
    exp_H_prime_sq_val = (E_n**4 - 4 * E_n**3 * exp_V + 6 * E_n**2 * exp_V2
                            - 4 * E_n * exp_V3 + exp_V4)
                            
    # --- Step 3: Apply the Unsöld approximation ---
    # Delta_E2 ~ (1/E_n) * [<H'^2> - <H'>^2]
    # Delta_E2_au ~ (1/E_n) * [alpha^4/4 * <(E_n-V)^4> - (alpha^2/2 * <(E_n-V)^2>)^2]
    # Delta_E2_au = alpha^4 / (4*E_n) * [<(E_n-V)^4> - <(E_n-V)^2>^2]
    
    delta_E2_coeff = (1 / (4 * E_n)) * (exp_H_prime_sq_val - exp_H_prime_val**2)

    # --- Step 4: Express final result in terms of fundamental constants ---
    # Energy in a.u. is converted to SI by multiplying by E_h = m*c^2*alpha^2
    # The final energy shift is Delta_E2_SI = Delta_E2_au * m*c^2*alpha^2
    # Delta_E2_SI = (delta_E2_coeff * alpha^4) * (m*c^2*alpha^2)
    # Delta_E2_SI = delta_E2_coeff * m*c^2*alpha^6
    
    # Print the result in a readable format
    # The final result is of the form: (numerator/denominator) * m*c^2*alpha^6
    
    final_coeff = sympy.ratsimp(delta_E2_coeff)
    
    print("Based on the Unsöld approximation, the second-order energy shift is:")
    print(f"\u0394E\u2082 \u2248 ({final_coeff.p}/{final_coeff.q}) \u00d7 m\u2091c\u00b2\u03b1\u2076")
    print("\nWhich evaluates to:")
    
    equation = f"\u0394E\u2082 \u2248 ({final_coeff.p} / {final_coeff.q}) \u00d7 m\u2091c\u00b2\u03b1\u2076"
    
    print(equation)
    
    # Return the coefficient for the final answer block
    return final_coeff

final_coefficient = solve_energy_shift()
# The required final answer format
answer_string = f"<<<{final_coefficient.p}/{final_coefficient.q}>>>"
# print(answer_string)