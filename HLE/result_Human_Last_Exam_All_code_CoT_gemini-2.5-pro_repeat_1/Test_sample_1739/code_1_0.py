import sympy

def calculate_frequency_correction_coefficient():
    """
    This function calculates the coefficient K that determines the second-order
    frequency correction for the given Rayleigh-Plesset equation.
    It follows the perturbation method analysis.
    """
    # Define gamma as a symbol
    gamma = sympy.Symbol('gamma')

    # Linear oscillation frequency squared
    omega0_sq = 3 * gamma

    # --- Step 1: Coefficients from O(epsilon) solution ---
    # The first order solution x1 contains terms with coefficients C0 and C2.
    # These arise from the quadratic nonlinearities.
    # C0 is the coefficient of the constant term (A*A_bar).
    # C2 is the coefficient of the second harmonic term (A^2).
    
    # Coeff of A*A_bar in the forcing term at O(epsilon)
    coeff_A_Abar = -omega0_sq + 3*gamma*(3*gamma + 1)
    # C0 = coeff / omega0_sq
    C0 = sympy.simplify(coeff_A_Abar / omega0_sq)

    # Coeff of A^2 in the forcing term at O(epsilon)
    coeff_A2 = (5/2)*omega0_sq + (3*gamma*(3*gamma+1))/2
    # C2 = coeff / (omega0_sq - (2*omega0)^2) = coeff / (-3*omega0_sq)
    C2 = sympy.simplify(coeff_A2 / (-3 * omega0_sq))

    # --- Step 2: Calculate the factor K ---
    # K is the sum of contributions from quadratic interactions (K_int)
    # and the cubic term (K_G).

    # Contribution from interaction of quadratic terms
    # K_int = (omega0_sq + 3*gamma*(3*gamma+1))*(C0 + C2) - 2*omega0_sq*C2
    term1_K_int = (omega0_sq + 3*gamma*(3*gamma+1)) * (C0 + C2)
    term2_K_int = -2 * omega0_sq * C2
    K_int = term1_K_int + term2_K_int

    # Contribution from the original cubic term in the expansion
    # The coefficient of x^3 is -gamma*(3*gamma+1)*(3*gamma+2)/2.
    # The contribution to K is 3 times this coefficient.
    K_G = 3 * (-gamma * (3*gamma + 1) * (3*gamma + 2) / 2)

    # Total K
    K = sympy.simplify(K_int + K_G)

    # --- Step 3: Extract coefficients of the polynomial K ---
    # Expand K into a polynomial in gamma
    K_poly = sympy.Poly(K, gamma)
    
    print("The factor K, which determines the frequency correction, is a polynomial in gamma:")
    print(f"K = {K_poly.as_expr()}")
    
    # Get the coefficients in descending order of power
    coeffs = K_poly.coeffs()
    
    print("\nThe coefficients of the polynomial K (from highest power to lowest) are:")
    for i, c in enumerate(coeffs):
        print(f"Coefficient of gamma^{K_poly.degree()-i}: {float(c)}")

    # The problem asks for the 3rd term's coefficient
    third_term_coeff = coeffs[2]
    
    return float(third_term_coeff)

# Run the calculation and get the final answer
final_answer = calculate_frequency_correction_coefficient()
print(f"\nThe 3rd term of the nonlinear correction is the coefficient of the third term in the polynomial K, which is: {final_answer}")
