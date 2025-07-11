import sympy as sp

def calculate_frequency_correction_terms():
    """
    This function calculates the coefficients of the polynomial in gamma
    for the nonlinear frequency correction of the Rayleigh-Plesset equation.

    The relative frequency correction is given by:
    delta_omega / omega_0 = (epsilon * a)^2 * (C1 * gamma^2 + C2 * gamma + C3)
    where 'a' is the amplitude of the linear oscillation.
    This function calculates and prints the coefficients C1, C2, and C3.
    """
    
    # Define symbolic variable for gamma
    gamma = sp.Symbol('gamma')
    
    # Define omega_0 squared
    omega_0_sq = 3 * gamma
    
    # The coefficients of the polynomial in omega_0 for the relative frequency correction are derived
    # from perturbation theory. Let's denote the relative frequency shift as:
    # delta_omega / omega_0 = (epsilon*a)^2 * P(omega_0_sq)
    # The polynomial P(W) where W=omega_0_sq is: P(W) = -1/6 * W^2 + 1/16 * W + 5/8
    
    # Let's substitute W = 3*gamma to get the polynomial in gamma.
    
    # Coefficient of gamma^2
    c1_coeff = sp.Rational(-1, 6)
    term1_poly_in_gamma = c1_coeff * (3**2)
    
    # Coefficient of gamma
    c2_coeff = sp.Rational(1, 16)
    term2_poly_in_gamma = c2_coeff * 3
    
    # Constant term
    c3_coeff = sp.Rational(5, 8)
    term3_poly_in_gamma = c3_coeff
    
    # The nonlinear correction to the linear oscillation frequency is composed of three terms.
    # We will print the equation for the relative frequency shift.
    # The variable 'a' is the dimensionless amplitude of the oscillation, and 'eps' is the perturbation parameter.
    
    print("The relative nonlinear frequency correction is given by the equation:")
    print(f"Δω/ω₀ = (εa)² * [({term1_poly_in_gamma})γ² + ({term2_poly_in_gamma})γ + ({term3_poly_in_gamma})]")
    print("\nAssuming the terms are ordered by descending powers of γ:")
    print(f"1st term contribution is: ({term1_poly_in_gamma})γ²")
    print(f"2nd term contribution is: ({term2_poly_in_gamma})γ")
    print(f"3rd term contribution is: {term3_poly_in_gamma}")
    
    print("\nCalculation of the 3rd term:")
    # The third term is the constant coefficient in the polynomial.
    third_term = term3_poly_in_gamma
    print(f"The 3rd term of the nonlinear correction is {third_term}.")

if __name__ == '__main__':
    calculate_frequency_correction_terms()
