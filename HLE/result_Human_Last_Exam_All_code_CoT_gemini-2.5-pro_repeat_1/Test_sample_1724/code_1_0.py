import sympy

def solve_frequency_correction():
    """
    Calculates the 3rd term of the nonlinear correction to the frequency
    for the Rayleigh-Plesset equation using the Poincare-Lindstedt method.
    """
    # Define the symbolic variable for the polytropic index
    gamma = sympy.Symbol('gamma')

    print("Step 1: Define coefficients from the O(epsilon^2) solution.")
    # From the solution of the second-order equation (for y_2), we get two key coefficients:
    # C0 is the constant term in y_2
    C0 = (3 * gamma) / 4
    # C22 is the amplitude of the cos(2*tau) term in y_2
    C22 = (2 + gamma) / 4
    print(f"Coefficient C0 = {C0}")
    print(f"Coefficient C22 = {C22}\n")

    print("Step 2: Define terms contributing to the frequency correction coefficient sigma_2.")
    # The frequency correction at O(epsilon^2) is found by analyzing the secular terms
    # in the O(epsilon^3) equation. The result for the correction coefficient sigma_2 is
    # derived from combining several terms. Let's calculate sigma_2.
    # The expression for sigma_2 is found from the relation:
    # 2*omega_0*omega_2 = (3*gamma/8) * (-6*gamma**2 + 3*gamma + 2)
    # And omega^2 is expanded as omega_0^2 * (1 + sigma_2 * epsilon^2 + ...)
    # which gives sigma_2 = (2 * omega_0 * omega_2) / omega_0**3 = (2 * omega_0 * omega_2) / (3*gamma*omega_0)
    # My derivation leads to the expression for (2*omega_0*omega_2)/(3*gamma):
    
    # Contribution from quadratic nonlinearities involving y_2
    term_A = -(C0 + C22 / 2) - (3 * gamma + 1) * (C0 - C22 / 2)
    
    # Contribution from the cubic nonlinearity y_1^3
    term_B = (3 * gamma + 1) * (3 * gamma + 2) / 8
    
    # The sum gives the core polynomial of sigma_2, divided by 8.
    sigma_2_poly_numerator = sympy.expand(term_A + term_B)
    
    print(f"The numerator of the sigma_2 coefficient polynomial is: {sigma_2_poly_numerator}\n")
    
    # The full expression for sigma_2 is:
    sigma_2 = sigma_2_poly_numerator / 8

    print("Step 3: Calculate the first nonlinear correction term for omega^2.")
    # The linear frequency squared is omega_0^2 = 3*gamma
    omega0_sq = 3 * gamma
    
    # The first nonlinear correction to omega^2 is delta_omega_sq = omega_0^2 * sigma_2
    delta_omega_sq_poly = sympy.expand(omega0_sq * sigma_2)
    
    print("The polynomial for the first nonlinear correction to omega^2 is:")
    sympy.pprint(delta_omega_sq_poly)
    print("")

    print("Step 4: Identify the third term of this correction polynomial.")
    # The terms are ordered by the power of gamma, from highest to lowest.
    terms = delta_omega_sq_poly.as_ordered_terms()
    
    if len(terms) < 3:
        print("The polynomial has fewer than 3 terms.")
        third_term = 0
    else:
        third_term = terms[2]

    print("The terms of the correction polynomial are:")
    for i, t in enumerate(terms):
        print(f"  Term {i+1}: ", end="")
        sympy.pprint(t)
    print("")

    print("The 3rd term of the nonlinear correction is:")
    sympy.pprint(third_term)
    
    # Extract the numbers from the final term as requested.
    # The term is of the form (n/d)*gamma.
    coeff = third_term / gamma
    n, d = coeff.as_numer_denom()

    print("\nFinal Answer:")
    print(f"The expression for the 3rd term is: ({n}/{d})*gamma")
    print(f"The numbers in the final equation are: {n}, {d}")

solve_frequency_correction()