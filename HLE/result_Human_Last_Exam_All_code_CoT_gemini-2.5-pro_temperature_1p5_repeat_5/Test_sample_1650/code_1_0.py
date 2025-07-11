import sympy as sp

def solve_overlap_integral():
    """
    This function calculates the overlap integral for two 2s orbitals in H2+
    using symbolic mathematics and prints the final analytical expression.
    """
    # Step 1: Define symbols
    zeta, R = sp.symbols('zeta R', positive=True, real=True)
    lam, mu = sp.symbols('lam mu', real=True)

    # The problem asks for an expression in R and zeta.
    # It is conventional to define rho = zeta * R for compactness.
    rho = sp.symbols('rho', positive=True, real=True)

    # Step 2: Define the integrand and volume element in elliptical coordinates
    # The combined 2s wave function product (unnormalized) in elliptical coordinates is:
    # (2 - zeta*r_a)*(2 - zeta*r_b) * exp(-zeta*(r_a+r_b)/2)
    # Substituting r_a = R/2*(lam+mu) and r_b = R/2*(lam-mu)
    integrand_part1 = (4 - zeta * R * lam + (zeta**2 * R**2 / 4) * (lam**2 - mu**2))
    integrand_part2 = sp.exp(-zeta * R * lam / 2)
    
    # Volume element part for integration over lam and mu (phi is integrated out)
    # dtau = (R^3/8) * (lam^2 - mu^2) dlam dmu dphi
    # The integration over dphi from 0 to 2*pi yields a factor of 2*pi.
    volume_factor = (2 * sp.pi) * (R**3 / 8) * (lam**2 - mu**2)

    # Full integrand for lambda and mu integration
    full_integrand = integrand_part1 * integrand_part2 * volume_factor
    
    # Step 3: Perform sequential integration
    # Integrate over mu from -1 to 1
    integral_after_mu = sp.integrate(full_integrand, (mu, -1, 1))
    integral_after_mu = sp.simplify(integral_after_mu)

    # Integrate over lambda from 1 to oo (sp.oo represents infinity)
    # The integration can be slow, so assumptions on zeta and R are important.
    # The 'noconds=False' can sometimes help with complex integrals.
    integral_after_lam = sp.integrate(integral_after_mu, (lam, 1, sp.oo), noconds=False)
    integral_after_lam = sp.simplify(integral_after_lam)

    # Step 4: Apply the normalization constant
    # The normalization constant N for a 2s orbital is (zeta^3 / (32*pi))^(1/2).
    # The overlap integral has N^2.
    N_squared = zeta**3 / (32 * sp.pi)
    
    S = N_squared * integral_after_lam
    # Simplify the final expression and substitute R with rho/zeta for a compact form.
    S_final = S.subs(R, rho / zeta)
    S_final = sp.simplify(S_final)
    
    # Step 5: Format and print the final equation
    # The expression should be of the form: exp(-rho/2) * (polynomial in rho)
    prefactor = sp.exp(-rho / 2)
    poly_part = sp.collect(sp.expand(S_final / prefactor), rho)

    # Build the output string for the polynomial part
    poly_in_rho = sp.Poly(poly_part, rho)
    terms_dict = poly_in_rho.as_dict()
    sorted_powers = sorted(terms_dict.keys(), key=lambda p: p[0], reverse=True)

    poly_str_terms = []
    for power_tuple in sorted_powers:
        power = power_tuple[0]
        coeff = terms_dict[power_tuple]
        
        if coeff == 0:
            continue
            
        sign = "+ " if coeff > 0 else "- "
        abs_coeff = abs(coeff)

        # First term should not have a leading sign unless negative
        if not poly_str_terms and sign == "+ ":
            sign = ""
        elif poly_str_terms:
            abs_coeff_str = f"({abs_coeff})"
        else:
            abs_coeff_str = f"{abs_coeff}"

        if power == 0:
            term_str = f"{sign}{abs_coeff}"
        elif power == 1:
            term_str = f"{sign}{abs_coeff_str}*(ζ*R)"
        else:
            term_str = f"{sign}{abs_coeff_str}*(ζ*R)**{power}"
        poly_str_terms.append(term_str)
        
    final_poly_str = " ".join(poly_str_terms)

    print("The analytical expression for the overlap integral S(R, ζ) is:")
    print(f"S(R, ζ) = exp(-ζ*R/2) * ( {final_poly_str} )")

# Run the calculation and print the result
solve_overlap_integral()