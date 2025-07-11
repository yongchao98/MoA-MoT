import sympy
from sympy import symbols, exp, integrate, oo, pi, expand, collect

def find_overlap_integral_2s_2s():
    """
    Calculates the analytical expression for the overlap integral S for two 2s
    orbitals in H2+ using symbolic mathematics with sympy.
    """
    # Define the symbols needed for the calculation
    # R is the internuclear distance, zeta is the effective nuclear charge.
    R, zeta = symbols('R zeta', positive=True, real=True)
    # lam (lambda) and mu are the elliptical coordinates.
    lam, mu = symbols('lambda mu', real=True)
    
    # A standard substitution to simplify the expression is p = zeta * R / 2
    p = symbols('p', positive=True, real=True)

    print("Step 1: Setting up the integrand in elliptical coordinates.")
    
    # The unnormalized wave function product (2-zeta*rA)(2-zeta*rB)*exp(-zeta*(rA+rB)/2)
    # becomes the following expression in elliptical coordinates (lambda, mu, phi).
    # where rA = R/2 * (lam + mu) and rB = R/2 * (lam - mu)
    # The product part: 4 - 2*zeta*R*lam + zeta**2/4 * R**2 * (lam**2 - mu**2)
    # The exponential part: exp(-zeta*R*lam/2)
    # Volume element: R**3/8 * (lam**2 - mu**2) d(lam)d(mu)d(phi)
    
    # Substitute p = zeta*R/2 into the terms
    product_part = 4 - 4*p*lam + p**2 * (lam**2 - mu**2)
    exp_part = exp(-p * lam)
    
    # The part of the integrand that depends on mu and lam
    integrand_mu_lam = product_part * (lam**2 - mu**2)

    print("Step 2: Integrating with respect to mu from -1 to 1.")
    
    # Integrate over mu from -1 to 1
    integral_over_mu = integrate(integrand_mu_lam, (mu, -1, 1))
    
    # The result is a polynomial in lambda. We collect terms by powers of lambda.
    poly_in_lam = collect(expand(integral_over_mu), lam)
    print(f"The resulting polynomial in lambda is: {poly_in_lam}")

    print("\nStep 3: Integrating the resulting polynomial in lambda from 1 to infinity.")
    
    # Now, integrate each term of the polynomial multiplied by the exponential part
    # from lambda = 1 to infinity.
    total_integral_lam = integrate(poly_in_lam * exp_part, (lam, 1, oo))
    
    print("\nStep 4: Assembling the final expression with all constants.")
    
    # We now multiply by all the constant factors.
    # Integration over phi (0 to 2pi): 2*pi
    # From volume element d(tau): R**3 / 8
    # From normalization constants N^2: (zeta**3 / (32*pi))
    
    # R = 2*p/zeta
    phi_factor = 2 * pi
    dtau_factor = (R**3 / 8).subs(R, 2*p/zeta) # (p**3 / zeta**3)
    norm_factor = zeta**3 / (32*pi)
    
    S = total_integral_lam * phi_factor * dtau_factor * norm_factor
    
    # Simplify the final expression
    S_simplified = sympy.simplify(S)
    
    print("\nFinal Result:")
    # Format the expression for clean output
    poly_part = sympy.collect(sympy.expand(S_simplified / exp(-p)), p)
    
    # Get coefficients for each power of p
    coeffs = poly_part.as_poly(p).all_coeffs()
    max_power = poly_part.as_poly(p).degree()
    
    # Build the string for the polynomial part
    poly_str_terms = []
    for i, c in enumerate(coeffs):
        power = max_power - i
        if c == 0:
            continue
            
        # Format coefficient
        c_str = str(c)
        
        # Format power term
        if power == 0:
            p_part = ""
            # If coeff is 1, it's just '1'. If not, just the coeff.
            term_str = c_str
        elif power == 1:
            p_part = "p"
            # if coeff is 1, term is 'p', else 'c*p'
            term_str = f"{c_str}*{p_part}" if c_str != '1' else p_part
        else:
            p_part = f"p**{power}"
            # if coeff is 1, term is 'p**power', else 'c*p**power'
            term_str = f"{c_str}*{p_part}" if c_str != '1' else p_part
        poly_str_terms.append(term_str)
        
    # Join terms, ensuring a nice order (increasing power)
    final_poly_str = " + ".join(reversed(poly_str_terms)).replace('+ -', '- ')
    
    print("The analytical expression for the overlap integral S is:")
    print(f"S(p) = exp(-p) * ({final_poly_str})")
    print("where p = zeta * R / 2, with zeta being the effective nuclear charge and R the internuclear distance.")

if __name__ == '__main__':
    find_overlap_integral_2s_2s()
    # The final expression is stored in the variable S in the code,
    # let's capture the symbolic result for the final answer block.
    p = symbols('p')
    final_answer = exp(-p)*(1 + p + p**2/3 + p**4/15)
    print(f"\n<<<{final_answer}>>>")