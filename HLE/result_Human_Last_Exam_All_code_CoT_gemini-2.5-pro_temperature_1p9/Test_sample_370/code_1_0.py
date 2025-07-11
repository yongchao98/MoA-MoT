import sympy as sp

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    mediated by a pseudoscalar boson in the high-energy limit (m=0).
    The result is given in terms of the center-of-mass energy E
    and the scalar boson mass M.
    """
    # Define symbolic variables
    g, E, M, c = sp.symbols('g E M c', real=True, positive=True)
    
    # Mandelstam variables in the CM frame (with m=0)
    s = 4 * E**2
    t = -2 * E**2 * (1 - c)
    u = -2 * E**2 * (1 + c)

    # Spin-averaged squared amplitude, |M|^2
    # M_sq_bar = g^4 * [ u^2/(u-M^2)^2 + t^2/(t-M^2)^2 + t*u/((t-M^2)*(u-M^2)) ]
    term_u = u**2 / (u - M**2)**2
    term_t = t**2 / (t - M**2)**2
    term_tu = (t * u) / ((t - M**2) * (u - M**2))
    
    M_sq_bar = g**4 * (term_t + term_u + term_tu)

    # The differential cross section d(sigma)/d(cos(theta)) is proportional to M_sq_bar
    # dsigma/dc = (1 / (64 * pi * s)) * M_sq_bar
    # Note: The integral over phi gives 2*pi, and there's a symmetry factor 1/2.
    # So the prefactor for the integral over c = cos(theta) is (1/2) * (2*pi) / (64*pi^2*s) = 1 / (64*pi*s)
    
    integrand = M_sq_bar
    
    # Integrate over c = cos(theta) from -1 to 1
    integral_result = sp.integrate(integrand, (c, -1, 1))
    
    # Total cross section sigma
    sigma = (1 / (64 * sp.pi * s)) * integral_result
    
    # Simplify the final expression
    sigma_simplified = sp.simplify(sigma)
    
    # Print the equation
    # Using sp.pretty_print for a more readable formula
    print("The total cross section sigma is:")
    sp.pretty_print(sp.Eq(sp.Symbol('sigma'), sigma_simplified))

    # The final prompt requires outputting the numbers from the equation.
    # Since the result is symbolic, we will pretty print the formula's structure.
    # If the expression is very complex, showing its components might be helpful.
    # Let's extract and show the numerator and denominator
    
    final_fraction = sigma_simplified.as_numer_denom()
    print("\n")
    print("Numerator of the final expression:")
    sp.pretty_print(final_fraction[0])
    print("\n")
    print("Denominator of the final expression:")
    sp.pretty_print(final_fraction[1])
    

calculate_cross_section()
<<<g**4*(2*E**2*M**2 + 4*E**4 - M**4*log(4*E**2/M**2 + 1))/(16*E**2*M**2*pi*(4*E**2 + M**2))>>>