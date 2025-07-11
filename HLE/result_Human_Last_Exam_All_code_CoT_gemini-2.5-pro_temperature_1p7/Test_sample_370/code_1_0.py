import sympy

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    mediated by a pseudoscalar particle in the high-energy limit (m=0).
    The Lagrangian for the interaction is given by:
    L = 1/2 (d(phi))^2 - 1/2 M^2 phi^2 + psi_bar(i*gamma*d - m)psi - g*psi_bar*gamma_5*psi*phi
    """
    # Define symbolic variables for the parameters
    g, E, M = sympy.symbols('g E M', real=True, positive=True)
    # Define the integration variable x = cos(theta)
    x = sympy.symbols('x')

    # In the Center-of-Mass (CM) frame, for massless fermions (m=0):
    # s is the square of the total CM energy.
    # E is the energy of a single incoming fermion.
    s = 4 * E**2
    
    # Mandelstam variables t and u in terms of E and x = cos(theta)
    # where theta is the scattering angle.
    t = -2 * E**2 * (1 - x)
    u = -2 * E**2 * (1 + x)

    # The spin-averaged invariant amplitude squared |M_bar|^2
    # This result is derived from the t-channel and u-channel Feynman diagrams,
    # including the interference term. For massless fermions interacting via a
    # pseudoscalar exchange, it is given by:
    # |M_bar|^2 = g^4 * [ t^2/(t-M^2)^2 + u^2/(u-M^2)^2 + t*u/((t-M^2)(u-M^2)) ]
    
    # Let's write the expression for the integrand
    term_t = t**2 / (t - M**2)**2
    term_u = u**2 / (u - M**2)**2
    term_interference = (t * u) / ((t - M**2) * (u - M**2))
    
    M_bar_squared = g**4 * (term_t + term_u + term_interference)

    # The total cross-section sigma is obtained by integrating the differential cross-section.
    # For two identical fermions in the final state, the formula is:
    # sigma = (1 / (256 * pi * E^2)) * Integral(|M_bar|^2, dx) from -1 to 1.
    # The pre-factor includes the flux factor, phase space, and a symmetry factor of 1/2.
    
    # Perform the definite integral of |M_bar|^2 with respect to x from -1 to 1
    integral_M_squared = sympy.integrate(M_bar_squared, (x, -1, 1))

    # Calculate the total cross section sigma
    sigma = (1 / (256 * sympy.pi * E**2)) * integral_M_squared

    # Simplify the final expression
    sigma_simplified = sympy.simplify(sigma)
    
    # Print the final result for the total cross section
    print("The total cross section for psi-psi scattering in the high-energy limit (m=0) is:")
    sympy.pprint(sigma_simplified)

    # The user asked for the final formula, so we return it in a string format as well.
    final_formula_str = str(sigma_simplified)
    # To conform to the output format, we extract the core expression
    final_formula_str = final_formula_str.replace("**", "^")
    return final_formula_str


if __name__ == '__main__':
    final_answer = calculate_cross_section()
    # The requested output format is <<<answer content>>>.
    # The full symbolic answer is very long, so we print it above and
    # provide a notification that the calculation is complete and the formula is printed.
    print(f"\nFinal formula has been computed and printed above.")
    print(f"<<<{final_answer}>>>")
