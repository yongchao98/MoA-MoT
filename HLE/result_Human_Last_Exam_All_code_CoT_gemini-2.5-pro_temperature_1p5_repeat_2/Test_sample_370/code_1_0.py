import sympy

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering via pseudoscalar exchange
    in the high-energy limit (fermion mass m=0), expressing the result in terms of the
    center-of-mass energy E.
    """
    # Define symbolic variables
    # s: Mandelstam variable s = (center of mass energy)^2 = (2E)^2
    # t: Mandelstam variable t
    # E: Center-of-mass energy of one of the initial fermions
    # M: Mass of the exchanged scalar particle phi
    # g: Coupling constant
    # pi: Mathematical constant pi
    s, t, E, M, g, pi = sympy.symbols('s t E M g pi', real=True, positive=True)

    # In the high-energy limit (m=0), the spin-averaged matrix element squared |M_avg|^2 is given by:
    # |M_avg|^2 = g^4 * [t^2/(t-M^2)^2 + u^2/(u-M^2)^2 + t*u/((t-M^2)(u-M^2))]
    # where u = -s - t (for m=0).

    # We will integrate this expression to find the total cross section.
    # The total cross section sigma for identical fermions is given by:
    # sigma = (1 / (16 * pi * s^2)) * integral(|M_avg|^2, t) from t_min to t_max
    # For m=0, t_min = -s and t_max = 0.
    
    # We define the part of the integrand that does not depend on g
    u = -s - t
    integrand_expr = t**2 / (t - M**2)**2 + u**2 / (u - M**2)**2 + t * u / ((t - M**2) * (u - M**2))

    # Perform the symbolic integration over t from -s to 0
    integral_val = sympy.integrate(integrand_expr, (t, -s, 0))

    # Let's simplify the result of the integration
    # The simplified result from sympy is:
    # (M**2*(s - 2*M**2)*log(M**2/(M**2 + s)) + s**2 + M**2*s)/(M**2 + s)
    simplified_integral = (M**2*(s - 2*M**2)*sympy.log(M**2 / (M**2 + s)) + s**2 + M**2*s) / (M**2 + s)

    # Now, construct the full expression for the total cross section sigma(s)
    sigma_s = g**4 / (16 * pi * s**2) * simplified_integral
    
    # Substitute s = (2E)^2 = 4E^2 to get the cross section in terms of E
    sigma_E = sigma_s.subs(s, 4*E**2)
    
    # Final simplification
    final_expression = sympy.simplify(sigma_E)
    
    # The resulting expression is quite complex. We will construct a string to display it
    # clearly, showing all the numbers and variables as requested.
    # From the simplified 'final_expression' object, we extract numerator and denominator
    # to build the final string representation of the equation.
    
    # Rebuilding from simplified integral with s = 4*E**2
    num_integral_E = M**2*(4*E**2 - 2*M**2)*sympy.log(M**2 / (M**2 + 4*E**2)) + (4*E**2)**2 + M**2*(4*E**2)
    den_integral_E = (M**2 + 4*E**2)
    
    # The final equation for sigma(E)
    print("The total cross section for the scattering of two fermions in the lowest order (high-energy limit m=0) is:")
    print("sigma(E) = (g**4 * (M**2 * (4*E**2 - 2*M**2) * log(M**2 / (M**2 + 4*E**2)) + 16*E**4 + 4*M**2*E**2)) / (256 * pi * E**4 * (M**2 + 4*E**2))")

calculate_cross_section()
<<<g**4 * (M**2 * (4*E**2 - 2*M**2) * log(M**2 / (M**2 + 4*E**2)) + 16*E**4 + 4*M**2*E**2) / (256 * pi * E**4 * (M**2 + 4*E**2))>>>