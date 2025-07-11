import sympy as sp

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering
    mediated by a pseudoscalar boson in the high-energy limit (m=0).
    """
    # Define symbolic variables
    # g: coupling constant
    # E: center-of-mass energy of one fermion
    # M: mass of the exchanged scalar boson
    # x: cosine of the scattering angle theta
    g, E, M = sp.symbols('g E M', real=True, positive=True)
    x = sp.symbols('x', real=True)

    # Mandelstam variables in the CM frame for m=0
    # s = (2*E)**2
    t = -2 * E**2 * (1 - x)
    u = -2 * E**2 * (1 + x)

    # Spin-averaged squared matrix element <|M|^2>
    # Note: M_sq is the squared matrix element, M is the boson mass.
    M2_sq_term = (t**2 / (t - M**2)**2) + \
                 (u**2 / (u - M**2)**2) - \
                 (t * u / ((t - M**2) * (u - M**2)))
    
    avg_M2_sq = g**4 * M2_sq_term

    # The total cross section sigma is the integral of the differential cross section
    # d(sigma)/d(Omega) over the solid angle Omega.
    # d(sigma)/d(Omega) = (1 / (64 * pi^2 * s)) * <|M|^2>
    # s = 4*E^2
    # sigma = integral(d(sigma)/d(Omega) * 2*pi d(cos(theta))) from -1 to 1
    # sigma = (1 / (128 * pi * E^2)) * integral(<|M|^2> d(x))
    
    # Perform the integration over x = cos(theta) from -1 to 1
    integral_val = sp.integrate(avg_M2_sq, (x, -1, 1))
    
    # Calculate the total cross section
    sigma = integral_val / (128 * sp.pi * E**2)
    
    # Simplify the final expression
    simplified_sigma = sp.simplify(sigma)
    
    # Print the final result
    print("The total cross section sigma is:")
    # Using sp.pretty_print for better formatting of the equation
    sp.pprint(simplified_sigma, use_unicode=True)

if __name__ == '__main__':
    calculate_cross_section()