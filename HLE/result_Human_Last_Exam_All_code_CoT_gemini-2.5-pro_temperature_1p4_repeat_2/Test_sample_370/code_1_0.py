import sympy as sp

def calculate_cross_section():
    """
    Calculates the total cross section for fermion-fermion scattering via pseudoscalar exchange
    in the high-energy limit using symbolic mathematics.
    """
    # Define symbols for the physical quantities
    # g: coupling constant
    # E: energy of one incoming fermion in the CM frame
    # M: mass of the exchanged scalar particle
    g, E, M = sp.symbols('g E M', real=True, positive=True)
    # x is the cosine of the scattering angle theta
    x = sp.Symbol('x')

    # Mandelstam variables in the CM frame (high-energy limit, m_fermion -> 0)
    s = 4 * E**2
    t = -2 * E**2 * (1 - x)
    u = -2 * E**2 * (1 + x)

    # Spin-averaged invariant amplitude squared
    # This is the sum of t-channel, u-channel, and interference terms
    term_t = t**2 / (t - M**2)**2
    term_u = u**2 / (u - M**2)**2
    term_int = t * u / ((t - M**2) * (u - M**2))
    
    M2_avg = g**4 * (term_t + term_u + term_int)

    # Differential cross section for identical fermions
    # The 1/2 factor accounts for identical final states.
    dsig_dOmega = (1 / (128 * sp.pi**2 * s)) * M2_avg

    # The integral over the azimuthal angle phi gives 2*pi
    # The integral over the polar angle theta (d(cos(theta))) is from -1 to 1
    integrand = dsig_dOmega * 2 * sp.pi

    # Perform the symbolic integration over x = cos(theta) from -1 to 1
    total_cross_section = sp.integrate(integrand, (x, -1, 1))

    # Simplify the final expression
    sigma_simplified = sp.simplify(total_cross_section)

    # Print the final result for the total cross section
    # The result shows how sigma depends on E, M, and g.
    print("The total cross section sigma is:")
    # Using sp.pretty_print for a more readable output format for the complex expression.
    sp.pretty_print(sigma_simplified)
    
    # Store the result for the final answer block. The expression is complex,
    # so we will represent it in a string format.
    global final_expression 
    final_expression = str(sigma_simplified)


if __name__ == '__main__':
    final_expression = ""
    calculate_cross_section()
    # The final answer is a complex symbolic expression.
    # It is provided below in the required format.
    # The Python script outputs a pretty-printed version for better readability.
    # print(f"<<<{final_expression}>>>")
