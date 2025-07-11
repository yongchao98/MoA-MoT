import sympy

def print_bare_greens_function_expression(omega, epsilon_k, particle_type='particle'):
    """
    Prints the expression for the Feynman bare Green's function G_0(k, omega).

    Args:
        omega (float or int): The frequency of the particle.
        epsilon_k (float or int): The single-particle energy eigenvalue for state k.
        particle_type (str): Either 'particle' (for states above the Fermi level)
                             or 'hole' (for states below the Fermi level).
    """
    # Use sympy for nice symbolic printing of eta
    eta = sympy.Symbol('η')
    i = sympy.I

    # The sign of the i*eta term depends on whether we are looking at a particle or a hole
    if particle_type == 'particle':
        sign = 1
        sign_str = '+'
    elif particle_type == 'hole':
        sign = -1
        sign_str = '-'
    else:
        raise ValueError("particle_type must be 'particle' or 'hole'")

    # Construct the symbolic expression for the denominator
    denominator = omega - epsilon_k + sign * i * eta

    # Print the explanation and the final equation
    print("The functional dependence of the bare Green's function G_0 on the single-particle energy ε_k is given by:")
    print(f"\nG_0(k, ω) = 1 / (ω - ε_k {sign_str} iη)\n")
    print("This shows G_0 has a simple pole at the energy eigenvalue ε_k.")
    print("\n--- Example Calculation ---")
    print(f"For a frequency ω = {omega} and a single-particle energy ε_k = {epsilon_k}:")
    print("\nThe equation is:")
    
    # Use sympy.pretty_print for a better looking formula
    # But to follow the prompt's instruction to output each number, an f-string is better.
    # We will build the string manually to show the substitution.
    
    final_equation = f"G_0 = 1 / ({omega} - {epsilon_k} {sign_str} iη)"
    print(final_equation)
    
    # Also show the result after subtracting the real numbers
    real_part_denom = omega - epsilon_k
    final_equation_simplified = f"G_0 = 1 / ({real_part_denom} {sign_str} iη)"
    print("\nWhich simplifies to:")
    print(final_equation_simplified)


# --- User-defined parameters ---
# Let's consider a particle (e.g., an electron added to a system)
example_omega = 5.0  # Example frequency in some energy units
example_epsilon_k = 3.5  # Example single-particle energy in the same units

# Run the function to print the result
print_bare_greens_function_expression(example_omega, example_epsilon_k, particle_type='particle')