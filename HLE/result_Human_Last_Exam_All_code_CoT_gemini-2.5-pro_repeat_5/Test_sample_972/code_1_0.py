import sympy

def calculate_amplitude_in_slab():
    """
    This function calculates and prints the symbolic expression for the amplitude
    of an electromagnetic wave at the exit of a slab with time-varying properties.
    """

    # Step 1: Define the symbols for the physical quantities.
    # A: Amplitude of the incident wave
    # L: Length of the slab in the x-direction
    # alpha: The rate of change of the slab's electromagnetic properties
    # c: The speed of light in vacuum
    A = sympy.Symbol('A')
    L = sympy.Symbol('L')
    alpha = sympy.Symbol('alpha')
    c = sympy.Symbol('c')

    # Step 2: State the key physical insights of the problem.
    # The problem states epsilon(t)/epsilon_0 = mu(t)/mu_0 = alpha*t + beta.
    # From this, we can derive two important properties of the slab:
    # a) The refractive index: n(t) = sqrt(epsilon*mu / (epsilon_0*mu_0)) = alpha*t + beta.
    # b) The wave impedance: Z(t) = sqrt(mu(t)/epsilon(t)) = sqrt(mu_0/epsilon_0) = Z_0.
    # Since the slab's impedance Z(t) is always equal to the vacuum impedance Z_0,
    # there are no reflections at the boundaries. The wave enters and exits the slab completely.

    # Step 3: Determine how the amplitude evolves inside the slab.
    # For a wave propagating in a medium with time-varying refractive index n(t) and
    # matched impedance, the electric field amplitude |E| is proportional to 1/sqrt(n(t)).
    # Therefore, the ratio of the output amplitude (A_out) to the input amplitude (A_in)
    # for a wavefront is given by: A_out / A_in = sqrt(n_in / n_out), where n_in and n_out
    # are the refractive indices at the time of entry and exit, respectively.

    # Step 4: Relate the refractive index ratio to the slab length L.
    # A wavefront travels at speed v(t) = c / n(t). The time it takes to cross the
    # slab of length L is found by integrating: L = integral(v(t) dt).
    # L = integral from t_in to t_out of (c / (alpha*t + beta)) dt
    # This integration yields: L = (c/alpha) * ln(n_out / n_in).
    # Rearranging for the ratio gives: n_out / n_in = exp(alpha * L / c).
    # So, the required ratio is: n_in / n_out = exp(-alpha * L / c).

    # Step 5: Combine the results to find the final amplitude.
    # Substitute the ratio from Step 4 into the amplitude relation from Step 3:
    # A_out / A_in = sqrt(exp(-alpha * L / c))
    # A_out / A_in = exp(-alpha * L / (2 * c))
    # The incident amplitude A_in is given as A.
    
    # We construct the final expression for the output amplitude A_out.
    # The exponent in the expression is -alpha*L/(2*c)
    exponent = -alpha * L / (2 * c)
    
    # The output amplitude is A * exp(exponent)
    A_out_expression = A * sympy.exp(exponent)

    # Step 6: Print the final answer in the required format.
    # The final equation gives the amplitude at x=L. We explicitly show all
    # variables and numbers in the equation.
    print("The amplitude of the electric field at the rightmost boundary of the slab is A_out, where:")
    # We use string formatting to build the equation clearly.
    # The number 2 from the derivation is explicitly shown. The negative sign represents the number -1.
    final_equation_str = f"A_out = {A} * exp(-({alpha} * {L}) / (2 * {c}))"
    print(final_equation_str)

# Execute the function to get the result.
calculate_amplitude_in_slab()