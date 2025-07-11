def print_polymer_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed chain polymer.

    The force F(x) is the force of attraction between the polymer ends.
    Variables:
    x: The separation of the polymer ends.
    l: The length of a single strut (segment).
    n: The number of segments in the polymer.
    E(0): The kinetic energy of the polymer at zero extension (x=0).
    """

    # The derived formula for the force law
    formula = "F(x) = (2 * E(0) * x / (n**2 * l**2)) * exp(x**2 / (n**2 * l**2))"

    # The numerical constants appearing in the final equation
    coefficient_in_prefactor = 2
    exponent_of_n_and_l = 2
    
    print("The force law between the polymer ends for a thermally isolated chain is:")
    print(formula)
    
    print("\nTo satisfy the request to output each number in the final equation, here are the numerical constants from the formula:")
    print(f"- The numerical coefficient in the term outside the exponential is {coefficient_in_prefactor}.")
    print(f"- The exponent for both 'n' and 'l' in the denominator and in the exponential is {exponent_of_n_and_l}.")

# Execute the function to display the result
print_polymer_force_law()
