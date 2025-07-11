def print_force_law_formula():
    """
    Prints the derived formula for the force between the ends of a thermally isolated,
    freely jointed polymer chain.

    The formula is expressed in terms of:
    - F: The force of attraction between the polymer ends.
    - E(0): The kinetic energy of the polymer at zero extension (x=0).
    - x: The separation distance between the polymer ends.
    - n: The number of segments in the polymer chain.
    - l: The length of each segment.
    """

    # Define the symbols used in the equation for clarity
    force_symbol = "F(x)"
    energy_at_zero_ext_symbol = "E(0)"
    extension_symbol = "x"
    num_segments_symbol = "n"
    segment_length_symbol = "l"

    # Construct the string representing the formula
    # Note: '**' denotes exponentiation (e.g., n**2 is n-squared)
    # The exponential part: exp(x**2 / (n**2 * l**2))
    # The linear part: - (2 * E(0) * x) / (n**2 * l**2)
    
    # We explicitly include the numbers from the derivation as requested
    term1_numerator = f"2 * {energy_at_zero_ext_symbol} * {extension_symbol}"
    term1_denominator = f"({num_segments_symbol}**2 * {segment_length_symbol}**2)"
    term1 = f"({term1_numerator} / {term1_denominator})"
    
    exp_numerator = f"{extension_symbol}**2"
    exp_denominator = f"({num_segments_symbol}**2 * {segment_length_symbol}**2)"
    exp_term = f"exp({exp_numerator} / {exp_denominator})"

    final_equation = f"{force_symbol} = - {term1} * {exp_term}"

    print("The derived force law between the polymer ends is:")
    print(final_equation)

if __name__ == "__main__":
    print_force_law_formula()