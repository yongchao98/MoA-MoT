def display_conductance_formula():
    """
    This function prints the derived formula for the four-terminal conductance G_12,34.
    """
    # Symbolic representation of the parameters
    M = "M"  # Total number of spin-degenerate edge states
    N = "N"  # Number of reflected edge states
    G0 = "e^2/h" # Quantum of conductance

    print("The derived four-terminal conductance G_12,34 is given by the formula:")
    
    # Using placeholders to format the final equation string
    # Numerator has two parts: M and (M-N)
    # Denominator is N
    numerator_part1 = M
    numerator_part2 = f"({M} - {N})"
    denominator = N
    
    print(f"\nG_12,34 = ( {numerator_part1} * {numerator_part2} ) / {denominator} * ({G0})")
    
    print("\nWhere the components of the formula are:")
    print(f"- {M}: Total number of spin-degenerate edge states.")
    print(f"- {N}: Number of edge states reflected by the Quantum Point Contact.")
    print(f"- ({G0}): The quantum of conductance.")

if __name__ == '__main__':
    display_conductance_formula()