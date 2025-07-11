def solve_polymer_force():
    """
    This function prints the derived force law for a thermally isolated polymer.
    """

    # Define the constants and variables as strings for the output equation
    # The numbers in the equation are 3 and 2.
    num_3 = 3
    num_2 = 2

    # Explain the variables used in the final equation
    print("The force law for a thermally isolated polymer is F(x), which represents the magnitude of the attractive force between the polymer ends.")
    print("\nThe variables are defined as:")
    print("  F(x): The magnitude of the force at extension x.")
    print("  x:   The separation distance between the polymer ends.")
    print("  l:   The length of a single segment of the polymer.")
    print("  n:   The number of segments in the polymer.")
    print("  E(0): The kinetic energy of the polymer when the extension x is zero.")
    
    # Construct and print the final equation as a formatted string
    equation_str = f"\nF(x) = ({num_3} * E(0) * x / (n**2 * l**2)) * exp(-{num_3} * x**2 / ({num_2} * n**2 * l**2))"
    
    print("\nThe derived force law is:")
    print(equation_str)

solve_polymer_force()