def print_partition_function():
    """
    This function prints the derived partition function Z in a formatted way.
    The equation includes numerical constants and symbolic representations of
    physical quantities.
    """
    # Define the integer constants that appear in the final equation.
    const_1 = 1
    const_2 = 2

    # Define unicode characters for the physical symbols for better readability.
    beta_symbol = "\u03B2"  # β for inverse temperature
    mu_symbol = "\u03BC"    # μ for chemical potential

    # Print the final equation for the partition function Z.
    print("The partition function Z for a system with Hamiltonian H = -μN is:")
    print(f"Z = {const_1} / ({const_1} - exp({const_2} * {beta_symbol} * {mu_symbol}))")
    
    # Add the condition for the validity of the result.
    print(f"\nNote: This result is valid for {mu_symbol} < 0.")

if __name__ == "__main__":
    print_partition_function()