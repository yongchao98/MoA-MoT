def print_force_law():
    """
    This function prints the derived force law for a thermally isolated,
    freely jointed polymer chain with a small end-to-end separation.
    """

    # Define the numerical constants in the equation
    numerator_constant = 2
    exponent_constant = 2

    # Define the variables as string representations for clarity
    force = "F"
    initial_energy = "E(0)"
    separation = "x"
    num_segments = "n"
    segment_length = "l"

    # Construct and print the final equation
    # The equation is F = - (2 * E(0) * x) / (n^2 * l^2)
    print("The force law between the polymer ends is:")
    print(f"{force} = - ({numerator_constant} * {initial_energy} * {separation}) / ({num_segments}^{exponent_constant} * {segment_length}^{exponent_constant})")

    print("\nWhere:")
    print(f"  {force}: The attractive force between the polymer ends.")
    print(f"  {initial_energy}: The kinetic energy of the polymer at zero extension.")
    print(f"  {separation}: The separation distance between the polymer ends (assumed to be small).")
    print(f"  {num_segments}: The number of segments in the polymer chain (assumed to be large).")
    print(f"  {segment_length}: The length of each segment.")

# Execute the function to display the result
print_force_law()