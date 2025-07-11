def print_force_law():
    """
    This function prints the derived force law for a thermally isolated polymer chain.
    """
    # Define the variables for clarity in the final expression
    force = "F(x)"
    initial_energy = "E(0)"
    extension = "x"
    num_segments = "n"
    segment_length = "l"
    
    # The derived force law for small extension x is F(x) = - (2 * E(0) * x) / (n^2 * l^2)
    # The code below will print this equation, showing all the components.
    
    print("The derived force law between the polymer ends is:")
    
    # We construct and print the equation string.
    # The numbers in the equation, like '2', are explicitly part of the string.
    equation = f"{force} = - (2 * {initial_energy} * {extension}) / ({num_segments}**2 * {segment_length}**2)"
    
    print(equation)

# Execute the function to display the result
print_force_law()
