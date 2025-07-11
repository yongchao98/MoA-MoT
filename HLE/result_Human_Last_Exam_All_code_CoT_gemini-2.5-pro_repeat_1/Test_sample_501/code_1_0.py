import math

def display_force_law():
    """
    This function derives and displays the force law for a thermally isolated
    freely jointed polymer chain.
    
    The derivation is outlined in the plan above. The final result is presented
    as a formatted string.
    """
    
    # The derived force law for large n is:
    # F(x) = - (3 * E(0) * x / (n^2 * l^2)) * exp(3 * x^2 / (2 * n^2 * l^2))
    # where:
    # F(x) = Force of attraction at extension x
    # E(0) = Kinetic energy at zero extension
    # x = Separation of the polymer ends
    # n = Number of identical mass points (large)
    # l = Length of each strut

    print("The derived force law for a thermally isolated polymer with a large number of links 'n' is:")
    
    # Printing the final equation with each number explicitly shown, as requested.
    final_equation = (
        "F(x) = - ( (3 * E(0) * x) / (n^2 * l^2) ) * exp( (3 * x^2) / (2 * n^2 * l^2) )"
    )
    
    print("\n" + final_equation + "\n")
    print("Where:")
    print("  F(x): The attractive force between the polymer ends.")
    print("  E(0): The kinetic energy of the polymer at zero extension (x=0).")
    print("  x:    The separation distance between the polymer ends.")
    print("  n:    The number of mass points in the chain.")
    print("  l:    The length of a single link between mass points.")
    print("  exp(): The exponential function.")

# Execute the function to display the result.
display_force_law()
