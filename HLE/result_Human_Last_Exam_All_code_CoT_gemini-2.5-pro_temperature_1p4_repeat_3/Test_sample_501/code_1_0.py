def polymer_force_law():
    """
    This function prints the derived force law for a thermally isolated
    freely jointed polymer chain.
    """
    
    # Define the variables for the equation symbolically
    F = "F"
    E0 = "E(0)"
    x = "x"
    n = "n"
    l = "l"
    
    # Print the final force law equation
    print("The force law for the internal attraction between the polymer ends is:")
    # The numbers in the equation (2, 2) are printed as required.
    print(f"{F} = - (2 * {E0} * {x}) / ({n}^2 * {l}^2)")

    print("\nWhere:")
    print(f"  {F} is the internal attractive force between the polymer ends.")
    print(f"  {E0} is the kinetic energy of the polymer at zero extension (x=0).")
    print(f"  {x} is the separation of the ends (assumed to be small).")
    print(f"  {n} is the number of segments in the polymer chain (assumed to be large).")
    print(f"  {l} is the length of a single segment.")

if __name__ == '__main__':
    polymer_force_law()