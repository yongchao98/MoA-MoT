def solve_joule_heat_dissipation():
    """
    This function calculates and displays the formula for the total Joule heat
    dissipated by a leaking and shrinking charged sphere.
    """
    
    # The final derived formula for the total Joule heat (H) is:
    # H = 2 * π * ε₀ * a * V^2
    
    # The numbers present in this final equation are the coefficient '2'
    # and the exponent '2'.
    
    coefficient = 2
    exponent = 2
    
    print("The total Joule heat (H) dissipated into the atmosphere is equal to the initial electrostatic energy of the sphere.")
    print("The derived formula is:")
    
    # We print the final equation, explicitly showing the numbers involved.
    print(f"H = {coefficient} * π * ε₀ * a * V^{exponent}")
    
    print("\nWhere:")
    print("  H = Total Joule heat dissipated")
    print("  π (pi) = The mathematical constant (~3.14159)")
    print("  ε₀ (epsilon_nought) = The permittivity of free space (a physical constant)")
    print("  a = The initial radius of the sphere")
    print("  V = The initial potential of the sphere")

if __name__ == "__main__":
    solve_joule_heat_dissipation()