def print_critical_exponent_nu():
    """
    Presents the value of the critical exponent ν for a G₄ (or φ⁴) theory
    in various spatial dimensions, d.
    
    The exponent ν governs the scaling of the correlation length ξ near criticality:
    ξ ~ |T - T_c|⁻ν. In finite-size scaling analysis at the critical temperature,
    the system size L sets the scale, and the correlation length is proportional
    to L.
    """
    print("The precise value of the critical exponent ν for a G₄ (φ⁴) theory is dependent on the spatial dimension d.")
    print("\nBelow are the established values for key dimensions:")
    
    # Case 1: d = 2
    # This is an exact result, famously derived from the solution of the 2D Ising model.
    val_d2 = 1
    print("\n--- For a spatial dimension of d = 2 ---")
    print("This value is known exactly and corresponds to the 2D Ising universality class.")
    # Printing the value as a simple equation.
    print(f"The critical exponent is: ν = {val_d2}")

    # Case 2: d = 3
    # This is the most common non-trivial case, determined with high precision from 
    # Monte Carlo simulations and renormalization group calculations.
    val_d3 = 0.6301
    print("\n--- For a spatial dimension of d = 3 ---")
    print("This corresponds to the 3D Ising universality class, relevant for many physical systems.")
    # Printing the value as a simple equation.
    print(f"The critical exponent is: ν = {val_d3}")

    # Case 3: d >= 4
    # For d at or above the upper critical dimension (d_c = 4), fluctuations are suppressed,
    # and mean-field theory becomes exact.
    val_d4_num = 1
    val_d4_den = 2
    val_d4_dec = 0.5
    print("\n--- For a spatial dimension of d >= 4 ---")
    print("This is the mean-field theory result, valid at and above the upper critical dimension.")
    # Printing the numbers in the final equation.
    print(f"The critical exponent is: ν = {val_d4_num}/{val_d4_den} = {val_d4_dec}")

# Run the function to display the results.
print_critical_exponent_nu()