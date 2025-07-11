def calculate_hydrogen_spectrum():
    """
    Calculates the wavelengths of the first four lines of the Hydrogen Balmer series
    using the Rydberg formula and identifies the element.
    """
    # Rydberg constant in m^-1
    R = 1.097373e7
    # Final energy level for Balmer series
    n1 = 2
    # Initial energy levels for the first four lines
    initial_levels = [3, 4, 5, 6]
    line_names = {
        3: "H-alpha (Red)",
        4: "H-beta (Cyan)",
        5: "H-gamma (Blue-Violet)",
        6: "H-delta (Violet)"
    }

    print("The spectral lines in the image correspond to the Balmer series of Hydrogen.")
    print("We can verify this using the Rydberg formula: 1/wavelength = R * (1/n1^2 - 1/n2^2)\n")

    for n2 in initial_levels:
        # Calculate 1/wavelength
        inv_lambda = R * ((1 / (n1**2)) - (1 / (n2**2)))
        
        # Calculate wavelength in meters
        lambda_m = 1 / inv_lambda
        
        # Convert wavelength to nanometers
        lambda_nm = lambda_m * 1e9

        print(f"Calculating the {line_names[n2]} line (transition from n={n2} to n={n1}):")
        
        # Displaying the calculation step-by-step
        term1 = n1**2
        term2 = n2**2
        print(f"  1/wavelength = {R:.4e} * (1/{term1} - 1/{term2})")
        
        fraction_val = (1/term1) - (1/term2)
        print(f"  1/wavelength = {R:.4e} * ({fraction_val:.6f})")
        
        print(f"  1/wavelength = {inv_lambda:.4e} m^-1")
        print(f"  wavelength = {lambda_m:.4e} m = {lambda_nm:.2f} nm\n")

    print("The calculated wavelengths match the prominent red, cyan, and blue-violet lines")
    print("in the provided spectrum. Therefore, the element is Hydrogen.")

calculate_hydrogen_spectrum()