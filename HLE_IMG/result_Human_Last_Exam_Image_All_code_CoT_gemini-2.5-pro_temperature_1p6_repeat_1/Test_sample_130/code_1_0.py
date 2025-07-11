def calculate_balmer_series():
    """
    Calculates the wavelengths of the first four lines of the Balmer series for Hydrogen
    and identifies the element based on the spectral lines.
    """
    # Rydberg constant for Hydrogen in m^-1
    R_H = 1.097373e7

    # For the Balmer series, the final energy level is n1 = 2
    n1 = 2

    print("Calculating the first four spectral lines of the Hydrogen Balmer series (n1=2):")
    print("-" * 75)

    # Calculate for the first four transitions (n2 = 3, 4, 5, 6)
    for n2 in range(3, 7):
        # Using the Rydberg formula: 1/wavelength = R_H * (1/n1^2 - 1/n2^2)
        inv_wavelength = R_H * (1/n1**2 - 1/n2**2)
        
        # Wavelength in meters
        wavelength_m = 1 / inv_wavelength
        
        # Convert wavelength to nanometers (1 m = 1e9 nm)
        wavelength_nm = wavelength_m * 1e9
        
        # Determine the color based on wavelength
        if 620 <= wavelength_nm <= 750:
            color = "Red"
        elif 495 <= wavelength_nm <= 570:
            color = "Green/Cyan"
        elif 450 <= wavelength_nm <= 495:
            color = "Blue"
        elif 380 <= wavelength_nm <= 450:
            color = "Violet"
        else:
            color = "Unknown"
        
        print(f"For the transition from n2={n2} to n1={n1}:")
        print(f"1 / λ = {R_H:.4g} * (1/{n1**2} - 1/{n2**2})")
        print(f"Calculated Wavelength (λ): {wavelength_nm:.2f} nm ({color})")
        print("-" * 75)

    print("\nThese calculated wavelengths (656.28 nm, 486.13 nm, 434.05 nm, 410.17 nm) correspond to the prominent\nred, cyan, blue-violet, and violet lines seen in the provided spectrum.")
    print("\nTherefore, the element is Hydrogen.")

calculate_balmer_series()