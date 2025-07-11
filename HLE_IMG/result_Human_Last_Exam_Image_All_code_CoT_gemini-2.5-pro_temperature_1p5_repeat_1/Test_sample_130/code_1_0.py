def identify_element_from_spectrum():
    """
    Identifies the element by calculating and comparing its spectral lines
    with the pattern shown in the image.
    """

    print("Analyzing the provided emission spectrum...")
    print("The spectrum shows a characteristic pattern of lines known as the Balmer series.")
    print("This series is unique to the hydrogen atom.")
    print("\nWe can verify this using the Rydberg formula for hydrogen:")
    print("1 / λ = R * (1/n₁² - 1/n₂²)\n")

    # Constants
    R = 1.097373e7  # Rydberg constant for Hydrogen (in m⁻¹)
    n1 = 2          # Final energy level for the Balmer series

    # --- Calculation for the H-alpha line (the prominent red line) ---
    n2 = 3          # Initial energy level for the H-alpha line
    
    # Calculate the inverse of the wavelength
    inv_lambda = R * (1/n1**2 - 1/n2**2)
    
    # Calculate the wavelength in meters, then convert to nanometers
    lambda_m = 1 / inv_lambda
    lambda_nm = lambda_m * 1e9
    
    print("Let's calculate the wavelength of the most prominent red line (H-alpha):")
    print(f"Here, R = {R:.6e} m⁻¹, n₁ = {n1}, and n₂ = {n2}.")
    print(f"1 / λ = {R:.6e} * (1/{n1}² - 1/{n2}²)")
    print(f"1 / λ = {R:.6e} * (1/{n1**2} - 1/{n2**2})")
    print(f"1 / λ = {R:.6e} * ({1/n1**2 - 1/n2**2})")
    print(f"1 / λ = {inv_lambda:.6e} m⁻¹")
    print(f"λ = 1 / {inv_lambda:.6e} = {lambda_m:.6e} m")
    print(f"λ = {lambda_nm:.1f} nm\n")
    
    print("This calculated value corresponds to the bright red line in the spectrum.")
    
    # Known Balmer series lines for comparison
    balmer_series = {
        'H-alpha (n=3)': (656.3, 'Red'),
        'H-beta (n=4)': (486.1, 'Cyan (Blue-Green)'),
        'H-gamma (n=5)': (434.1, 'Blue-Violet'),
        'H-delta (n=6)': (410.2, 'Violet')
    }

    print("The primary visible lines for Hydrogen (Balmer series) are:")
    for line, (wavelength, color) in balmer_series.items():
        print(f"- {line}: ~{wavelength} nm ({color})")

    print("\nThe pattern of lines in the image perfectly matches this series.")
    print("\nConclusion: The element is Hydrogen.")

identify_element_from_spectrum()