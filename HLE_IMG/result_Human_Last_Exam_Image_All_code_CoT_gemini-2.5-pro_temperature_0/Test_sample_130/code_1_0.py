def identify_element_from_spectrum():
    """
    Identifies an element by its characteristic spectral lines.
    The provided spectrum image matches the known emission spectrum of Helium.
    This function prints the prominent spectral lines of Helium as confirmation.
    """

    # Prominent spectral lines of Helium in the visible spectrum (wavelengths in nanometers)
    helium_spectral_lines = {
        "Red": [706.5, 667.8],
        "Yellow": [587.6],
        "Green": [504.8, 501.6, 492.2],
        "Blue": [471.3],
        "Violet": [447.1, 438.8, 412.1, 402.6]
    }

    print("The element with the given spectral lines is Helium (He).")
    print("\nHere are its prominent characteristic spectral lines (wavelengths in nm):")

    for color, wavelengths in helium_spectral_lines.items():
        for wavelength in wavelengths:
            # The prompt asks to output each number in the final equation.
            # We will format this as "Color: Wavelength nm"
            print(f"{color}: {wavelength} nm")

identify_element_from_spectrum()