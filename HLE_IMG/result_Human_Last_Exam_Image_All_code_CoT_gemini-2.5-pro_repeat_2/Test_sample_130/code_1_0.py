def identify_element_from_spectrum():
    """
    Identifies the element by its characteristic spectral lines and prints the evidence.
    The visual pattern in the image, especially the prominent green line and yellow doublet,
    is a classic fingerprint of Mercury (Hg).
    """

    # A dictionary of prominent visible spectral lines for Mercury (Hg)
    # Wavelengths are given in nanometers (nm).
    mercury_spectral_lines = {
        "Violet": 404.7,
        "Blue": 435.8,
        "Green": 546.1,
        "Yellow 1": 577.0,
        "Yellow 2": 579.1,
        "Red": 690.7
    }

    print("The element with the given spectral lines is Mercury (Hg).")
    print("\nThis identification is based on its unique and well-known emission spectrum.")
    print("Here are the wavelengths of its most prominent visible lines:")

    for line_name, wavelength in mercury_spectral_lines.items():
        # The instruction was "output each number in the final equation".
        # While there's no equation, we can format the output to highlight each number.
        print(f"-> Prominent {line_name} line: {wavelength} nm")
        
    print("\nThe bright green line at 546.1 nm and the yellow doublet at 577.0 nm and 579.1 nm are particularly strong identifiers.")

if __name__ == "__main__":
    identify_element_from_spectrum()