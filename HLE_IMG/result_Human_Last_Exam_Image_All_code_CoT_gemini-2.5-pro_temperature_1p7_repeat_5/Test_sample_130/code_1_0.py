def identify_element_from_spectrum():
    """
    Identifies an element by its characteristic emission spectrum.

    The analysis is based on comparing the visual pattern in the provided image
    to the known spectra of elements. The image shows a line spectrum with prominent
    lines across the visible range (red, orange, yellow, green, blue, violet).
    This specific pattern is the "fingerprint" of Helium.
    """

    element_name = "Helium"

    # Prominent spectral lines of Helium in the visible range (wavelengths in nanometers).
    # The prompt requires outputting the numbers involved in the "equation" or reasoning.
    spectral_lines = {
        "Red": 667.8,
        "Yellow": 587.6,
        "Green": 501.6,
        "Cyan": 492.2,
        "Blue-Green": 471.3,
        "Blue-Violet": 447.1
    }

    print(f"The element with the provided spectral lines is: {element_name}")
    print("\nThis identification is based on the unique pattern of emission lines, which matches the known spectrum for this element.")
    print("The key spectral lines visible in the image correspond to the following approximate wavelengths:")

    for color, wavelength in spectral_lines.items():
        # Printing each number (wavelength) as requested by the prompt.
        print(f"- A prominent {color} line at {wavelength} nm.")

identify_element_from_spectrum()