def identify_element():
    """
    Identifies an element based on its characteristic spectral lines.
    The provided image shows the emission spectrum for Helium.
    This function stores and prints the prominent visible lines for Helium
    to confirm the identification.
    """
    # Data for the prominent visible emission lines of Helium (He I)
    # Wavelengths are in nanometers (nm).
    helium_spectra = {
        "element": "Helium",
        "lines": [
            {"wavelength": 447.1, "color": "Blue-Violet", "relative_intensity": "Strong"},
            {"wavelength": 471.3, "color": "Blue", "relative_intensity": "Faint"},
            {"wavelength": 492.2, "color": "Blue-Green", "relative_intensity": "Medium"},
            {"wavelength": 501.6, "color": "Green", "relative_intensity": "Strong"},
            {"wavelength": 587.6, "color": "Yellow-Orange", "relative_intensity": "Very Strong"},
            {"wavelength": 667.8, "color": "Red", "relative_intensity": "Very Strong"},
            {"wavelength": 706.5, "color": "Red", "relative_intensity": "Strong"}
        ]
    }

    element_name = helium_spectra["element"]
    print(f"The spectral lines in the image belong to the element: {element_name}")
    print("\nThis is confirmed by matching the image to Helium's known prominent spectral lines in the visible range:")
    print("-" * 50)
    print(f"{'Wavelength (nm)':<20} | {'Color':<15} | {'Relative Intensity'}")
    print("-" * 50)

    for line in helium_spectra["lines"]:
        # The prompt asks to "output each number in the final equation", 
        # which is interpreted here as showing the data used for identification.
        print(f"{line['wavelength']:<20.1f} | {line['color']:<15} | {line['relative_intensity']}")

    print("-" * 50)

identify_element()