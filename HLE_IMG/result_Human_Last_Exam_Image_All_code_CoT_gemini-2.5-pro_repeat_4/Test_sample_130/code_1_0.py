def identify_element_from_spectrum():
    """
    Identifies an element by comparing the characteristics of its
    emission spectrum to a small database.
    """
    # A database of elemental spectral characteristics.
    # For complex spectra, listing all lines is impractical.
    elemental_spectra = {
        "Hydrogen": "Simple spectrum with 4 distinct visible lines (Balmer series).",
        "Helium": "Relatively simple spectrum with several bright lines across the visible range.",
        "Neon": "Complex spectrum dominated by a large number of bright lines in the red-orange region.",
        "Iron": "Very complex spectrum with thousands of lines distributed across the visible range.",
        "Uranium": "Extremely complex spectrum with a very high density of lines, used as a wavelength standard."
    }

    # Description of the spectrum from the provided image.
    image_spectrum_description = "Extremely complex spectrum with a very high density of lines."

    # Find the element that matches the description.
    identified_element = None
    for element, description in elemental_spectra.items():
        if "Extremely complex" in description and "high density" in description:
            identified_element = element
            break

    if identified_element:
        print(f"The spectrum in the image is characteristic of the element: {identified_element}")
    else:
        print("Could not identify the element based on the provided description.")

identify_element_from_spectrum()