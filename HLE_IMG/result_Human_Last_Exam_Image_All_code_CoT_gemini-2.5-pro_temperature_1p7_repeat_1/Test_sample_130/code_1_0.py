def identify_element_from_spectrum():
    """
    Identifies an element by comparing the characteristics of its
    emission spectrum to a database of known spectra.
    """

    # A simplified database of spectral characteristics for common elements.
    spectral_database = {
        "Hydrogen": "A simple spectrum with four prominent lines in the visible range (red, blue-green, blue, violet).",
        "Helium": "A spectrum with several bright, distinct lines spread out across the visible range.",
        "Sodium": "A spectrum dominated by a very strong and close pair of lines in the yellow region.",
        "Neon": "A spectrum with numerous, very bright lines in the red-orange part of the spectrum.",
        "Iron": "A very complex spectrum with thousands of lines distributed across the entire visible range."
    }

    # Description of the observed spectrum from the image.
    observed_spectrum_description = "A very complex spectrum with thousands of lines distributed across the entire visible range."

    # Find the element that matches the description.
    identified_element = None
    for element, description in spectral_database.items():
        if description == observed_spectrum_description:
            identified_element = element
            break

    if identified_element:
        print(f"The observed spectrum is characteristic of the element: {identified_element}")
    else:
        print("Could not identify the element based on the provided description.")

identify_element_from_spectrum()