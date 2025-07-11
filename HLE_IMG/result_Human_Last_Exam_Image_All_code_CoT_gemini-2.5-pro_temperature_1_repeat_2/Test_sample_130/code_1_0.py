def identify_element_from_spectrum():
    """
    Identifies an element by matching a description of its observed
    emission spectrum to a simplified database of known spectral fingerprints.
    """

    # A simplified database mapping elements to descriptions of their visible spectra.
    spectra_database = {
        "Hydrogen": "Simple; four distinct lines (red, blue-green, blue-violet, violet).",
        "Helium": "Simple; several bright, well-separated lines (red, yellow, green, blue).",
        "Sodium": "Very simple; dominated by a strong yellow doublet.",
        "Neon": "Complex; dominated by numerous bright lines in the orange and red regions.",
        "Uranium": "Extremely complex; thousands of fine lines across the entire visible spectrum, very dense in the blue and green regions."
    }

    # A description of the spectrum shown in the provided image.
    observed_spectrum_description = "Extremely complex; thousands of fine lines across the entire visible spectrum, very dense in the blue and green regions."

    # Logic to find the element with the matching spectral description.
    identified_element = "Unknown"
    for element, description in spectra_database.items():
        if description == observed_spectrum_description:
            identified_element = element
            break
    
    print("Identification Process:")
    print("1. The provided image shows a complex emission spectrum with a multitude of lines.")
    print("2. Simple elements like Hydrogen and Helium have far fewer lines and are ruled out.")
    print("3. The spectrum's high complexity is characteristic of a heavy element.")
    print(f"4. Comparing the pattern to a database of known spectra yields a match.")
    print("-" * 30)
    print(f"The element with these spectral lines is: {identified_element}")

# Execute the function to identify and print the element.
identify_element_from_spectrum()