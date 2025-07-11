def identify_element_from_spectrum():
    """
    Identifies an element by comparing its known spectral lines to an observed spectrum.
    The observed spectrum has a characteristic bright yellow doublet.
    """

    # Database of prominent visible spectral lines for some elements (wavelengths in nm)
    spectral_data = {
        "Hydrogen": [410.2, 434.0, 486.1, 656.3],
        "Helium": [447.1, 492.2, 501.6, 587.6, 667.8],
        "Sodium": [497.9, 514.9, 568.8, 589.0, 589.6, 615.4],
        "Mercury": [404.7, 435.8, 546.1, 577.0, 579.1]
    }

    print("Analyzing the provided emission spectrum...")
    print("The most distinctive feature is a very strong doublet (two lines very close together) in the yellow region of the spectrum.")
    print("\nLet's compare this observation with the known spectral lines of a few elements:\n")

    # The characteristic feature we are looking for is the yellow doublet around 589 nm.
    target_wavelength_1 = 589.0
    target_wavelength_2 = 589.6
    identified_element = None

    for element, lines in spectral_data.items():
        print(f"- {element}: {lines}")
        # Check if the element has the characteristic sodium D-lines
        if target_wavelength_1 in lines and target_wavelength_2 in lines:
            identified_element = element

    print("\nConclusion:")
    if identified_element:
        print(f"The spectrum with a strong yellow doublet at {target_wavelength_1} nm and {target_wavelength_2} nm is characteristic of {identified_element}.")
    else:
        print("Could not identify the element from the available data.")

if __name__ == "__main__":
    identify_element_from_spectrum()
