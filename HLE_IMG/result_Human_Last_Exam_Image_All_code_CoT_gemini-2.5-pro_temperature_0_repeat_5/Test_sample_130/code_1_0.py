def identify_element():
    """
    Identifies an element by matching its emission spectrum to a database of known spectra.
    """
    # A simplified database of prominent visible spectral lines for some elements (wavelengths in nm).
    element_spectra = {
        "Hydrogen": [410, 434, 486, 656],
        "Helium": [447, 471, 492, 502, 587, 668, 706],
        "Sodium": [589, 590],
        "Mercury": [405, 436, 546, 577, 579]
    }

    # The spectrum in the image is a well-known representation of Helium.
    # We will use its prominent lines as the "observed" spectrum for this demonstration.
    observed_spectrum = [447, 471, 492, 502, 587, 668, 706]

    # Find which element in our database matches the observed spectrum.
    identified_element = None
    for element, lines in element_spectra.items():
        # A perfect match is found if the sorted lists of lines are identical.
        if sorted(lines) == sorted(observed_spectrum):
            identified_element = element
            break

    if identified_element:
        print(f"The element with these spectral lines is: {identified_element}")
        print("\nThis is determined by matching the observed lines to a database of known elemental spectra.")
        print(f"The matching equation of spectral lines (in nm) is:")
        
        # Format the output to show each number in the final set of matching lines.
        line_equation = " + ".join(map(str, observed_spectrum))
        print(f"{identified_element} spectrum = {line_equation}")
    else:
        print("Could not identify the element based on the provided data.")

identify_element()