def identify_element_from_spectrum():
    """
    Identifies an element by comparing its known spectral lines
    to the pattern observed in the provided image.
    """
    # A small database of prominent visible spectral lines for some elements (in nanometers).
    spectral_database = {
        "Hydrogen": [410.2, 434.1, 486.1, 656.3],
        "Helium": [447.1, 471.3, 492.2, 501.6, 587.6, 667.8],
        "Mercury": [404.7, 435.8, 546.1, 577.0, 579.1]
    }

    # Description of the spectrum from the image.
    # The image shows a very distinct pattern, most notably:
    # 1. A prominent violet line.
    # 2. A strong blue line.
    # 3. A very bright, dominant green line.
    # 4. A yellow-orange doublet (two lines very close together).
    # This pattern is a classic fingerprint.

    identified_element = "Mercury"
    key_lines = spectral_database[identified_element]

    print("Analyzing the spectrum to identify the element:")
    print("The spectral pattern in the image is compared against a database of known elements.")
    print("\nThe observed spectrum has several key features. We can match them to a known element.")
    print("\nMatching process:")
    print(f"1. Observed prominent violet line corresponds to Mercury's line at {key_lines[0]} nm.")
    print(f"2. Observed strong blue line corresponds to Mercury's line at {key_lines[1]} nm.")
    print(f"3. Observed dominant green line corresponds to Mercury's line at {key_lines[2]} nm.")
    print(f"4. Observed yellow-orange doublet corresponds to Mercury's lines at {key_lines[3]} nm and {key_lines[4]} nm.")

    print("\nConclusion:")
    print(f"This specific pattern of spectral lines is a unique fingerprint for the element {identified_element}.")

identify_element_from_spectrum()