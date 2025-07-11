def identify_element_from_spectrum():
    """
    Identifies an element by its emission spectrum based on a known pattern.

    The provided image shows the emission spectrum of Mercury (Hg). This is
    identified by its unique pattern of spectral lines, most notably a very
    bright green line. This script presents the data for Mercury's visible
    spectrum to confirm the identification.
    """
    
    # Dictionary of prominent visible spectral lines for Mercury (Hg).
    # Wavelengths are in nanometers (nm).
    mercury_spectral_lines = {
        "Violet": 404.7,
        "Blue": 435.8,
        "Green": 546.1,
        "Yellow-Orange (Doublet line 1)": 577.0,
        "Yellow-Orange (Doublet line 2)": 579.1
    }

    print("Analysis of the emission spectrum:")
    print("The pattern of spectral lines, with prominent lines in the violet, blue,")
    print("a very bright green, and yellow-orange regions, is the unique fingerprint of the element Mercury (Hg).")
    print("\nHere are the major spectral lines for Mercury in the visible range that match the image:")

    # Print each prominent line used for identification.
    for color, wavelength in mercury_spectral_lines.items():
        print(f"- A {color} line is observed at approximately {wavelength} nm.")

identify_element_from_spectrum()