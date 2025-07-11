def identify_element():
    """
    This function identifies the element based on its known spectral lines
    that match the provided image.
    The image displays the emission spectrum of Mercury (Hg).
    """

    # The element identified based on the spectral pattern
    element_name = "Mercury"
    element_symbol = "Hg"

    # Known prominent visible spectral lines for Mercury (wavelengths in nanometers)
    spectral_lines = {
        "Violet": 404.7,
        "Blue-Violet": 435.8,
        "Green": 546.1,
        "Yellow 1": 577.0,
        "Yellow 2": 579.1
    }

    print(f"The spectrum shown belongs to the element: {element_name} ({element_symbol})")
    print("\nThis is identified by its unique pattern of spectral lines. The prominent visible lines for Mercury are:")

    # Print each prominent line's color and wavelength
    for color, wavelength in spectral_lines.items():
        print(f"- {color}: {wavelength} nm")

identify_element()