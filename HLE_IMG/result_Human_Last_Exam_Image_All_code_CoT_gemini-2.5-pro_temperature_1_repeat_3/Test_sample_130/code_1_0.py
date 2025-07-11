def identify_element():
    """
    Identifies the element based on its characteristic spectral lines visible in the image.
    The pattern strongly matches the emission spectrum of Mercury (Hg).
    """

    element = "Mercury (Hg)"

    # Prominent spectral lines of Mercury in the visible spectrum (in nanometers)
    spectral_lines = {
        404.7: "Violet",
        435.8: "Blue",
        546.1: "Green",
        577.0: "Yellow",
        579.1: "Yellow-Orange"
    }

    print(f"The element with the spectral lines shown in the image is: {element}")
    print("\nThis identification is based on the following characteristic emission lines which are visible in the spectrum:")

    # Print each prominent line as evidence, matching it to the image
    print("\nImage Spectrum => Mercury's Known Lines")
    print("-" * 40)
    for wavelength, color in spectral_lines.items():
        print(f"A prominent {color} line is visible => Corresponds to the {wavelength} nm line of Mercury.")

identify_element()