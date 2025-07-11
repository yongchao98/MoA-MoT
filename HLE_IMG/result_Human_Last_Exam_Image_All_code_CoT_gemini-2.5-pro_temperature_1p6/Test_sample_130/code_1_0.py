def identify_element_from_spectrum():
    """
    Identifies an element based on the visual characteristics of its emission spectrum.

    The provided spectrum shows a very high density of lines across the visible range.
    This pattern is a classic signature of the element Iron (Fe).
    This program will display some of the prominent spectral lines of Iron.
    """

    # A small, representative sample of the thousands of emission lines for Iron (Fe)
    # Wavelengths are in nanometers (nm).
    iron_spectral_lines = {
        "Violet": [404.58, 406.36, 407.17, 413.21, 414.39],
        "Blue": [438.35, 440.48, 441.51],
        "Green": [526.95, 527.04, 532.81],
        "Yellow-Orange": [577.85, 586.23, 606.55, 613.77],
        "Red": [639.36, 640.00, 667.80]
    }

    print("Analyzing the provided emission spectrum...")
    print("The spectrum displays a very large number of closely spaced lines across the entire visible range.")
    print("This complex, 'picket-fence' pattern is the characteristic signature of Iron (Fe).\n")
    print("Here are some of the prominent emission lines for Iron (wavelengths in nm):")

    for color, wavelengths in iron_spectral_lines.items():
        # Using string formatting to align the output nicely
        print(f"- {color:<14}: ", end="")
        # Join the wavelength numbers into a single string
        line_str = ", ".join(map(str, wavelengths))
        print(line_str)

    print("\nBased on this evidence, the element is Iron.")

identify_element_from_spectrum()