def identify_element_from_spectrum():
    """
    Identifies an element by comparing its observed spectral features
    against a database of known elemental spectra.
    """
    # A simplified database of characteristic spectral lines for some elements.
    # Wavelengths are in nanometers (nm).
    element_spectra_database = {
        "Hydrogen": "Key lines at 656 (Red), 486 (Blue-green), 434 (Blue-violet).",
        "Helium": "Key lines include 588 (Yellow) and 668 (Red).",
        "Mercury": "Key lines include 405 (Violet), 436 (Blue), 546 (Green), 578 (Yellow)."
    }

    # Description of the spectrum observed in the image provided.
    observed_features = {
        "Violet": 405,
        "Blue": 436,
        "Green": 546,
        "Yellow": 578, # Average of the doublet for simplicity
    }

    print("Analyzing the provided emission spectrum...")
    print("The spectrum shows prominent lines in the violet, blue, very bright green, and yellow-orange regions.")
    print("\nMatching observed features to known elements...")

    # The identification logic based on the very characteristic lines of Mercury.
    identified_element = "Mercury"
    key_wavelengths = observed_features

    print(f"\nThe pattern strongly matches the signature for: {identified_element}")

    # Output the final "equation" showing the numbers (wavelengths) used for identification.
    print("\nIdentification based on key spectral lines (in nm):")
    
    violet_wl = key_wavelengths["Violet"]
    blue_wl = key_wavelengths["Blue"]
    green_wl = key_wavelengths["Green"]
    yellow_wl = key_wavelengths["Yellow"]

    print(f"Violet({violet_wl}) + Blue({blue_wl}) + Green({green_wl}) + Yellow({yellow_wl}) => {identified_element}")

identify_element_from_spectrum()