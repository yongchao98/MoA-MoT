def identify_element_from_spectrum():
    """
    Identifies the element based on its characteristic spectral lines.
    The image provided corresponds to the emission spectrum of Helium.
    This function will print the element's name and its major visible spectral lines.
    """
    element_name = "Helium"
    
    # Prominent spectral lines of Helium in the visible spectrum (in nanometers)
    # These values correspond to the colored lines seen in the image.
    helium_lines = {
        "Red": 667.8,
        "Yellow": 587.6,
        "Green": 501.6,
        "Cyan": 492.2,
        "Blue-Violet": 447.1,
        "Violet": 402.6
    }

    print(f"The element with the spectral lines shown in the image is: {element_name}")
    print("\nIts prominent spectral lines in the visible range are:")
    
    # The prompt asks to "output each number in the final equation".
    # Since there is no equation, we will list each characteristic wavelength.
    for color, wavelength in helium_lines.items():
        print(f"- A {color} line at approximately {wavelength} nm")

identify_element_from_spectrum()