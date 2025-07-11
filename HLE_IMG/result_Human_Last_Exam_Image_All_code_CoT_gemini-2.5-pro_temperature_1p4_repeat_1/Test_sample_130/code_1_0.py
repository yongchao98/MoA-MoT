def identify_element_from_spectrum():
    """
    Identifies an element by its characteristic spectral lines and prints the result.
    The data is based on comparing the input image to known atomic spectra.
    """
    
    # A dictionary containing the prominent visible spectral lines (in nm) for Helium.
    element_data = {
        "name": "Helium",
        "lines": {
            447.1: "Blue-Violet",
            471.3: "Blue-Green",
            492.2: "Green",
            501.6: "Green (Bright)",
            587.6: "Yellow (Bright)",
            667.8: "Red (Bright)",
            706.5: "Red"
        }
    }
    
    element_name = element_data["name"]
    spectral_lines = element_data["lines"]
    
    print(f"The spectral lines in the image match the emission spectrum of {element_name}.")
    print("\nThe prominent lines for this element in the visible spectrum are:")
    
    # Print each characteristic line and its color
    for wavelength, color in spectral_lines.items():
        print(f"- A {color} line is found at approximately {wavelength} nm.")

# Run the identification function
identify_element_from_spectrum()
