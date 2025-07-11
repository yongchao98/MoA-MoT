def identify_element_from_spectrum():
    """
    Identifies an element by comparing its spectral features to known elements.
    """
    
    # Known prominent visible spectral lines (wavelengths in nanometers)
    hydrogen_balmer_lines = {
        'red': 656.3,
        'cyan': 486.1,
        'blue-violet': 434.1,
        'violet': 410.2
    }
    
    helium_lines = {
        'red': 667.8,
        'yellow': 587.6,
        'green': 501.6,
        'cyan': 492.2,
        'blue-violet': 447.1
    }

    # Step 1: Describe the features observed in the provided spectrum image.
    print("Step 1: Analyzing the observed spectrum from the image.")
    print("The spectrum shows multiple bright lines including:")
    print("- Prominent red lines")
    print("- A prominent yellow-orange line")
    print("- Multiple green and cyan lines")
    print("- Several blue-violet lines")
    print("-" * 30)

    # Step 2: Compare with Hydrogen
    print("Step 2: Comparing with the spectrum of Hydrogen.")
    print(f"Hydrogen's visible (Balmer) spectrum has four main lines at approximately:")
    print(f"Red: {hydrogen_balmer_lines['red']} nm")
    print(f"Cyan: {hydrogen_balmer_lines['cyan']} nm")
    print(f"Blue-Violet: {hydrogen_balmer_lines['blue-violet']} nm")
    print(f"Violet: {hydrogen_balmer_lines['violet']} nm")
    print("\nConclusion for Hydrogen: The observed spectrum is much more complex and contains a bright yellow line, which is absent in Hydrogen's visible spectrum. Therefore, the element is not Hydrogen.")
    print("-" * 30)

    # Step 3: Compare with Helium
    print("Step 3: Comparing with the spectrum of Helium.")
    print("Helium's visible spectrum has several prominent lines, including:")
    print(f"Red: {helium_lines['red']} nm")
    print(f"Yellow: {helium_lines['yellow']} nm")
    print(f"Green: {helium_lines['green']} nm")
    print(f"Cyan: {helium_lines['cyan']} nm")
    print(f"Blue-Violet: {helium_lines['blue-violet']} nm")
    print("\nConclusion for Helium: The pattern of lines observed in the image, including the bright red, yellow, green, and blue-violet lines, is a characteristic fingerprint of Helium.")
    print("-" * 30)
    
    # Final Conclusion
    print("Final Conclusion: Based on the comparison, the element that produces these spectral lines is Helium.")

identify_element_from_spectrum()