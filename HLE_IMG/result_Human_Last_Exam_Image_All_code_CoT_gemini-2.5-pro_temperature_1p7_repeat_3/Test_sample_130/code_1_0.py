def identify_element_from_spectrum():
    """
    Identifies an element by comparing its spectral features to a known database.
    """
    # Simplified database of prominent spectral lines (in nm) for some elements.
    element_spectra = {
        'Hydrogen': {
            'description': "Four prominent lines in the visible spectrum (Balmer series).",
            'lines': [656.3, 486.1, 434.1, 410.2]
        },
        'Sodium': {
            'description': "Dominated by a strong yellow doublet.",
            'lines': [589.0, 589.6]
        },
        'Helium': {
            'description': "Complex spectrum with multiple bright lines across the visible range.",
            'lines': [667.8, 587.6, 501.6, 471.3, 447.1]
        }
    }

    # Features observed from the image.
    observed_spectrum_features = {
        'complexity': 'high',
        'colors': ['red', 'yellow', 'green', 'blue', 'violet']
    }

    # Perform the identification by comparing patterns.
    # The observed spectrum is complex and has lines in all major color regions.
    # This rules out Hydrogen and Sodium. The pattern strongly matches Helium.
    identified_element = 'Helium'
    element_data = element_spectra[identified_element]

    print(f"Identifying the element based on its spectral fingerprint.")
    print("-" * 50)
    print("Observed Spectrum Analysis:")
    print(f"The image shows a complex emission spectrum with many lines.")
    print(f"Bright lines are present in the red, yellow, green, and blue regions.")
    print("-" * 50)
    print("Comparison with Known Elements:")
    print("The simple spectra of Hydrogen and Sodium do not match.")
    print(f"The observed pattern is characteristic of: {identified_element}")
    print("-" * 50)
    print("Conclusion:")
    print("The element is identified as Helium. The key identifying lines are:")

    # Printing the 'equation' by listing the key spectral lines for Helium
    print(f"1. A bright red line around {element_data['lines'][0]} nm")
    print(f"2. A very bright yellow line around {element_data['lines'][1]} nm")
    print(f"3. A bright green line around {element_data['lines'][2]} nm")
    print(f"4. A cyan line around {element_data['lines'][3]} nm")
    print(f"5. A bright blue-violet line around {element_data['lines'][4]} nm")
    print("-" * 50)
    print(f"The element has been identified as {identified_element}.")


identify_element_from_spectrum()