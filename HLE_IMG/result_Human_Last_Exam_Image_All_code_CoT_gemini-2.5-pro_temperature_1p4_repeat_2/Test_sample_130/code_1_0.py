def identify_element_from_spectrum():
    """
    Identifies an element by comparing the visual characteristics of its
    emission spectrum with known data. The provided spectrum is very complex
    and dense with lines, a key characteristic of certain elements.
    """
    # A small subset of prominent visible spectral lines for comparison (wavelengths in nm).
    # Note: Iron has thousands of lines; this is just a sample.
    spectral_data = {
        "Hydrogen": [410, 434, 486, 656],
        "Helium": [447, 492, 502, 588, 668],
        "Iron": [438.3, 440.4, 495.7, 527.0, 537.1, 613.7, 639.3, 640.0]
    }

    print("Step 1: Analyze the visual properties of the given spectrum.")
    print("The spectrum image shows an extremely high number of emission lines across the entire visible range (blue, green, red).")
    print("-" * 30)

    print("Step 2: Compare these properties with known elemental spectra.")
    print(f"\n- A simple spectrum like Hydrogen's (lines at {spectral_data['Hydrogen']} nm) is not a match.")
    print(f"- A moderately complex spectrum like Helium's (lines at {spectral_data['Helium']} nm) is also not dense enough.")
    print("\n- An element known for a very dense and complex spectrum is Iron (Fe).")
    print("-" * 30)

    print("Step 3: Conclude based on the evidence.")
    print("The high density of lines is a classic fingerprint of Iron.")
    print("Here is a small sample of the thousands of spectral lines for Iron that create this pattern:")
    
    fe_lines = spectral_data['Iron']
    print(f"Sample Iron Lines (nm): {fe_lines[0]}, {fe_lines[1]}, {fe_lines[2]}, {fe_lines[3]}, {fe_lines[4]}, {fe_lines[5]}, {fe_lines[6]}, {fe_lines[7]}")

    print("\nTherefore, the element with these spectral lines is Iron.")

identify_element_from_spectrum()