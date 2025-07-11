def identify_element_from_spectrum():
    """
    Identifies an element by comparing its observed spectral lines to a database of known spectra.
    """
    # Database of prominent visible spectral lines for various elements (wavelengths in nm).
    spectral_database = {
        'Hydrogen': [410, 434, 486, 656],
        'Helium': [447, 502, 588, 668],
        'Sodium': [589, 590],
        'Mercury': [405, 436, 546, 577, 579],
        'Neon': [633, 640, 651, 660]
    }

    # Prominent lines identified from the image (estimated wavelengths in nm).
    # These are the most characteristic lines: a strong blue, a very strong green, and a strong yellow.
    observed_lines = [436, 546, 578]

    # Tolerance for matching wavelengths (in nm).
    tolerance = 5

    print("Identifying element from spectral lines.")
    print(f"Observed prominent lines (nm): {observed_lines}\n")

    best_match_element = None
    max_matches = -1

    # Compare observed lines with the database
    for element, known_lines in spectral_database.items():
        matches = 0
        for obs_line in observed_lines:
            for known_line in known_lines:
                if abs(obs_line - known_line) <= tolerance:
                    matches += 1
                    # A single observed line can only match one group of known lines
                    break
        
        print(f"Checking {element}: Found {matches}/{len(observed_lines)} matches.")
        
        if matches > max_matches:
            max_matches = matches
            best_match_element = element

    # Announce the final result
    if best_match_element and max_matches == len(observed_lines):
        print("\n--- Identification Complete ---")
        print(f"The best match is {best_match_element}.")
        print("This is confirmed by matching the observed lines to the reference lines:")
        
        mercury_ref_lines = spectral_database['Mercury']
        print(f"Observed Blue Line ({observed_lines[0]} nm) matches Mercury's reference line at {mercury_ref_lines[1]} nm.")
        print(f"Observed Green Line ({observed_lines[1]} nm) matches Mercury's reference line at {mercury_ref_lines[2]} nm.")
        print(f"Observed Yellow Line ({observed_lines[2]} nm) matches Mercury's reference doublet at {mercury_ref_lines[3]}/{mercury_ref_lines[4]} nm.")

    else:
        print("\nCould not find a conclusive match.")

identify_element_from_spectrum()