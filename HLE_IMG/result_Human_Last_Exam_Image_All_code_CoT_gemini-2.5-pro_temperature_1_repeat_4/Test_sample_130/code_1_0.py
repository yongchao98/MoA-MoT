def identify_element():
    """
    Identifies an element by comparing its observed spectral lines
    against a database of known spectra.
    """
    # A database of prominent spectral lines for some elements (wavelengths in nanometers)
    element_spectra = {
        "Hydrogen": [410, 434, 486, 656],
        "Helium": [447, 501, 587, 667],
        "Sodium": [589],  # Strong yellow doublet, often seen as one line
        "Mercury": [404, 436, 546, 578] # Includes the yellow doublet as one average value
    }

    # Prominent lines estimated from the provided image (in nm).
    # We can clearly see a blue-violet, a very strong green, and a yellow-orange line.
    observed_lines = [436, 546, 578]

    print("Starting element identification process...")
    print(f"Observed prominent spectral lines (estimated wavelengths in nm): {observed_lines}")
    print("="*40)

    best_match_element = "Unknown"
    max_matches = 0

    # Iterate through our database to find the best match
    for element, known_lines in element_spectra.items():
        matches = 0
        matched_lines_info = []
        
        # Check each observed line against the element's known lines
        for obs_line in observed_lines:
            # We'll allow a small tolerance (e.g., +/- 5 nm) for our estimation
            for known_line in known_lines:
                if abs(obs_line - known_line) <= 5:
                    matches += 1
                    matched_lines_info.append(f"Observed ~{obs_line}nm -> {element} line at {known_line}nm")
                    break # Found a match for this observed line, move to the next one
        
        print(f"Comparing with {element}:")
        if matched_lines_info:
            for info in matched_lines_info:
                print(f"  - Match found: {info}")
        else:
            print("  - No significant matches found.")
        print(f"Total matches for {element}: {matches}")
        print("-"*40)

        # Update the best match if this element is a better fit
        if matches > max_matches:
            max_matches = matches
            best_match_element = element

    print("\n--- Conclusion ---")
    print(f"The element with the most matches ({max_matches}) is: {best_match_element}")
    print("The distinct pattern of a strong blue-violet line, a very prominent green line, and bright yellow-orange lines is characteristic of Mercury.")

identify_element()