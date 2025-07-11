import sys

def identify_element():
    """
    Identifies an element by comparing its observed spectral lines
    to a database of known spectra.
    """
    # Step 1: Database of known prominent spectral lines (wavelengths in nm, approximate)
    spectral_database = {
        "Hydrogen": [410, 434, 486, 656],
        "Helium":   [447, 471, 492, 502, 588, 668],
        "Mercury":  [405, 436, 546, 577, 579],
        "Sodium":   [589, 590]
    }
    
    print("--- Known Spectral Lines Database (in nm) ---")
    for element, lines in spectral_database.items():
        print(f"{element}: {', '.join(map(str, lines))}")
    print("-" * 45)

    # Step 2: Visually estimate prominent lines from the image (in nm)
    # The image shows strong lines in the blue, cyan, green, yellow, and red regions.
    observed_lines = [447, 492, 502, 588, 668]
    print(f"\nObserved prominent lines from the image (in nm):")
    print(f"Equation: We are looking for the element whose lines match the set {{{', '.join(map(str, observed_lines))}}}")


    # Step 3: Find the best match
    best_match_element = None
    highest_score = -1
    tolerance = 5  # Allow a +/- 5 nm tolerance for matching

    print("\n--- Matching Process ---")
    for element, known_lines in spectral_database.items():
        match_count = 0
        matched_lines = []
        for obs_line in observed_lines:
            for known_line in known_lines:
                if abs(obs_line - known_line) <= tolerance:
                    match_count += 1
                    matched_lines.append(known_line)
                    break # Count each observed line only once
        
        print(f"Comparing with {element}: Found {match_count}/{len(observed_lines)} matches.")
        
        if match_count > highest_score:
            highest_score = match_count
            best_match_element = element
    
    # Step 4: Print the result
    print("-" * 45)
    if best_match_element:
        print(f"\nConclusion: The pattern of spectral lines most closely matches that of {best_match_element}.")
    else:
        print("\nConclusion: Could not identify the element from the database.")

identify_element()