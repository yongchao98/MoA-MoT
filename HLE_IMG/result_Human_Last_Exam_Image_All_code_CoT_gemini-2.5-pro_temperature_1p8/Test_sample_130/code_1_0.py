import collections

def identify_element():
    """
    Identifies an element by comparing its observed spectral lines
    with a database of known spectral lines.
    """
    # Database of prominent visible spectral lines for some elements (in nanometers)
    spectral_lines_data = {
        'Hydrogen': [410.2, 434.0, 486.1, 656.3],
        'Helium': [447.1, 501.6, 587.6, 667.8],
        'Mercury': [404.7, 435.8, 546.1, 577.0, 579.1],
        'Sodium': [589.0, 589.6] # Strong doublet
    }

    # Prominent spectral lines estimated from the provided image (in nanometers)
    # The image shows strong lines in violet, blue, green, and a yellow-orange doublet.
    observed_spectrum_approx_nm = [405, 436, 546, 578]

    print("Attempting to identify the element from its spectral lines.")
    print(f"Observed prominent lines (approx. nm): {observed_spectrum_approx_nm}\n")

    best_match_element = None
    best_match_score = -1
    match_details = collections.defaultdict(list)

    # Tolerance for matching wavelengths (in nm)
    tolerance = 8.0

    # Compare the observed spectrum with the known spectra
    for element, known_lines in spectral_lines_data.items():
        score = 0
        current_matches = []
        for obs_line in observed_spectrum_approx_nm:
            for known_line in known_lines:
                if abs(obs_line - known_line) <= tolerance:
                    score += 1
                    current_matches.append(f"{obs_line}nm ~ {known_line}nm")
                    break # Count each observed line only once per element
        
        match_details[element] = current_matches
        
        if score > best_match_score:
            best_match_score = score
            best_match_element = element

    print("--- Matching Results ---")
    for element in spectral_lines_data:
        score = len(match_details[element])
        print(f"Element: {element}, Match Score: {score}/{len(observed_spectrum_approx_nm)}")
        if match_details[element]:
            print(f"  Matched lines: {', '.join(match_details[element])}")
    print("------------------------\n")
    
    if best_match_element:
        print(f"The element with the best match is {best_match_element}.")
    else:
        print("Could not identify the element.")

identify_element()