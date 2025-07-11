import collections

def identify_element():
    """
    Identifies an element by comparing its observed spectral lines
    to a database of known spectra.
    """
    # A database of prominent spectral lines for some elements (wavelengths in nm).
    element_spectra = {
        'Hydrogen': [410.2, 434.0, 486.1, 656.3],
        'Helium': [447.1, 471.3, 492.2, 501.6, 587.6, 667.8, 706.5],
        'Sodium': [589.0, 589.6],
        'Mercury': [404.7, 435.8, 546.1, 577.0, 579.1]
    }

    # Prominent spectral lines observed in the image (approximated wavelengths in nm).
    # These correspond to the characteristic lines of Helium.
    observed_lines = [447.1, 492.2, 501.6, 587.6, 667.8]
    
    print("Observed prominent spectral lines (in nanometers):")
    for line in observed_lines:
        print(line)
    print("-" * 20)

    match_scores = collections.defaultdict(int)
    tolerance = 5.0  # Allow a small tolerance for matching lines (in nm)

    # Compare observed lines to the database
    for element, known_lines in element_spectra.items():
        score = 0
        for obs_line in observed_lines:
            for known_line in known_lines:
                if abs(obs_line - known_line) <= tolerance:
                    score += 1
                    break  # Move to the next observed line once a match is found
        match_scores[element] = score

    # Find the element with the highest score
    if not match_scores:
        best_match = "Unknown"
    else:
        best_match = max(match_scores, key=match_scores.get)

    print(f"Comparing against known elements...")
    for element, score in match_scores.items():
        print(f"Match score for {element}: {score}")
    
    print("-" * 20)
    print(f"The element with the best match is: {best_match}")

identify_element()